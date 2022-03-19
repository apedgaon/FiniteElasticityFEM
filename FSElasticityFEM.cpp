#include <iostream>
#include <fstream>
#include <iomanip>
#include "FSElasticityFEM.h"

template class FSElasticityFEM<3>;

template<unsigned int dim>
inline void FSElasticityFEM<dim>::generate_mesh()
{
    mesh.generate(geom);
}

template<unsigned int dim>
void FSElasticityFEM<dim>::initialize()
{
    // create map for unconstrained dofs to global dofs
    create_dirichlet_map();

    // intialize global and dirichlet matrix and vectors
    d = Eigen::VectorXd::Zero(gDofs);       // initial zero guess
    R.resize(gDofs);
    Kj.resize(gDofs, gDofs);
    dd.resize(ucDofs);
    Rd.resize(ucDofs);
    Kjd.resize(ucDofs, ucDofs);
    analysis_complete = false;

    // get initial dirichlet displacements
    for (unsigned int idx = 0; idx < ucDofs; ++idx)
        dd(idx) = d(dirmap[idx]);
}

template<unsigned int dim>
void FSElasticityFEM<dim>::assemble()
{
    // Assign zero value to Residual and System Stiffness
    R = Eigen::VectorXd::Zero(gDofs);
    Kj = Eigen::MatrixXd::Zero(gDofs, gDofs);

    // number of local dofs
    unsigned int loc_dofs = mesh.Nen * dim;

    // Initialize local matrices and vectors
    Eigen::VectorXd Rl(loc_dofs);
    Eigen::MatrixXd Kjl(loc_dofs, loc_dofs);
    Eigen::MatrixXd xe(dim, mesh.Nen);
    Eigen::VectorXui gidx(loc_dofs);
    Eigen::VectorXd de(mesh.Nen * dim);
    Eigen::MatrixXd jac(dim, dim), jac_inv(dim, dim), dNdX(mesh.Nen, dim);
    double det_jac;
    Eigen::MatrixXd dfgrd(dim, dim), E_gl(dim, dim), S_pk2(dim, dim);
    Eigen::VectorXd Rq(loc_dofs);
    Eigen::MatrixXd Kjq(loc_dofs, loc_dofs);

    // Element Loop
    for (auto& conn : mesh.elem_conn)
    {
        // Set zero values to residual and system stiffness  
        Rl.setZero();
        Kjl.setZero();

        // Compute elemental coordinates and displacements
        compute_elem_quantities(conn, gidx, xe, de);

        // Quadrature loop
        for (unsigned int iq = 0; iq < mesh.Nq; ++iq)
        {
            // Compute quadrature jacobian and derivatives
            compute_quad_derivatives(iq, xe, de, jac, jac_inv, dNdX, det_jac);

            // Compute quadrature deformation gradient, strain and stress
            compute_quad_deformation(de, dNdX, dfgrd, E_gl, S_pk2);

            // Compute quadrature Residual
            compute_quad_res(iq, dNdX, dfgrd, S_pk2, det_jac, Rq);
            
            // Update elemental Residual 
            Rl += Rq;

            // Compute quadrature NR Stiffness
            compute_quad_stiff(iq, dNdX, dfgrd, S_pk2, det_jac, Kjq);

            // Update elemental NR Stiffness
            Kjl += Kjq;
        }

        // Global Assembly Rl -> R, Kjl -> Kj
        for (unsigned int jdx = 0; jdx < loc_dofs; ++jdx)
        {
            R(gidx(jdx)) += Rl(jdx);
            for (unsigned int idx = 0; idx < loc_dofs; ++idx)
            {
                Kj(gidx(idx), gidx(jdx)) += Kjl(idx, jdx);
            }
        }
    }
}

template<unsigned int dim>
void FSElasticityFEM<dim>::solve()
{
    // Define delta d
    Eigen::VectorXd delta_d(ucDofs);

    // Load step loop
    for (unsigned int istep = 0; istep < sol_ctrls.nsteps; ++istep)
    {
        // Print load step info
        std::cout << "\n";
        std::cout << "Load Step = " << istep + 1 << "\n";
        bool converged = false;

        // Apply load using a multiplier based on load step
        double step_mult = 1.0 * (istep + 1) / sol_ctrls.nsteps;
        apply_step_BC(step_mult);

        // Newton-Raphson (NR) iterations
        for (unsigned int inc = 0; inc < sol_ctrls.ninc_max; ++inc)
        {
            // Assemble Residual and System Stiffness
            assemble();

            // Get Dirchlet Residual and System Stiffness
            apply_dirichlet_BC();

            // Check convergence by computing norm of the residual and comparing with tolerance
            check_convergence(inc, converged);
            if (converged)
                break;

            // If not converged, solve for delta d (currently uses Cholesky)
            delta_d = Kjd.ldlt().solve(-Rd);

            // Update unconstrained displacements
            dd += delta_d;

            // Update full displacements
            reform_full_sol();
        }

        // Message if not converged and return from this method
        if (!converged)
        {
            std::cout << "*** Convergence Failed in max NR iterations!\n"
                << "*** Analysis Incomplete!\n";
            return;
        }
    }

    // Set Analysis complete and print message
    std::cout << "\n*** Analysis Complete!\n";
    analysis_complete = true;
}

template<unsigned int dim>
void FSElasticityFEM<dim>::post_print(std::string file)
{
    if (analysis_complete)
    {
        // raw_sol_print();

        // Write vtk file with mesh and solutions
        fem_to_vtk_vector(file, mesh.elem_conn, mesh.nodal_coords, d);
    }
}

template<unsigned int dim>
void FSElasticityFEM<dim>::compute_elem_quantities(connectivity& conn, Eigen::VectorXui& gidx,
                                                   Eigen::MatrixXd& xe, Eigen::VectorXd& de)
{
    for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            auto loc_idx = dim * ind + idim;
            gidx(loc_idx) = dim * conn[ind] + idim;
            xe(idim, ind) = mesh.nodal_coords[conn[ind]][idim];
            de(loc_idx) = d(gidx(loc_idx));
        }
    }
}

template<unsigned int dim>
void FSElasticityFEM<dim>::compute_quad_derivatives(unsigned int iq, Eigen::MatrixXd& xe, Eigen::VectorXd& de, Eigen::MatrixXd& jac,
                                                    Eigen::MatrixXd& jac_inv, Eigen::MatrixXd& dNdX, double& det_jac)
{
    // Compute isoparametric Jacobian and Inverse
    jac.setZero();
    for (unsigned int jdim = 0; jdim < dim; ++jdim)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
                jac(idim, jdim) += xe(idim, ind) * mesh.dNdxi[jdim](ind, iq);
        }
    }

    det_jac = jac.determinant();
    jac_inv = jac.inverse();

    // Compute dN/dX
    Eigen::MatrixXd dNdxi(mesh.Nen, dim);
    for (unsigned int idim = 0; idim < dim; ++idim)
    {
        for (unsigned int ishp = 0; ishp < mesh.Nen; ++ishp)
            dNdxi(ishp, idim) = mesh.dNdxi[idim](ishp, iq);
    }

    dNdX = dNdxi * jac_inv;
}

template<unsigned int dim>
void FSElasticityFEM<dim>::compute_quad_deformation(Eigen::VectorXd& de, Eigen::MatrixXd& dNdX, Eigen::MatrixXd& dfgrd,
                                                    Eigen::MatrixXd& E_gl, Eigen::MatrixXd& S_pk2)
{
    // Compute du/dX
    Eigen::MatrixXd ugrd = Eigen::MatrixXd::Zero(dim, dim);
    for (unsigned int jdim = 0; jdim < dim; ++jdim)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            for (unsigned int ishp = 0; ishp < mesh.Nen; ++ishp)
                ugrd(idim, jdim) += de(dim * ishp + idim) * dNdX(ishp, jdim);
        }
    }

    // Compute Deformation Gradient, Green Lagrange Strain
    Eigen::MatrixXd del = Eigen::MatrixXd::Identity(dim, dim);
    dfgrd = del + ugrd;
    E_gl = 0.5 * (dfgrd.transpose() * dfgrd - del);

    // Compute PK2 Stress using St. Venant Kirchoff Relations
    double tr_E_gl = E_gl.trace();
    for (unsigned int jdim = 0; jdim < dim; ++jdim)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
            S_pk2(idim, jdim) = 2 * mat.mu * E_gl(idim, jdim) + mat.lambda * del(idim, jdim) * tr_E_gl;
    }

    // Alternate way by C_ijkl E_kl
    //S_pk2 = Eigen::MatrixXd::Zero(dim, dim);
    //for (unsigned int jdim = 0; jdim < dim; ++jdim)
    //{
    //    for (unsigned int idim = 0; idim < dim; ++idim)
    //    {
    //        for (unsigned int ldim = 0; ldim < dim; ++ldim)
    //        {
    //            for (unsigned int kdim = 0; kdim < dim; ++kdim)
    //                S_pk2(idim, jdim) += mat.C[idim][jdim][kdim][ldim] * E_gl(kdim, ldim);
    //        }
    //    }
    //}
}

template<unsigned int dim>
void FSElasticityFEM<dim>::compute_quad_res(unsigned int iq, Eigen::MatrixXd& dNdX, Eigen::MatrixXd& dfgrd,
                                            Eigen::MatrixXd& S_pk2, double det_jac, Eigen::VectorXd& Rq)
{
    unsigned int loc_dofs = dim * mesh.Nen;

    // Contributions from PK2 stress
    // Contributions from body forces (not implemented)
    Eigen::VectorXd rq(loc_dofs);
    for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            unsigned int loc_idx = dim * ashp + idim;
            double temp = 0.0;
            for (unsigned int Jdim = 0; Jdim < dim; ++Jdim)
            {
                for (unsigned int Kdim = 0; Kdim < dim; ++Kdim)
                    temp += dNdX(ashp, Kdim) * dfgrd(idim, Jdim) * S_pk2(Jdim, Kdim);
            }

            rq(loc_idx) = temp * mesh.wts(iq) * det_jac;
        }
    }

    // Contribuitions from tractions (Not implemented)
    Eigen::VectorXd sq = Eigen::VectorXd::Zero(loc_dofs);

    // Total Residual
    Rq = rq - sq;
}

template<unsigned int dim>
void FSElasticityFEM<dim>::compute_quad_stiff(unsigned int iq, Eigen::MatrixXd& dNdX, Eigen::MatrixXd& dfgrd,
                                              Eigen::MatrixXd& S_pk2, double det_jac, Eigen::MatrixXd& Kjq)
{
    unsigned int loc_dofs = dim * mesh.Nen;

    // Compute system geometric stiffness
    Eigen::MatrixXd Kg = Eigen::MatrixXd::Zero(loc_dofs, loc_dofs);
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                unsigned int loc_idx = dim * ashp + idim;
                unsigned int loc_jdx = dim * bshp + idim;
                double temp = 0.0;
                for (unsigned int Idim = 0; Idim < dim; ++Idim)
                {
                    for (unsigned int Jdim = 0; Jdim < dim; ++Jdim)
                        temp += dNdX(ashp, Jdim) * S_pk2(Idim, Jdim) * dNdX(bshp, Idim);
                }

                Kg(loc_idx, loc_jdx) = temp * mesh.wts(iq) * det_jac;
            }
        }
    }

    // Compute system material stiffness
    Eigen::MatrixXd Km = Eigen::MatrixXd::Zero(loc_dofs, loc_dofs);
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int jdim = 0; jdim < dim; ++jdim)
        {
            unsigned int loc_jdx = dim * bshp + jdim;
            for (unsigned int ashp = 0; ashp < mesh.Nen; ++ashp)
            {
                for (unsigned int idim = 0; idim < dim; ++idim)
                {
                    unsigned int loc_idx = dim * ashp + idim;
                    double temp = 0.0;
                    for (unsigned int Idim = 0; Idim < dim; ++Idim)
                    {
                        for (unsigned int Jdim = 0; Jdim < dim; ++Jdim)
                        {
                            for (unsigned int Kdim = 0; Kdim < dim; ++Kdim)
                            {
                                for (unsigned int Ldim = 0; Ldim < dim; ++Ldim)
                                {
                                    temp += dNdX(ashp, Jdim) * dfgrd(idim, Idim) *
                                        mat.C[Idim][Jdim][Kdim][Ldim] * dfgrd(jdim, Ldim) * dNdX(bshp, Kdim);
                                }
                            }
                        }
                    }

                    Km(loc_idx, loc_jdx) = temp * mesh.wts(iq) * det_jac;
                }
            }
        }
    }

    // Compute system total stiffness
    Kjq = Kg + Km;
}

template<unsigned int dim>
void FSElasticityFEM<dim>::create_dirichlet_map()
{
    // dirBCs should be sorted for this to work currently

    gDofs = dim * mesh.Nnd;
    unsigned int dirBCsize = static_cast<unsigned int>(dirBCs.size());
    ucDofs = gDofs - dirBCsize;       // unconstrained dofs
    dirmap.reserve(ucDofs);
    unsigned int offset = 0;
    auto dir_itr = dirBCs.begin();
    for (unsigned int idx = 0; idx < mesh.Nnd; ++idx)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            if (dir_itr != dirBCs.end())
            {
                if (idx == dir_itr->node && idim == dir_itr->dof)
                {
                    ++dir_itr;
                    continue;
                }
            }

            dirmap.push_back(dim * idx + idim);
        }
    }
}

template<unsigned int dim>
void FSElasticityFEM<dim>::apply_step_BC(double step_mult)
{
    for (auto& dirBC : dirBCs)
    {
        unsigned gidx = dirBC.node * dim + dirBC.dof;
        d(gidx) = dirBC.val * step_mult;
    }
}

template<unsigned int dim>
void FSElasticityFEM<dim>::apply_dirichlet_BC()
{
    for (unsigned int jdx = 0; jdx < ucDofs; ++jdx)
    {
        for (unsigned int idx = 0; idx < ucDofs; ++idx)
            Kjd(idx, jdx) = Kj(dirmap[idx], dirmap[jdx]);

        Rd(jdx) = R(dirmap[jdx]);
    }
}

template<unsigned int dim>
void FSElasticityFEM<dim>::check_convergence(unsigned int nr_iter, bool& converged)
{
    converged = false;
    double norm_res = Rd.norm();
    std::cout << "NR iteration number = " << nr_iter + 1 << "\n";
    std::cout << "Residual norm = " << norm_res << "\n";
    if (norm_res < sol_ctrls.nr_tol)
    {
        converged = true;
        std::cout << "*** Convergence Successful!\n";
    }
}

template<unsigned int dim>
void FSElasticityFEM<dim>::reform_full_sol()
{
    for (unsigned int idx = 0; idx < ucDofs; ++idx)
        d(dirmap[idx]) = dd(idx);
}

template<unsigned int dim>
void FSElasticityFEM<dim>::raw_sol_print()
{
    std::ofstream nc("nodal_coords.txt");
    if (nc.is_open())
    {
        for (unsigned int idx = 0; idx < mesh.Nnd - 1; ++idx)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
                nc << std::setprecision(16) << mesh.nodal_coords[idx][idim] << "\t";

            nc << "\n";
        }

        for (unsigned int idim = 0; idim < dim; ++idim)
            nc << std::setprecision(16) << mesh.nodal_coords[mesh.Nnd - 1][idim] << "\t";
        nc.close();
    }

    std::ofstream ec("elem_conn.txt");
    if (ec.is_open())
    {
        for (unsigned int idx = 0; idx < mesh.Nel - 1; ++idx)
        {
            for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
                ec << std::setprecision(16) << mesh.elem_conn[idx][ind] << "\t";

            ec << "\n";
        }

        for (unsigned int ind = 0; ind < mesh.Nen; ++ind)
            ec << std::setprecision(16) << mesh.elem_conn[mesh.Nel - 1][ind] << "\t";
        ec.close();
    }

    std::ofstream sol("sol.txt");
    if (sol.is_open())
    {
        for (unsigned int ind = 0; ind < mesh.Nnd - 1; ++ind)
        {
            for (unsigned int idim = 0; idim < dim; ++idim)
            {
                unsigned int gidx = dim * ind + idim;
                sol << std::setprecision(16) << d(gidx) << "\t";
            }

            sol << "\n";
        }

        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            unsigned int gidx = dim * (mesh.Nnd - 1) + idim;
            sol << std::setprecision(16) << d(gidx) << "\t";
        }

        sol.close();
    }
}
