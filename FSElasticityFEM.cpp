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
    create_dirichlet_map();

    d = Eigen::VectorXd::Zero(gDofs);
    for (auto& dirBC : dirBCs)
    {
        unsigned gidx = dirBC.node * dim + dirBC.dof;
        d(gidx) = dirBC.val;
    }

    dd.resize(ucDofs);
    Rd.resize(ucDofs);
    Kjd.resize(ucDofs, ucDofs);

    for (unsigned int idx = 0; idx < ucDofs; ++idx)
        dd(idx) = d(dirmap[idx]);
}

template<unsigned int dim>
void FSElasticityFEM<dim>::assemble()
{
    R = Eigen::VectorXd::Zero(gDofs);
    Kj = Eigen::MatrixXd::Zero(gDofs, gDofs);
    unsigned int loc_dofs = mesh.Nen * dim;

    // Element Loop
    for (auto& conn : mesh.elem_conn)
    {
        Eigen::VectorXd Rl = Eigen::VectorXd::Zero(loc_dofs);
        Eigen::MatrixXd Kjl = Eigen::MatrixXd::Zero(loc_dofs, loc_dofs);

        // Compute elemental coordinates and displacements
        Eigen::MatrixXd xe(dim, mesh.Nen);
        Eigen::VectorXui gidx(loc_dofs);
        Eigen::VectorXd de(mesh.Nen * dim);
        compute_elem_quantities(conn, gidx, xe, de);

        // Quadrature loop
        for (unsigned int iq = 0; iq < mesh.Nq; ++iq)
        {
            // Compute quadrature jacobian and derivatives
            Eigen::MatrixXd jac, jac_inv, dNdX;
            double det_jac;
            compute_quad_derivatives(iq, xe, de, jac, jac_inv, dNdX, det_jac);

            // Compute quadrature deformation gradient, strain and stress
            Eigen::MatrixXd dfgrd, E_gl, S_pk2;
            compute_quad_deformation(de, dNdX, dfgrd, E_gl, S_pk2);

            // Compute quadrature Residual
            Eigen::MatrixXd Rq;
            compute_quad_res(iq, dNdX, dfgrd, S_pk2, det_jac, Rq);
            
            // Update elemental Residual 
            Rl += Rq;

            // Compute quadrature NR Stiffness
            Eigen::MatrixXd Kjq;
            compute_quad_stiff(iq, dNdX, dfgrd, S_pk2, det_jac, Kjq);

            // Update elemental NR Stiffness
            Kjl += Kjq;
        }

        // Global Assembly Rl -> R
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
    Eigen::VectorXd delta_d(ucDofs);
    double norm_res;

    // NR iterations
    for (unsigned int nr_idx = 0; nr_idx < 6; ++nr_idx)
    {
        assemble();
        apply_dirichlet_BC();
        norm_res = Rd.norm();
        delta_d = Kjd.ldlt().solve(-Rd);
        dd += delta_d;
        reform_full_sol();
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
    // Isoparametric Jacobian and Inverse
    jac = Eigen::MatrixXd::Zero(dim, dim);
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
    S_pk2.resize(dim, dim);
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
                                            Eigen::MatrixXd& S_pk2, double det_jac, Eigen::MatrixXd& Rq)
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
            for (unsigned int kdim = 0; kdim < dim; ++kdim)
            {
                for (unsigned int ldim = 0; ldim < dim; ++ldim)
                    temp += dfgrd(idim, kdim) * S_pk2(kdim, ldim) * dNdX(ashp, ldim);
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
                        temp += dNdX(bshp, Idim) * S_pk2(Idim, Jdim) * dNdX(ashp, Jdim);
                }

                Kg(loc_idx, loc_jdx) = temp * mesh.wts(iq) * det_jac;
            }
        }
    }

    // Compute system material stiffness
    Eigen::MatrixXd Km = Eigen::MatrixXd::Zero(loc_dofs, loc_dofs);
    for (unsigned int bshp = 0; bshp < mesh.Nen; ++bshp)
    {
        for (unsigned int kdim = 0; kdim < dim; ++kdim)
        {
            unsigned int loc_kdx = dim * bshp + kdim;
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
                                    temp += dfgrd(idim, Idim) * mat.C[Idim][Jdim][Kdim][Ldim] *
                                            dfgrd(kdim, Kdim) * dNdX(bshp, Ldim) * dNdX(ashp, Jdim);
                                }
                            }
                        }
                    }

                    Km(loc_idx, loc_kdx) = temp * mesh.wts(iq) * det_jac;
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
    // dirBCs should be sorted
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
void FSElasticityFEM<dim>::reform_full_sol()
{
    for (unsigned int idx = 0; idx < ucDofs; ++idx)
        d(dirmap[idx]) = dd(idx);
}
