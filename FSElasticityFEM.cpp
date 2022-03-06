#include "FSElasticityFEM.h"

template class FSElasticityFEM<3>;



template<unsigned int dim>
inline void FSElasticityFEM<dim>::generate_mesh()
{
    mesh.generate(geom);
}


template<unsigned int dim>
void FSElasticityFEM<dim>::assemble()
{
    // Initialize
    unsigned int num_dofs = mesh.Nnd * dim;
    unsigned int loc_dofs = mesh.Nen * dim;
    R = Eigen::VectorXd::Zero(num_dofs);
    d = Eigen::VectorXd::Zero(num_dofs);
    for (auto& dirBC : dirBCs)
    {
        unsigned gidx = dirBC.node * dim + dirBC.dof;
        d(gidx) = dirBC.val;
    }

    // Element Loop
    for (auto& conn : mesh.elem_conn)
    {
        Eigen::VectorXd Rl = Eigen::VectorXd::Zero(loc_dofs);

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
        }

        // Global Assembly Rl -> R
        for (unsigned int idx = 0; idx < loc_dofs; ++idx)
        {
            R(gidx(idx)) += Rl(idx);
        }
    }

    auto bb = R.norm();
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

    // Compute PK2 Stress
    S_pk2 = Eigen::MatrixXd::Zero(dim, dim);
    for (unsigned int jdim = 0; jdim < dim; ++jdim)
    {
        for (unsigned int idim = 0; idim < dim; ++idim)
        {
            for (unsigned int ldim = 0; ldim < dim; ++ldim)
            {
                for (unsigned int kdim = 0; kdim < dim; ++kdim)
                    S_pk2(idim, jdim) += mat.C[idim][jdim][kdim][ldim] * E_gl(kdim, ldim);
            }
        }
    }
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

    Rq = rq - sq;
}
