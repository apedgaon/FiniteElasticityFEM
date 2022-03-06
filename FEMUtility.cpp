#include "FEMUtility.h"

template class Mesh<3>;

double LinearLagrange::N(unsigned int index, double x)
{
    double val = 0.0;
    if (index == 0)
        val = 0.5 * (1 - x);
    else if (index == 1)
        val = 0.5 * (1 + x);

    return val;
}

double LinearLagrange::dNdx(unsigned int index, double x)
{
    double val = 0.0;
    if (index == 0)
        val = -0.5;
    else if (index == 1)
        val = 0.5;

    return val;
}

double LinearLagrange::N_quad(unsigned int index, double xi, double eta)
{
    double val = 0.0;
    if (index == 0)
        val = 0.25 * (1 - xi) * (1 - eta);
    else if (index == 1)
        val = 0.25 * (1 + xi) * (1 - eta);
    else if (index == 2)
        val = 0.25 * (1 + xi) * (1 + eta);
    else if (index == 3)
        val = 0.25 * (1 - xi) * (1 + eta);

    return val;
}

double LinearLagrange::dNdxi_quad(unsigned int index, double xi, double eta)
{
    double val = 0.0;
    if (index == 0)
        val = -0.25 * (1 - eta);
    else if (index == 1)
        val = 0.25 * (1 - eta);
    else if (index == 2)
        val = 0.25 * (1 + eta);
    else if (index == 3)
        val = -0.25 * (1 + eta);

    return val;
}

double LinearLagrange::dNdeta_quad(unsigned int index, double xi, double eta)
{
    double val = 0.0;
    if (index == 0)
        val = -0.25 * (1 - xi);
    else if (index == 1)
        val = -0.25 * (1 + xi);
    else if (index == 2)
        val = 0.25 * (1 + xi);
    else if (index == 3)
        val = 0.25 * (1 - xi);

    return val;
}

double LinearLagrange::N_tri(unsigned int index, double xi, double eta)
{
    double val = 0.0;
    if (index == 0)
        val = 1.0 - xi - eta;
    if (index == 1)
        val = xi;
    else if (index == 2)
        val = eta;

    return val;
}

double LinearLagrange::dNdxi_tri(unsigned int index, double xi, double eta)
{
    double val = 0.0;
    if (index == 0)
        val = -1.0;
    else if (index == 1)
        val = 1.0;
    else if (index == 2)
        val = 0.0;

    return val;
}

double LinearLagrange::dNdeta_tri(unsigned int index, double xi, double eta)
{
    double val = 0.0;
    if (index == 0)
        val = -1.0;
    else if (index == 1)
        val = 0.0;
    else if (index == 2)
        val = 1.0;

    return val;
}

double LinearLagrange::N_hex(unsigned int index, double xi, double eta, double kappa)
{
    double val = 0.0;
    if (index == 0)
        val = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 - kappa);
    else if (index == 1)
        val = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 - kappa);
    else if (index == 2)
        val = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 - kappa);
    else if (index == 3)
        val = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 - kappa);
    else if (index == 4)
        val = 0.125 * (1.0 - xi) * (1.0 - eta) * (1.0 + kappa);
    else if (index == 5)
        val = 0.125 * (1.0 + xi) * (1.0 - eta) * (1.0 + kappa);
    else if (index == 6)
        val = 0.125 * (1.0 + xi) * (1.0 + eta) * (1.0 + kappa);
    else if (index == 7)
        val = 0.125 * (1.0 - xi) * (1.0 + eta) * (1.0 + kappa);

    return val;
}

double LinearLagrange::dNdxi_hex(unsigned int index, unsigned int diff_idx, double xi, double eta, double kappa)
{
    double val = 0.0;
    if (diff_idx == 0) // dN/dxi
    {
        if (index == 0)
            val = -0.125 * (1.0 - eta) * (1.0 - kappa);
        else if (index == 1)
            val = 0.125 * (1.0 - eta) * (1.0 - kappa);
        else if (index == 2)
            val = 0.125 * (1.0 + eta) * (1.0 - kappa);
        else if (index == 3)
            val = -0.125 * (1.0 + eta) * (1.0 - kappa);
        else if (index == 4)
            val = -0.125 * (1.0 - eta) * (1.0 + kappa);
        else if (index == 5)
            val = 0.125 * (1.0 - eta) * (1.0 + kappa);
        else if (index == 6)
            val = 0.125 * (1.0 + eta) * (1.0 + kappa);
        else if (index == 7)
            val = -0.125 * (1.0 + eta) * (1.0 + kappa);
    }
    else if (diff_idx == 1) // dN/dEta
    {
        if (index == 0)
            val = -0.125 * (1.0 - xi) * (1.0 - kappa);
        else if (index == 1)
            val = -0.125 * (1.0 + xi) * (1.0 - kappa);
        else if (index == 2)
            val = 0.125 * (1.0 + xi) * (1.0 - kappa);
        else if (index == 3)
            val = 0.125 * (1.0 - xi) * (1.0 - kappa);
        else if (index == 4)
            val = -0.125 * (1.0 - xi) * (1.0 + kappa);
        else if (index == 5)
            val = -0.125 * (1.0 + xi) * (1.0 + kappa);
        else if (index == 6)
            val = 0.125 * (1.0 + xi) * (1.0 + kappa);
        else if (index == 7)
            val = 0.125 * (1.0 - xi) * (1.0 + kappa);
    }
    else if (diff_idx == 2) // dN/dKappa
    {
        if (index == 0)
            val = -0.125 * (1.0 - xi) * (1.0 - eta);
        else if (index == 1)
            val = -0.125 * (1.0 + xi) * (1.0 - eta);
        else if (index == 2)
            val = -0.125 * (1.0 + xi) * (1.0 + eta);
        else if (index == 3)
            val = -0.125 * (1.0 - xi) * (1.0 + eta);
        else if (index == 4)
            val = 0.125 * (1.0 - xi) * (1.0 - eta);
        else if (index == 5)
            val = 0.125 * (1.0 + xi) * (1.0 - eta);
        else if (index == 6)
            val = 0.125 * (1.0 + xi) * (1.0 + eta);
        else if (index == 7)
            val = 0.125 * (1.0 - xi) * (1.0 + eta);
    }

    return val;
}

double QuadraticLagrange::N(unsigned int index, double x)
{
    double val = 0.0;
    if (index == 0)
        val = -0.5 * x * (1 - x);
    else if (index == 1)
        val = (1 - x) * (1 + x);
    else if (index == 2)
        val = 0.5 * x * (1 + x);

    return val;
}

double QuadraticLagrange::dNdx(unsigned int index, double x)
{
    double val = 0.0;
    if (index == 0)
        val = -0.5 + x;
    else if (index == 1)
        val = -2.0 * x;
    else if (index == 2)
        val = 0.5 + x;

    return val;
}

double GsQuad::quadrature_pt(unsigned int num_pts, unsigned int index)
{
    double retval = 0.0;
    if (num_pts == 1)
        retval = 0.0;
    else if (num_pts == 2)
    {
        if (index == 0)
            retval = -1.0 / sqrt(3.0);
        else if (index == 1)
            retval = 1.0 / sqrt(3.0);
    }
    else if (num_pts == 3)
    {
        if (index == 0)
            retval = -sqrt(3.0 / 5.0);
        else if (index == 1)
            retval = 0.0;
        else if (index == 2)
            retval = sqrt(3.0 / 5.0);
    }

    return retval;
}

double GsQuad::weights(unsigned int num_pts, unsigned int index)
{
    double retval = 0.0;
    if (num_pts == 1)
        retval = 2.0;
    else if (num_pts == 2)
    {
        if (index == 0)
            retval = 1.0;
        else if (index == 1)
            retval = 1.0;
    }
    else if (num_pts == 3)
    {
        if (index == 0)
            retval = 5.0 / 9.0;
        else if (index == 1)
            retval = 8.0 / 9.0;
        else if (index == 2)
            retval = 5.0 / 9.0;
    }

    return retval;
}


template<unsigned int dim>
inline void Mesh<dim>::generate(Geometry<dim> geom)
{
    switch (quadrature_type)
    {
    case QuadratureType::two_point:
    {
        if (mesh_type == MeshType::hexahedral)
        {
            Nq = static_cast<unsigned int>(pow(2, dim));
            Qpts.resize(dim, Nq);
            Qpts(0, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(2, 0) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 1) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 1) = GsQuad::quadrature_pt(2, 0);
            Qpts(2, 1) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 2) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 2) = GsQuad::quadrature_pt(2, 1);
            Qpts(2, 2) = GsQuad::quadrature_pt(2, 0);
            Qpts(0, 3) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 3) = GsQuad::quadrature_pt(2, 1);
            Qpts(2, 3) = GsQuad::quadrature_pt(2, 0);

            Qpts(0, 4) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 4) = GsQuad::quadrature_pt(2, 0);
            Qpts(2, 4) = GsQuad::quadrature_pt(2, 1);
            Qpts(0, 5) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 5) = GsQuad::quadrature_pt(2, 0);
            Qpts(2, 5) = GsQuad::quadrature_pt(2, 1);
            Qpts(0, 6) = GsQuad::quadrature_pt(2, 1);
            Qpts(1, 6) = GsQuad::quadrature_pt(2, 1);
            Qpts(2, 6) = GsQuad::quadrature_pt(2, 1);
            Qpts(0, 7) = GsQuad::quadrature_pt(2, 0);
            Qpts(1, 7) = GsQuad::quadrature_pt(2, 1);
            Qpts(2, 7) = GsQuad::quadrature_pt(2, 1);

            wts.resize(Nq);
            wts(0) = GsQuad::weights(2, 0) * GsQuad::weights(2, 0) * GsQuad::weights(2, 0);
            wts(1) = GsQuad::weights(2, 1) * GsQuad::weights(2, 0) * GsQuad::weights(2, 0);
            wts(2) = GsQuad::weights(2, 1) * GsQuad::weights(2, 1) * GsQuad::weights(2, 0);
            wts(3) = GsQuad::weights(2, 0) * GsQuad::weights(2, 1) * GsQuad::weights(2, 0);

            wts(4) = GsQuad::weights(2, 0) * GsQuad::weights(2, 0) * GsQuad::weights(2, 1);
            wts(5) = GsQuad::weights(2, 1) * GsQuad::weights(2, 0) * GsQuad::weights(2, 1);
            wts(6) = GsQuad::weights(2, 1) * GsQuad::weights(2, 1) * GsQuad::weights(2, 1);
            wts(7) = GsQuad::weights(2, 0) * GsQuad::weights(2, 1) * GsQuad::weights(2, 1);
        }

    }break;
    }

    switch (shp_func_type)
    {
    case ShapeFcnType::lagrange_linear:
    {
        if (mesh_type == MeshType::hexahedral)
        {
            Nen = static_cast<unsigned int>(pow(2, dim));
            unsigned int Nnd_dim[dim];
            Nel = 1;  Nnd = 1;
            double h_dim[dim];
            for (unsigned int idx = 0; idx < dim; ++idx)
            {
                Nel *= divs[idx];
                Nnd_dim[idx] = divs[idx] + 1;
                Nnd *= Nnd_dim[idx];
                h_dim[idx] = geom.length[idx] / divs[idx];
            }

            nodal_coords.reserve(Nnd);
            for (unsigned int kdx = 0; kdx < Nnd_dim[2]; ++kdx)
            {
                for (unsigned int jdx = 0; jdx < Nnd_dim[1]; ++jdx)
                {
                    for (unsigned int idx = 0; idx < Nnd_dim[0]; ++idx)
                    {
                        coords coordinates(3);
                        coordinates[0] = idx * h_dim[0];
                        coordinates[1] = jdx * h_dim[1];
                        coordinates[2] = kdx * h_dim[2];
                        nodal_coords.push_back(coordinates);
                    }
                }
            }

            elem_conn.reserve(Nel);
            for (unsigned int kdx = 0; kdx < divs[2]; ++kdx)
            {
                for (unsigned int jdx = 0; jdx < divs[1]; ++jdx)
                {
                    for (unsigned int idx = 0; idx < divs[0]; ++idx)
                    {
                        connectivity conn(8);
                        conn[0] = kdx * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + idx;
                        conn[1] = kdx * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + (idx + 1);
                        conn[2] = kdx * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + (idx + 1);
                        conn[3] = kdx * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + idx;
                        conn[4] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + idx;
                        conn[5] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + jdx * Nnd_dim[0] + (idx + 1);
                        conn[6] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + (idx + 1);
                        conn[7] = (kdx + 1) * Nnd_dim[1] * Nnd_dim[0] + (jdx + 1) * Nnd_dim[0] + idx;
                        elem_conn.push_back(conn);
                    }
                }
            }

            N.resize(Nen, Nq);
            for (unsigned int idx = 0; idx < dim; ++idx)
                dNdxi[idx].resize(Nen, Nq);

            for (unsigned int jdx = 0; jdx < Nq; ++jdx)
            {
                for (unsigned int idx = 0; idx < Nen; ++idx)
                {
                    N(idx, jdx) = LinearLagrange::N_hex(idx, Qpts(0, jdx), Qpts(1, jdx), Qpts(2, jdx));
                    for (unsigned int idim = 0; idim < dim; ++idim)
                        dNdxi[idim](idx, jdx) = LinearLagrange::dNdxi_hex(idx, idim, Qpts(0, jdx), Qpts(1, jdx), Qpts(2, jdx));
                }
            }
        }
    }break;

    }
}

void Material::init_elasticity_tensor()
{
    Eigen::Matrix3d del = Eigen::Matrix3d::Identity();
    for (unsigned int idx = 0; idx < 3; ++idx)
    {
        for (unsigned int jdx = 0; jdx < 3; ++jdx)
        {
            for (unsigned int kdx = 0; kdx < 3; ++kdx)
            {
                for (unsigned int ldx = 0; ldx < 3; ++ldx)
                {
                    C[idx][jdx][kdx][ldx] = lambda * del(idx, jdx) * del(kdx, ldx) +
                        mu * (del(idx, kdx) * del(jdx, ldx) + del(idx, ldx) * del(jdx, kdx));
                }
            }
        }
    }
}
