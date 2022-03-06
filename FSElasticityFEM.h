#pragma once
#include "FEMUtility.h"

// Finite Strain Elasticity Finite Element Method Class
template <unsigned int dim>
class FSElasticityFEM
{
public:
    void generate_mesh();
    void assemble();

    Geometry<dim> geom;
    Material mat;
    Mesh<dim> mesh;
    std::vector<DirichletBC> dirBCs;

private:
    void compute_elem_quantities(connectivity& conn, Eigen::VectorXui& gidx, Eigen::MatrixXd& xe, Eigen::VectorXd& de);
    void compute_quad_derivatives(unsigned int iq, Eigen::MatrixXd& xe, Eigen::VectorXd& de, Eigen::MatrixXd& jac,
                                  Eigen::MatrixXd& jac_inv, Eigen::MatrixXd& dNdX, double& det_jac);
    void compute_quad_deformation(Eigen::VectorXd& de, Eigen::MatrixXd& dNdX, Eigen::MatrixXd& dfgrd, Eigen::MatrixXd& E_gl,
                                  Eigen::MatrixXd& S_pk2);
    void compute_quad_res(unsigned int iq, Eigen::MatrixXd& dNdX, Eigen::MatrixXd& dfgrd, Eigen::MatrixXd& S_pk2,
                          double det_jac, Eigen::MatrixXd& Rq);

private:
    Eigen::VectorXd d;
    Eigen::VectorXd dd;

    Eigen::VectorXd R;
    Eigen::VectorXd Rd;
};