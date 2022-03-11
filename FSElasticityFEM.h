#pragma once
#include <vector>
#include <Eigen/Dense>
#include "FEMUtility.h"

// Finite Strain Elasticity Finite Element Method Class
template <unsigned int dim>
class FSElasticityFEM
{
public:
    void generate_mesh();
    void initialize();
    void assemble();
    void solve();
    void post_print();

    Geometry<dim> geom;
    Material mat;
    Mesh<dim> mesh;
    std::vector<DirichletBC> dirBCs;
    SolverControls sol_ctrls;

private:
    void compute_elem_quantities(connectivity& conn, Eigen::VectorXui& gidx, Eigen::MatrixXd& xe, Eigen::VectorXd& de);
    void compute_quad_derivatives(unsigned int iq, Eigen::MatrixXd& xe, Eigen::VectorXd& de, Eigen::MatrixXd& jac,
                                  Eigen::MatrixXd& jac_inv, Eigen::MatrixXd& dNdX, double& det_jac);
    void compute_quad_deformation(Eigen::VectorXd& de, Eigen::MatrixXd& dNdX, Eigen::MatrixXd& dfgrd, Eigen::MatrixXd& E_gl,
                                  Eigen::MatrixXd& S_pk2);
    void compute_quad_res(unsigned int iq, Eigen::MatrixXd& dNdX, Eigen::MatrixXd& dfgrd, Eigen::MatrixXd& S_pk2,
                          double det_jac, Eigen::MatrixXd& Rq);
    void compute_quad_stiff(unsigned int iq, Eigen::MatrixXd& dNdX, Eigen::MatrixXd& dfgrd, Eigen::MatrixXd& S_pk2,
                            double det_jac, Eigen::MatrixXd& Kjq);
    void create_dirichlet_map();
    void apply_step_BC(double step_mult);
    void apply_dirichlet_BC();
    void check_convergence(unsigned int nr_iter, bool& converged);
    void reform_full_sol();

private:
    Eigen::VectorXd d;
    Eigen::VectorXd dd;

    Eigen::VectorXd R;
    Eigen::VectorXd Rd;

    Eigen::MatrixXd Kj;
    Eigen::MatrixXd Kjd;

    unsigned int gDofs;
    unsigned int ucDofs;
    std::vector<unsigned int> dirmap;
};