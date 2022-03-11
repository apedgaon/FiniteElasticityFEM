#include <iostream>
#include "FSElasticityFEM.h"


void Problem1()
{
    FSElasticityFEM<3> fem;
    fem.geom.length[0] = 10.0;  // cm
    fem.geom.length[1] = 3.0;   // cm
    fem.geom.length[2] = 3.0;   // cm

    fem.mat.lambda = 6.0e6;     // N cm^-2
    fem.mat.mu = 2.0e6;         // N cm^-2
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::hexahedral;
    fem.mesh.quadrature_type = QuadratureType::two_point;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 20;
    fem.mesh.divs[1] = 6;
    fem.mesh.divs[2] = 6;

    fem.generate_mesh();
    for (unsigned int idx = 0; idx < fem.mesh.Nnd; ++idx)
    {
        auto& nd_coords = fem.mesh.nodal_coords[idx];
        if (fabs(nd_coords[0] - 0.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);

            temp.dof = 2;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }
        else if (fabs(nd_coords[0] - 10.0) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 1.0;     // cm
            fem.dirBCs.push_back(temp);
        }
        
        if ((fabs(nd_coords[2] - 3.0) < 1.0e-8) && (fabs(nd_coords[0] - 0.0) > 1.0e-8))
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }
    }

    fem.sol_ctrls.antype = AnalysisType::_static;
    fem.sol_ctrls.nsteps = 1;
    fem.sol_ctrls.ninc_max = 10;
    fem.sol_ctrls.nr_tol = 1.0e-8;

    fem.initialize();
    fem.solve();
    fem.post_print();
}

int main(int argc, char* argv[])
{
    Problem1();
    return 0;
}