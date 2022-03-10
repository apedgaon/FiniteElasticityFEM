#include <iostream>
#include "FSElasticityFEM.h"


void Problem1()
{
    FSElasticityFEM<3> fem;
    fem.geom.length[0] = 0.1;
    fem.geom.length[1] = 0.03;
    fem.geom.length[2] = 0.03;

    fem.mat.lambda = 6.0e10;
    fem.mat.mu = 2.0e10;
    fem.mat.init_elasticity_tensor();

    fem.mesh.mesh_type = MeshType::hexahedral;
    fem.mesh.quadrature_type = QuadratureType::two_point;
    fem.mesh.shp_func_type = ShapeFcnType::lagrange_linear;
    fem.mesh.divs[0] = 10;
    fem.mesh.divs[1] = 3;
    fem.mesh.divs[2] = 3;

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
        else if (fabs(nd_coords[0] - 0.1) < 1.0e-8)
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 0;
            temp.val = 0.01;
            fem.dirBCs.push_back(temp);
        }
        
        if ((fabs(nd_coords[2] - 0.03) < 1.0e-8) && (fabs(nd_coords[0] - 0.0) > 1.0e-8))
        {
            DirichletBC temp;
            temp.node = idx;
            temp.dof = 1;
            temp.val = 0.0;
            fem.dirBCs.push_back(temp);
        }
    }

    fem.solve();
}

int main(int argc, char* argv[])
{
    Problem1();
    return 0;
}