This is C++ Finite Element code for solving 3D Finite Strain Elasticity. The material is assumed to be hyperelastic with
St. Venant Kirchoff Model.

This code is dependent on Eigen 3.4.0 which needs to be installed before building.
	 https://eigen.tuxfamily.org/index.php?title=Main_Page

A path to downloaded eigen directory should be provided to cmake. 

BUILD INTRUCTIONS:

mkdir build
cd build
cmake -DEIGEN_INCLUDE=/path/to/eigen -DCMAKE_BUILD_TYPE=Release ..
cmake --build .



CODE STRUCTURE:
- FSProblem.cpp - Contains MAIN function and Problem1 and Problem2 definitions
- FEMUtility.h, FEMUtility.cpp - Contains Shape fcn, quadrature, elasticty tensor,
				 mesh generation, printing and vtk utility
- FSElasticity.h, FSElasticity.cpp - Contains Finite Strain Elasticity Class declaration
				     and implementations 
				     (Load Steps, Newton Raphson, Assembly, Green Lagrange Strain,
				      PK2 stress, residual, system stiffness, etc.)
