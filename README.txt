This C++ code is dependent on Eigen 3.4.0 and will need it downloaded.
	 https://eigen.tuxfamily.org/index.php?title=Main_Page

A path to downloaded eigen directory should be provided to cmake. 

BUILD INTRUCTIONS:

mkdir build
cd build
cmake -DEIGEN_INCLUDE=/path/to/eigen -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
