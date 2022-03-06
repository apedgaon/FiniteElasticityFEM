cmake_minimum_required(VERSION 3.10)

project(FSElasticityFEM VERSION 1.0)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(HEADER_FILES FSElasticityFEM.h FEMUtility.h)
add_executable(FSElasticityFEM FSProblem.cpp FSElasticityFEM.cpp FEMUtility.cpp)
target_include_directories(FSElasticityFEM PUBLIC
							"${EIGEN_INCLUDE}")