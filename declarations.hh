#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <mpi.h>
#include <stdexcept>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>

#define PI 3.14159265358979323846

using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

void jacobi(int n, double tol);

void save_vtk_xml(const std::string& filename, const Matrix& u, double h);