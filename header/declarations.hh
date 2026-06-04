/**
 * @file declarations.hh
 * @brief Header file containing includes, type definitions, and function declarations.
 */
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

/**
 * @brief Solves the Laplace equation using a sequential Jacobi iteration.
 * * @param n Number of grid points along one dimension.
 * @param tol Tolerance for the stopping criterion based on the L2 norm of the difference.
 * @param f The forcing function of the Laplace equation.
 * @return Matrix containing the computed solution over the grid.
 */
Matrix sequential_jacobi(int n, double tol, std::function<double(double,double)>);

/**
 * @brief Solves the Laplace equation using a hybrid parallel (MPI + OpenMP) Jacobi iteration.
 * * @param n Number of grid points along one dimension.
 * @param tol Tolerance for the stopping criterion based on the L2 norm of the difference.
 * @param f The forcing function of the Laplace equation.
 * @return Matrix containing the gathered solution on rank 0, or an empty/partial matrix on other ranks.
 */
Matrix jacobi(int n, double tol, std::function<double(double,double)>);

/**
 * @brief Saves the solution matrix to a VTK file.
 * * @param filename The name of the output VTK file.
 * @param u The solution matrix to be saved.
 * @param n Number of grid points along one dimension.
 */
void save_vtk_xml(const std::string& filename, const Matrix& u, int n);