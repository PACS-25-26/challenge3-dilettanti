/**
 * @file main.cpp
 *
 * This file initializes the MPI environment, sets up the problem parameters,
 * and calls both the sequential and parallel Jacobi solvers to compare their
 * execution times. It also handles exporting the solutions to VTK format.
 */
#include "declarations.hh"
#include "chrono.hpp"


int main (int argc, char *argv[]){
     MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Timings::Chrono timer;

    int n = 100;                 // Number of grid points along one dimension
    double tolleranza = 1e-5;    // Stopping criterion tolerance for the Jacobi solver
    int test_lenght = 20;        // Number of iterations to average the execution time
  
    /**
     * @brief The forcing function f(x,y) for the Laplace equation.
     */
    auto f = [](double x, double y){ 
        return 8*PI*PI*sin(2*PI*x)*sin(2*PI*y);
    };

    timer.start();
    Matrix U_seq;
    for(int i=0; i<test_lenght; ++i)
        U_seq = sequential_jacobi(n, tolleranza, f);
    
    timer.stop();
    if(rank == 0){
        std::cout << "Tempo medio Sequenziale: " << timer.wallTime() / test_lenght << " microsec" << std::endl;  
        save_vtk_xml("sequential_solution.vti", U_seq, n);  
    }

    timer.start();
    Matrix U_par;
    for(int i=0; i<test_lenght; ++i)
        U_par = jacobi(n, tolleranza, f);

    timer.stop();
    if(rank == 0){
        std::cout << "Tempo medio Parallelo: " << timer.wallTime() / test_lenght << " microsec" << std::endl;
        save_vtk_xml("parallel_solution.vti", U_par, n);  
    }

    MPI_Finalize();

    return 0;
}
