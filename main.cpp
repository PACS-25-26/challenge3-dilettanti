#include "declarations.hh"
#include "chrono.hpp"


int main (int argc, char *argv[]){

    Timings::Chrono timer;
    int dimensione_griglia = 100;
    double tolleranza = 1e-5;

    timer.start();

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    jacobi(dimensione_griglia, tolleranza);

    MPI_Finalize();

    timer.stop();
    std::cout << timer;
    return 0;
}
