int main (int argc, char *argv[]){
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Lettura da riga di comando (default a 100 se non fornito)
    int dimensione_griglia = (argc > 1) ? std::stoi(argv[1]) : 100;
    double tolleranza = (argc > 2) ? std::stod(argv[2]) : 1e-5;

    Timings::Chrono timer;
    if(rank == 0) timer.start();

    jacobi(dimensione_griglia, tolleranza);

    if(rank == 0) {
        timer.stop();
        std::cout << timer;
    }

    MPI_Finalize();
    return 0;
}
