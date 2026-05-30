#include "declarations.hh"


void jacobi(int n, double tol){

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto f = [](double x, double y){
        return 8 * PI * PI * std::sin(2 * PI * x) * std::sin(2 * PI * y);
    };

    int local_n = n / size;
    int resto = n % size;
    if(rank < resto)
        ++local_n;
    int local_dim = local_n * n;

    std::vector<int> displ_vector(size);
    std::vector<int> sendcount_vector(size);

    displ_vector[0] = 0;
    for(int i = 1; i < size; ++i){
        int loc_n = n / size;
        int r = n % size;
        if(i - 1 < r)
            ++loc_n;
        displ_vector[i] = n * loc_n + displ_vector[i-1];
    }
    for(int i = 0; i < size; ++i){
        int loc_n = n / size;
        int r = n % size;
        if(i < r)
            ++loc_n;
        sendcount_vector[i] = n * loc_n;
    }

    Matrix U;

    if(rank == 0)
        U = Matrix::Zero(n, n);

    Matrix local_U(local_n + 2, n);
    local_U.setZero();

    MPI_Scatterv(U.data(), sendcount_vector.data(), displ_vector.data(), MPI_DOUBLE, local_U.row(1).data(), local_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double h = 1.0 / (n - 1);
    Matrix old_local_U = local_U;
    double eps = 1.0;
    int check = 0;

    int top_neigbor, bottom_neigbor;
    if(rank == 0){
        top_neigbor = MPI_PROC_NULL;
        bottom_neigbor = 1;
    } else if(rank == size - 1){
        top_neigbor = size - 2;
        bottom_neigbor = MPI_PROC_NULL;
    } else {
        top_neigbor = rank - 1;
        bottom_neigbor = rank + 1;
    }

    int iter = 0;
    // Modificato il controllo inserendo un limite massimo di iterazioni per sicurezza (es. 100000)
    while(check != 1 && iter < 100000){
        check = 0;

        // step 1: comunicazione
        MPI_Sendrecv(local_U.row(1).data(), n, MPI_DOUBLE, top_neigbor, 0,
                     old_local_U.row(0).data(), n, MPI_DOUBLE, top_neigbor, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(local_U.row(local_n).data(), n, MPI_DOUBLE, bottom_neigbor, 0,
                     old_local_U.row(local_n+1).data(), n, MPI_DOUBLE, bottom_neigbor, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // step 2: calcolo
        #pragma omp parallel for collapse(2)
        for(int i = 1; i <= local_n; ++i){
            for(int j = 1; j < n - 1; ++j){
                double x = h * j;
                int global_i = (displ_vector[rank] / n + (i - 1));
                if(global_i == 0 || global_i == n - 1)
                    continue;
                double y = h * global_i;
                local_U(i,j) = 0.25 * (old_local_U(i-1,j) + old_local_U(i+1,j) + old_local_U(i,j-1) + old_local_U(i,j+1) + h * h * f(x,y));
            }
        }

        // step 3: check condizione
        double local_sq_err = 0.0;

        #pragma omp parallel for collapse(2) reduction(+:local_sq_err)
        for(int i = 1; i <= local_n; ++i){
            for(int j = 0; j < n; ++j){
                double diff = local_U(i,j) - old_local_U(i,j);
                local_sq_err += diff * diff;
            }
        }

        double global_sq_err = 0.0;
        MPI_Allreduce(&local_sq_err, &global_sq_err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        eps = std::sqrt(h * global_sq_err);

        if(eps < tol){
            check = 1;
        }

        old_local_U = local_U;
        iter++;
    }

    MPI_Gatherv(local_U.row(1).data(), local_dim, MPI_DOUBLE,
                U.data(), sendcount_vector.data(), displ_vector.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    // Solo il rank 0 salva il file VTK e calcola l'errore
    if (rank == 0) {
        // Calcolo errore esatto L2
        auto u_exact = [](double x, double y){
            return std::sin(2 * PI * x) * std::sin(2 * PI * y);
        };

        double err_L2 = 0.0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double diff = U(i, j) - u_exact(h * j, h * i); // Attenzione: j è la coordinata x, i è la coordinata y
                err_L2 += diff * diff;
            }
        }
        err_L2 = std::sqrt(h * err_L2);

        std::cout << "Convergenza raggiunta in " << iter << " iterazioni." << std::endl;
        std::cout << "Dimensione griglia (n): " << n << std::endl;
        std::cout << "Errore L2 rispetto alla soluzione esatta: " << std::scientific << err_L2 << std::endl;

        save_vtk_xml("poisson_mpi_solution.vti", U, h);
    }
}


void save_vtk_xml(const std::string& filename, const Matrix& u, double h) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Errore nell'apertura del file VTK!" << std::endl;
        return;
    }

    // Numero di punti lungo X e Y (u è N x N)
    int nX = u.rows();
    int nY = u.cols();

    // Scrittura dell'intestazione XML per griglie regolari
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    // WholeExtent definisce la griglia index-based: da 0 a N-1 per X e Y, 0 per Z (2D)
    file << "  <ImageData WholeExtent=\"0 " << nX-1 << " 0 " << nY-1 << " 0 0\" ";
    file << "Origin=\"0.0 0.0 0.0\" Spacing=\"" << h << " " << h << " 0.0\">\n";
    file << "    <Piece Extent=\"0 " << nX-1 << " 0 " << nY-1 << " 0 0\">\n";
    
    file << "      <PointData Scalars=\"Soluzione\">\n";
    file << "        <DataArray type=\"Float64\" Name=\"Soluzione\" format=\"ascii\">\n";

    file << std::scientific << std::setprecision(6);
    
    for (int i = 0; i < nX; ++i) {
        for (int j = 0; j < nY; ++j) {
            file << u(i, j) << " ";
        }
        file << "\n";
    }

    file << "        </DataArray>\n";
    file << "      </PointData>\n";
    file << "    </Piece>\n";
    file << "  </ImageData>\n";
    file << "</VTKFile>\n";

    file.close();
    std::cout << "File VTK salvato con successo: " << filename << std::endl;
}





