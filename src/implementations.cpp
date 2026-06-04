/**
 * @file implementations.cpp
 * @brief Source file containing the implementations of the Jacobi solvers and I/O functions.
 */
#include "declarations.hh"

Matrix sequential_jacobi(int n, double tol, std::function<double(double,double)> f){
    Matrix U = Matrix::Zero(n, n);
    Matrix old_U = U;
    
    double h = 1.0/(n - 1);
    double eps = 1.0;
    int check = 0;

    while(check != 1){
        for(int i = 1; i < n - 1; ++i){
            for(int j = 1; j < n - 1; ++j){
                double x = h*j;
                double y = h*i;
                U(i, j) = 0.25*(old_U(i - 1, j) + old_U(i + 1, j) + old_U(i, j - 1) + old_U(i, j + 1) + h*h*f(x, y));
            }
        }

        eps = 0.0;
        for(int i = 1; i < n - 1; ++i){
            for(int j = 1; j < n - 1; ++j){
                double diff = U(i, j) - old_U(i, j);
                eps += diff*diff;
            }
        }
        eps = std::sqrt(h*eps);

        if(eps < tol)
            check = 1;

        old_U = U;
    }
    return U;
}


Matrix jacobi(int n, double tol, std::function<double(double,double)> f){

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int local_n = n/size;
    int resto = n%size;
    if(rank < resto)
        ++local_n;
    int local_dim = local_n*n;

    std::vector<int> displ_vector(size);
    std::vector<int> sendcount_vector(size);

    displ_vector[0] = 0;
    for(int i=1; i<size; ++i){
        int local_n = n/size;
        int resto = n%size;
        if(i-1 < resto)
            ++local_n;
        displ_vector[i] = n*local_n + displ_vector[i-1];
    }
    for(int i=0; i<size; ++i){
        int local_n = n/size;
        int resto = n%size;
        if(i < resto)
            ++local_n;
        sendcount_vector[i] = n*local_n;
    }

    Matrix U;

    if( rank == 0)
        U = Matrix::Zero(n, n);
 
    Matrix local_U(local_n + 2, n); // le local_U le facciamo più grandi di 2 righe così dopo è più comodo
    local_U.setZero();   // la prima e ultima riga sono zeri, le righe interne sono la U scatterata
   
    MPI_Scatterv(U.data(), sendcount_vector.data(), displ_vector.data(), MPI_DOUBLE, local_U.row(1).data(), local_dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    double h = 1.0/(n-1);
    Matrix old_local_U = local_U;
    double eps = 1;
    int check = 0;
    int local_check = 0;

    int top_neigbor, bottom_neigbor;
    if(rank == 0){
        top_neigbor = MPI_PROC_NULL;
        bottom_neigbor = 1;
    }
    if(rank == size-1){
        top_neigbor = size-2;
        bottom_neigbor = MPI_PROC_NULL;
    }
    if(rank != 0 and rank != size-1){
        top_neigbor = rank - 1;
        bottom_neigbor = rank + 1;
    }

    int iter = 0;
    while(check != 1 && iter < 100000){
        check = 0;
        local_check = 0;

        if(size>1){
        // step 1: comunicazione
            MPI_Sendrecv(local_U.row(1).data(), n, MPI_DOUBLE, top_neigbor, 0, old_local_U.row(0).data(), n, MPI_DOUBLE, top_neigbor, 
                        0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // invio la mia riga 1 a top, ricevo la riga ghost 0 da top
            MPI_Sendrecv(local_U.row(local_n).data(), n, MPI_DOUBLE, bottom_neigbor, 0, old_local_U.row(local_n+1).data(), n, MPI_DOUBLE, bottom_neigbor, 
                        0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // invio la mia riga local_n a bottom, ricevo la riga ghost local_n+1 da bottom
        }
        // step 2: calcolo
        #pragma omp parallel for collapse(2)
        for(int i = 1; i<=local_n; ++i){
            for(int j=1; j<n-1; ++j){
                double x = h*j;
                int global_i = (displ_vector[rank]/n + (i - 1));
                if(global_i == 0 || global_i == n - 1) 
                    continue;
                double y = h*global_i;
                local_U(i,j) = 0.25*(old_local_U(i-1,j) + old_local_U(i+1,j) + old_local_U(i,j-1) + old_local_U(i,j+1) + h*h*f(x,y));
            }
        }

        // step 3: check condizione
        eps = 0.0;
        #pragma omp parallel for collapse(2) reduction(+:eps)
        for(int i = 1; i<=local_n; ++i){
            for(int j=0; j<n; ++j){
                eps += (local_U(i,j)-old_local_U(i,j))*(local_U(i,j)-old_local_U(i,j));
            }
        }
        eps = h*eps;
        eps = std::sqrt(eps);

        if(eps<tol)
            local_check = 1;
        MPI_Allreduce(&local_check, &check, 1, MPI_INT, MPI_PROD, MPI_COMM_WORLD);

        old_local_U = local_U;
        iter++;
        }

    MPI_Gatherv(local_U.row(1).data(), local_dim, MPI_DOUBLE,
                    U.data(), sendcount_vector.data(), displ_vector.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

    return U;
}


void save_vtk_xml(const std::string& filename, const Matrix& u, int n) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Errore nell'apertura del file VTK!" << std::endl;
        return;
    }

    double h = 1.0/(n-1);

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
    // Output the matrix data
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





