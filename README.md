[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/tKSbaXxd)
# Challenge 3: Parallel Solver for the Laplace Equation

This repository contains the implementation of a hybrid parallel solver (MPI + OpenMP) for the Laplace equation on a square domain $\Omega = (0,1)^2$ with homogeneous Dirichlet boundary conditions, solved using the Jacobi iteration method.

## Project Structure

The code is divided into the following main folders:
- `test`: Contains `main.cpp` that execute the time benchmarks and print the results. It also contains `chrono.hpp` used for computational time measurement.
- `header.hh`: Contains file  `declarations.hh` and the `Eigen` library.
- `src`: Contains file `implementations.cpp` with the source code of the implementation.


## Compilation

The project includes a `Makefile` to automate the process. To compile the code inside the `test` folder run:

```bash
make
```

To remove object files and the generated executable:

```bash
make clean
```

## Execution

The program can be launched by specifying the number of MPI parallel tasks using `mpirun`. 

For example, to run the solver with **4 MPI processes**:
```bash
mpirun -np 4 ./main
```

### Hybrid Parallelism Management (MPI + OpenMP)
To fully exploit the hybrid architecture, you can configure the number of OpenMP threads assigned to each MPI process by setting the `OMP_NUM_THREADS` environment variable:

```bash
export OMP_NUM_THREADS=2
mpirun -np 4 ./main
```

## Output and Visualization

At the end of the computation, the program exports the grid and the numerical solution in XML ImageData format (`.vti`):
- `sequential_solution.vti`
- `parallel_solution.vti`

The generated files are ready to be imported and visually analyzed within scientific post-processing tools like **ParaView**.