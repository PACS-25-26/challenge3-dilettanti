CXX = mpicxx

CXXFLAGS = -Wall -O3 -I.

SRCS = main.cpp implementations.cpp
OBJS = main.o implementations.o

all: main

main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o main

main.o: main.cpp declarations.hh
	$(CXX) $(CXXFLAGS) -c main.cpp

implementations.o: implementations.cpp declarations.hh
	$(CXX) $(CXXFLAGS) -c implementations.cpp

run: main
	mpirun -np 4 ./main 

clean:
	rm -f *.o main poisson_mpi_solution.vti

.PHONY: all run clean