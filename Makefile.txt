CXX = mpicxx

INCLUDES = -I. -I$(PACS_ROOT)/include/eigen3 -I/usr/include/eigen3

CXXFLAGS = -Wall -O3 $(INCLUDES) -fopenmp

SRCS = main.cpp implementations.cpp
OBJS = $(SRCS:.cpp=.o)

all: main

main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o main

main.o: main.cpp declarations.hh
	$(CXX) $(CXXFLAGS) -c main.cpp

implementations.o: implementations.cpp declarations.hh
	$(CXX) $(CXXFLAGS) -c implementations.cpp

run: main
	mpirun -np 4 ./main 100 1e-5

clean:
	rm -f *.o main poisson_mpi_solution.vti

.PHONY: all run clean