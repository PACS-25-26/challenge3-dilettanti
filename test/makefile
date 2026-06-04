CXX = mpicxx

# -I../header permette al compilatore di trovare sia i tuoi file .hh 
# sia la cartella Eigen (risolvendo <Eigen/Dense> come ../header/Eigen/Dense)
CXXFLAGS = -Wall -O3 -I../header

# I sorgenti e gli oggetti
SRCS = main.cpp ../src/implementations.cpp
OBJS = main.o implementations.o

all: main

main: $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o main

# Compilazione di main.o
main.o: main.cpp ../header/declarations.hh
	$(CXX) $(CXXFLAGS) -c main.cpp -o main.o

# Compilazione di implementations.o
implementations.o: ../src/implementations.cpp ../header/declarations.hh
	$(CXX) $(CXXFLAGS) -c ../src/implementations.cpp -o implementations.o

run: main
	mpirun -np 4 ./main 

clean:
	rm -f *.o main .vti

.PHONY: all run clean