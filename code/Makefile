#Compiler chosen
CC=gcc
#MPI Compiler
MPICC=mpicc
#CFLAGS. To be expirimented with later
CFLAGS=-O3 
#LFLAGS. Ising requires GNU Scientific Library and libmath
LFLAGS=-lm -lgsl -lgslcblas -openmp -pg
#PFLAGS. Flags for parallel computing
################################################################
#TODO: implement make check for unit tests
.PHONY: all clean

all: ising

ising: ising.c
	$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

clean: 
	rm -r -f ising
	
