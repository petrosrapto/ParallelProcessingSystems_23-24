CC=mpicc
CFLAGS=-Wall -Wextra -Wno-unused -O3
# uncomment the following lines in order not to define the following macros

all: jacobi seidelsor redblacksor

jacobi: Jacobi_mpi.c utils.c
	$(CC) $(CFLAGS) Jacobi_mpi.c utils.c -o jacobi

seidelsor: GaussSeidelSOR_mpi.c utils.c
	$(CC) $(CFLAGS) GaussSeidelSOR_mpi.c utils.c -o seidelsor -lm -lmpi

redblacksor: RedBlackSOR_mpi.c utils.c
	$(CC) $(CFLAGS) RedBlackSOR_mpi.c utils.c -o redblacksor -lm -lmpi

clean:
	rm jacobi seidelsor redblacksor


