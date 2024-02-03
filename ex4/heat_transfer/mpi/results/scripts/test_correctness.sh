#!/bin/bash

## Give the Job a descriptive name
#PBS -N testcorrectness

## Output and error files
#PBS -o test_correctness.out
#PBS -e test_correctness.err

## How many nodes:processors_per_node should we get?
## Run on parlab
#PBS -l nodes=8:ppn=8
 
## NOTE: Fix the path to show to your serial executables 

module load openmpi/1.8.3
cd /home/parallel/parlab10/ex4_p/heat_transfer/serial

for execfile in jacobi seidelsor
##jacobi seidelsor redblacksor
do
	./${execfile} 256 256
done

## NOTE: Fix the path to show to your MPI executables
cd /home/parallel/parlab10/ex4_p/heat_transfer/mpi

for execfile in jacobi seidelsor
## jacobi seidelsor redblacksor
do
	mpirun -np 32 --map-by node --mca btl self,tcp ${execfile} 256 256 8 4
done

## Make sure you enable convergence testing and printing
