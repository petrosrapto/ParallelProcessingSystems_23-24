#!/bin/bash

## Give the Job a descriptive name
#PBS -N testing

## Output and error files
#PBS -o testing.out
#PBS -e testing.err

## How many nodes:processors_per_node should we get?
## Run on parlab
#PBS -l nodes=8:ppn=8

## Run the job (use full paths to make sure we execute the correct thing) 
## NOTE: Fix the path to show to your executable! 

module load openmpi/1.8.3
cd /home/parallel/parlab10/ex4_p/heat_transfer/mpi

mpirun -np 64 --map-by node --mca btl self,tcp jacobi 1024 1024 16 4

## Make sure you disable convergence testing and printing

