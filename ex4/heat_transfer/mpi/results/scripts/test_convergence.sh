#!/bin/bash

## Give the Job a descriptive name
#PBS -N testconvergence

## Output and error files
#PBS -o test_convergence_4_16_3rd.out
#PBS -e test_convergence_4_16_3rd.err

## How many machines should we get? 
#PBS -l nodes=8:ppn=8

## Start 

module load openmpi/1.8.3
cd /home/parallel/parlab10/ex4_p/heat_transfer/mpi

mpirun -np 64 --mca btl self,tcp ./jacobi 1024 1024 4 16

## make sure that you have enabled convergence
