#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_mpi_helloworld

## Output and error files
#PBS -o run_mpi_helloworld.out
#PBS -e run_mpi_helloworld.err

## How many machines should we get? 
#PBS -l nodes=8:ppn=8

## Start 
## Run make in the src folder (modify properly)

module load openmpi/1.8.3
cd /home/parallel/parlab10/ex4_p/kmeans

processes='1 2 4 8 16 32 64'
for proc in $processes; do
	## mpirun -np $proc --map-by node --mca btl self,tcp ./kmeans_mpi -s 256 -n 16 -c 16 -l 10
	mpirun -np $proc --map-by node --mca btl self,tcp ./kmeans_mpi -s 256 -n 16 -c 16 -l 10
done


