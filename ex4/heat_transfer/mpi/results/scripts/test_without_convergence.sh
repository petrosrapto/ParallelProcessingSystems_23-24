#!/bin/bash

## Give the Job a descriptive name
#PBS -N test_without_convergence

## Output and error files
#PBS -o test_without_convergence.out
#PBS -e test_without_convergence.err

## How many nodes:processors_per_node should we get?
## Run on parlab
#PBS -l nodes=8:ppn=8

## Run the job (use full paths to make sure we execute the correct thing) 
## NOTE: Fix the path to show to your executable! 

module load openmpi/1.8.3
cd /home/parallel/parlab10/ex4_p/heat_transfer/mpi

for size in 2048 4096 6144
do
	for execfile in jacobi ##seidelsor redblacksor
	do
		mpirun  -np 1 --map-by node --mca btl self,tcp ${execfile} ${size} ${size} 1 1 >> ${execfile}_s${size}_c1
		mpirun  -np 2 --map-by node --mca btl self,tcp ${execfile} ${size} ${size} 2 1 >> ${execfile}_s${size}_c2
		mpirun  -np 4 --map-by node --mca btl self,tcp ${execfile} ${size} ${size} 2 2 >> ${execfile}_s${size}_c4
		mpirun  -np 8 --map-by node --mca btl self,tcp ${execfile} ${size} ${size} 4 2 >> ${execfile}_s${size}_c8
		mpirun  -np 16 --map-by node --mca btl self,tcp ${execfile} ${size} ${size} 4 4 >> ${execfile}_s${size}_c16
		mpirun  -np 32 --map-by node --mca btl self,tcp ${execfile} ${size} ${size} 8 4 >> ${execfile}_s${size}_c32
                mpirun  -np 64 --map-by node --mca btl self,tcp ${execfile} ${size} ${size} 8 8 >> ${execfile}_s${size}_c64
	done
done

## Make sure you disable convergence testing and printing
