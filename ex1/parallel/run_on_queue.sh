#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_gameOfLife

## Output and error files
#PBS -o run_gameOfLife.out
#PBS -e run_gameOfLife.err

## How many machines should we get? 
#PBS -l nodes=1:ppn=8

##How long should the job run for?
#PBS -l walltime=00:10:00

## Start 
## Run make in the src folder 

module load openmp
cd /home/parallel/parlab10/ex1/parallel

for thr in 1 2 4 6 8
do
	export OMP_NUM_THREADS=$thr
	./gameOfLife 64 1000
	./gameOfLife 1024 1000
	./gameOfLife 4096 1000
done

