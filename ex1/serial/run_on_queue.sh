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
## Run make in the src folder (modify properly)

module load openmp
cd /home/parallel/parlab10/ex1/serial
export OMP_NUM_THREADS=8
./gameOfLife 64 1000

