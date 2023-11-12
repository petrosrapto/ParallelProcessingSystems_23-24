#!/bin/bash

## Give the Job a descriptive name
#PBS -N run_kmeans

## Output and error files
#PBS -o run_kmeans.out
#PBS -e run_kmeans.err

## How many machines should we get? 
#PBS -l nodes=1:ppn=8

##How long should the job run for?
#PBS -l walltime=00:10:00

## Start 
## Run make in the src folder (modify properly)

module load openmp
cd /home/parallel/parlab10/ex2_p/kmeans/build

# the name of the executable is given as the environmental
# variable "EXECUTABLE".
if [ -z "$EXECUTABLE" ]; then
    echo "No executable name provided. Submit to torque as follows: qsub -q parlab -v EXECUTABLE=seq_kmeans run_on_queue.sh">&2
    exit 1
fi

./kmeans_seq -s 256 -n 16 -c 16 -l 10

for thr in 1 2 4 8 16 32 64
do
	# the file run_kmeans.out is specified as stdout 
	export OMP_NUM_THREADS=$thr
	./$EXECUTABLE -s 256 -n 16 -c 16 -l 10
done

#./kmeans_seq -s <SIZE> -n <COORDS> -c <CLUSTERS> -l <LOOPS>
