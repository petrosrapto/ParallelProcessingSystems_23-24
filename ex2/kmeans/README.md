The makefile takes the environmental variable SRC_FILE as argument, which file to "make".
So, the make_on_queue.sh should be executed providing the proper file as the SRC_FILE var.
Similarly, the run_on_queue.sh takes the environmental variable EXECUTABLE as argument and runs that exec with different thread number and the given configuration.
Examples submitting the above to sandman: 
$ qsub -q serial -l nodes=sandman:ppn=64 -v SRC_FILE=seq_kmeans make_on_queue.sh
$ qsub -q serial -l nodes=sandman:ppn=64 -v EXECUTABLE=omp_naive_kmeans run_on_queue.sh
