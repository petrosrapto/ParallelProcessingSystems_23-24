#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "mpi.h"
#include "utils.h"

//#define TEST_CONV
//#define DEBUGGING
//#define PRINT_RESULTS

void RedSOR(double ** u_previous, double ** u_current, int X_min, int X_max, int Y_min, int Y_max, double omega) {
	int i,j;
	for (i=X_min;i<X_max;i++)
		for (j=Y_min;j<Y_max;j++)
			if ((i+j)%2==0) 
				u_current[i][j]=u_previous[i][j]+(omega/4.0)*(u_previous[i-1][j]+u_previous[i+1][j]+u_previous[i][j-1]+u_previous[i][j+1]-4*u_previous[i][j]);		         
}

void BlackSOR(double ** u_previous, double ** u_current, int X_min, int X_max, int Y_min, int Y_max, double omega) {
	int i,j;
	for (i=X_min;i<X_max;i++)
		for (j=Y_min;j<Y_max;j++)
			if ((i+j)%2==1) 
				u_current[i][j]=u_previous[i][j]+(omega/4.0)*(u_current[i-1][j]+u_current[i+1][j]+u_current[i][j-1]+u_current[i][j+1]-4*u_previous[i][j]); 
}


int main(int argc, char ** argv) {
    int rank,size;
    int global[2],local[2]; //global matrix dimensions and local matrix dimensions (2D-domain, 2D-subdomain)
    int global_padded[2];   //padded global matrix dimensions (if padding is not needed, global_padded=global)
    int grid[2];            //processor grid dimensions
// how many processes there are in each dimension. if grid[0]=16 then 16 processes share the 0th
// dimension. if the global matrix has 'x' cells in the 0th dimension, then each process will get
// 'x'/16 cells in that dimension
    int i,j,t;
    int global_converged=0,converged=0; //flags for convergence, global and per process
    MPI_Datatype dummy;     //dummy datatype used to align user-defined datatypes in memory
    double omega;                       //relaxation factor - useless for Jacobi

    struct timeval tts,ttf,tcs,tcf, tcvs, tcvf;   //Timers: total-> tts,ttf, computation -> tcs,tcf
    double ttotal=0,tcomp=0,tconv=0,total_time,comp_time,conv_time;

    double ** U, ** u_current, ** u_previous, ** swap; //Global matrix, local current and previous matrices, pointer to swap between current and previous
    

    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
/*
    #ifdef DEBUGGING
    printf("Rank %d out of %d\n", rank, size);
    fflush(stdout);
    #endif
*/

    //----Read 2D-domain dimensions and process grid dimensions from stdin----//

    if (argc!=5) {
        fprintf(stderr,"Usage: mpirun .... ./exec X Y Px Py");
        exit(-1);
    }
    else {
        global[0]=atoi(argv[1]);
        global[1]=atoi(argv[2]);
        grid[0]=atoi(argv[3]);
        grid[1]=atoi(argv[4]);
    }


    //----Create 2D-cartesian communicator----//
    //----Usage of the cartesian communicator is optional----//

    MPI_Comm CART_COMM;         //CART_COMM: the new 2D-cartesian communicator
    int periods[2]={0,0};       //periods={0,0}: the 2D-grid is non-periodic
    int rank_grid[2];           //rank_grid: the position of each process on the new communicator

    MPI_Cart_create(MPI_COMM_WORLD,2,grid,periods,0,&CART_COMM);    //communicator creation
    MPI_Cart_coords(CART_COMM,rank,2,rank_grid);                  //rank mapping on the new communicator
/*
    #ifdef DEBUGGING
    printf("Rank %d: Cartesian coordinates [%d, %d]\n", rank, rank_grid[0], rank_grid[1]);
    fflush(stdout);
    #endif
*/


    //----Compute local 2D-subdomain dimensions----//
    //----Test if the 2D-domain can be equally distributed to all processes----//
    //----If not, pad 2D-domain----//

    for (i=0;i<2;i++) {
        if (global[i]%grid[i]==0) {
            local[i]=global[i]/grid[i];
            global_padded[i]=global[i];
        }
        else {
// add one row/column in order to cover all cells of the global matrix
            local[i]=(global[i]/grid[i])+1;
            global_padded[i]=local[i]*grid[i];
        }
    }

        //Initialization of omega
    omega=2.0/(1+sin(3.14/global[0]));


    //----Allocate global 2D-domain and initialize boundary values----//
    //----Rank 0 holds the global 2D-domain----//
    if (rank==0) {
        U=allocate2d(global_padded[0],global_padded[1]);
        init2d(U,global[0],global[1]);
    }

    //----Allocate local 2D-subdomains u_current, u_previous----//
    //----Add a row/column on each size for ghost cells----//

    u_previous=allocate2d(local[0]+2,local[1]+2);
    u_current=allocate2d(local[0]+2,local[1]+2);

    if (u_previous == NULL || u_current == NULL) {
        fprintf(stderr, "Rank %d: Memory allocation failed.\n", rank);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    //----Distribute global 2D-domain from rank 0 to all processes----//

    //----Appropriate datatypes are defined here----//
    /*****The usage of datatypes is optional*****/

    //----Datatype definition for the 2D-subdomain on the global matrix----//

    MPI_Datatype global_block;
    MPI_Type_vector(local[0],local[1],global_padded[1],MPI_DOUBLE,&dummy);
    MPI_Type_create_resized(dummy,0,sizeof(double),&global_block);
    MPI_Type_commit(&global_block);

    //----Datatype definition for the 2D-subdomain on the local matrix----//

    MPI_Datatype local_block;
    MPI_Type_vector(local[0],local[1],local[1]+2,MPI_DOUBLE,&dummy);
    MPI_Type_create_resized(dummy,0,sizeof(double),&local_block);
    MPI_Type_commit(&local_block);

/*
Defining both global_block and local_block datatypes allows for more precise control over the memory layout and data distribution in a distributed computing environment. It ensures that both the scattering (distribution) and gathering (collection) of data among multiple

processes are handled correctly, accounting for any differences in memory layout between the global matrix and the local submatrices with potential ghost cells. */


    //----Rank 0 defines positions and counts of local blocks (2D-subdomains) on global matrix----//
    int * scatteroffset, * scattercounts;
    if (rank==0) {
        scatteroffset=(int*)malloc(size*sizeof(int));
        scattercounts=(int*)malloc(size*sizeof(int));
        for (i=0;i<grid[0];i++)
            for (j=0;j<grid[1];j++) {
                scattercounts[i*grid[1]+j]=1;
                scatteroffset[i*grid[1]+j]=(local[0]*local[1]*grid[1]*i+local[1]*j);
            }
    }


    //----Rank 0 scatters the global matrix----//

    //----Rank 0 scatters the global matrix----//

        //*************TODO*******************//

        /*Fill your code here*/

        // we use MPI_COMM_WORLD and not CART_COMM cause we need to broadcast to all processes and
        // we dont use the grid layout for the operation

        // in C &(U[0][0]) and U are NOT the same.
        // U can be thought of as a pointer to the start of a "list of rows."
        // &(U[0][0]) is a pointer to the first element of the first row in the 2D array.


    MPI_Scatterv(&(U[0][0]), scattercounts, scatteroffset, global_block, &(u_previous[1][1]), 1, local_block, 0, MPI_COMM_WORLD);
    MPI_Scatterv(&(U[0][0]), scattercounts, scatteroffset, global_block, &(u_current[1][1]), 1, local_block, 0, MPI_COMM_WORLD);

    /* Make sure u_current and u_previous are both initialized */

    //************************************//


    if (rank==0)
        free2d(U);


    //----Define datatypes or allocate buffers for message passing----//

        //*************TODO*******************//



        /*Fill your code here*/


    /* Define a datatype that corresponds to a collumn of a local block */
    MPI_Datatype column;
    MPI_Type_vector(local[0], 1,local[1]+2,MPI_DOUBLE,&dummy);
    MPI_Type_create_resized(dummy,0,sizeof(double),&column);
    MPI_Type_commit(&column);


    //************************************//


    //----Find the 4 neighbors with which a process exchanges messages----//
       
	//*************TODO*******************//
    int north, south, east, west;

	/*Fill your code here*/
    MPI_Cart_shift(CART_COMM, 0, 1, &north, &south);
    MPI_Cart_shift(CART_COMM, 1, 1, &west, &east);


    /*Make sure you handle non-existing neighbors appropriately*/

    /*
MPI_Cart_shift(CART_COMM, 0, 1, &north, &south);
This function call is used to determine the ranks of the neighboring processes along the first dimension (dimension 0) of the Cartesian grid.
The second argument (0) specifies the dimension along which the neighbors are to be found. In this case, it is the first dimension (which could be thought of as the "row" if you visualize the grid).
The third argument (1) is the displacement. It indicates how far from the current process the neighbors are to be found. A displacement of 1 means the immediate neighbors.
&north and &south are pointers to integers where the ranks of the northern and southern neighbors will be stored. After this call, north will contain the rank of the process directly above the current process in the grid, and south will contain the rank of the process directly below.
*/
/*
if there isnt a neighbor to a specific direction, that variable will get value 'MPI_PROC_NULL'
*/


    //************************************//

    //---Define the iteration ranges per process-----//
        //*************TODO*******************//

    int i_min, i_max, j_min, j_max;


        /*Fill your code here*/


    /* internal process (ghost cell only) */
    i_min = 1;
    i_max = local[0] + 1;

    /* boundary process - no possible padding */
    if (north == MPI_PROC_NULL) {
        i_min = 2;  // ghost cell + boundary
    }

    /* boundary process and padded global array */
    if (south == MPI_PROC_NULL){
        i_max -= (global_padded[0] - global[0]) + 1;
    }

    /* internal process (ghost cell only) */
    j_min = 1;
    j_max = local[1] + 1;

    /* boundary process - no possible padding */
    if (west == MPI_PROC_NULL) {
        j_min = 2;  //ghost cell + boundary
    }

    /* boundary process and padded global array */
    if (east == MPI_PROC_NULL){
        j_max -= (global_padded[1] - global[1]) + 1;
    }



    /*Three types of ranges:
      -internal processes
      -boundary processes
      -boundary processes and padded global array
     */
/*
-Internal Process (Only Ghost Cells)
        i_min = 1 and i_max = local[0] + 1:
For an internal process (one that is not on the boundary of the grid), the iteration starts from 1 and
goes up to local[0] + 1. The +1 is to account for the extra row/column of ghost cells. This means that the actual computation starts one row/column inside of the local domain to avoid modifying the ghost cells, which are used for communication with neighboring processes.

-Boundary Process (No Padding)
        if (north == MPI_PROC_NULL) { i_min = 2; }:
If the process has no northern neighbor (i.e., it's on the north boundary), the iteration in the i-direction should start from 2 instead of 1. This adjustment skips the first row, which is a ghost cell row, and also skips the actual boundary row because it's a fixed boundary that shouldn't be updated. U is one-based through X-1.
        if (west == MPI_PROC_NULL) { j_min = 2; }:
Similarly, if the process is on the west boundary, it starts iterating from the second column, skipping the ghost cell column and the boundary column.
-Boundary Process (With Padding)
        if (south == MPI_PROC_NULL) { i_max -= (global_padded[0] - global[0]) + 1; }:
If the process is on the south boundary, i_max is reduced to exclude the padded rows (if any). The subtraction (global_padded[0] - global[0]) calculates the number of padded rows, and +1 adjusts for the boundary row that should not be updated.
        if (east == MPI_PROC_NULL) { j_max -= (global_padded[1] - global[1]) + 1; }:
For a process on the east boundary, the iteration in the j-direction is adjusted similarly to exclude padded columns and the boundary column.
*/


    //************************************//


    //----Computational core----//
    gettimeofday(&tts, NULL);
#ifdef TEST_CONV
    for (t=0;t<T && !global_converged;t++) {
#endif
#ifndef TEST_CONV
#undef T
#define T 256
    for (t=0;t<T;t++) {
#endif


                //*************TODO*******************//


        /*Fill your code here*/
        swap=u_previous;
        u_previous=u_current;
        u_current=swap;

        /*Fill your code here*/

        /*Compute and Communicate*/

	// communication for the red phase
        if (north != MPI_PROC_NULL){
            MPI_Sendrecv(&u_previous[1][1], local[1], MPI_DOUBLE, north, 0, &u_previous[0][1], local[1], MPI_DOUBLE, north, 0, MPI_COMM_WORLD, &status );
        }

        // Communicate with south
        if (south != MPI_PROC_NULL){
            MPI_Sendrecv(&u_previous[local[0]][1], local[1], MPI_DOUBLE, south, 0, &u_previous[local[0]+1][1], local[1], MPI_DOUBLE, south, 0, MPI_COMM_WORLD, &status );
        }

        // Communicate with east
        if (east != MPI_PROC_NULL){
            MPI_Sendrecv(&u_previous[1][local[1]], 1, column, east, 0, &u_previous[1][local[1]+1], 1, column, east, 0, MPI_COMM_WORLD, &status );
        }

        // Communicate with west
        if (west != MPI_PROC_NULL){
            MPI_Sendrecv(&u_previous[1][1], 1, column, west, 0, &u_previous[1][0], 1, column, west, 0, MPI_COMM_WORLD, &status );
        }

        /*Add appropriate timers for computation*/
        gettimeofday(&tcs, NULL);

       /*#ifdef DEBUGGING
        printf("Rank %d: Starting Jacobi iteration %d with ranges X[%d, %d), Y[%d, %d)\n", rank, t, i_min, i_max, j_min, j_max);
        #endif*/

	RedSOR(u_previous, u_current, i_min, i_max, j_min, j_max, omega);

        gettimeofday(&tcf, NULL);

        tcomp += (tcf.tv_sec-tcs.tv_sec)+(tcf.tv_usec-tcs.tv_usec)*0.000001;


	// communication for the black phase
        // Communicate with north
        if (north != MPI_PROC_NULL){
            MPI_Sendrecv(&u_current[1][1], local[1], MPI_DOUBLE, north, 0, &u_current[0][1], local[1], MPI_DOUBLE, north, 0, MPI_COMM_WORLD, &status );
        }

        // Communicate with south
        if (south != MPI_PROC_NULL){
            MPI_Sendrecv(&u_current[local[0]][1], local[1], MPI_DOUBLE, south, 0, &u_current[local[0]+1][1], local[1], MPI_DOUBLE, south, 0, MPI_COMM_WORLD, &status );
        }

        // Communicate with east
        if (east != MPI_PROC_NULL){
            MPI_Sendrecv(&u_current[1][local[1]], 1, column, east, 0, &u_current[1][local[1]+1], 1, column, east, 0, MPI_COMM_WORLD, &status );
        }

        // Communicate with west
        if (west != MPI_PROC_NULL){
            MPI_Sendrecv(&u_current[1][1], 1, column, west, 0, &u_current[1][0], 1, column, west, 0, MPI_COMM_WORLD, &status );
        }

        gettimeofday(&tcs, NULL);

        BlackSOR(u_previous, u_current, i_min, i_max, j_min, j_max, omega);

        gettimeofday(&tcf, NULL);

        tcomp += (tcf.tv_sec-tcs.tv_sec)+(tcf.tv_usec-tcs.tv_usec)*0.000001;


#ifdef TEST_CONV
        if (t%C==0) {
            //*************TODO**************//

            /*Test convergence*/
	    
	    // i_max - 1 cause in the utils.c the loops of the converge() goes up to <=i_max,j_max
                // while the loops of the Jacobi() goes up to < i_max,j_max

	   
	    gettimeofday(&tcvs, NULL);

            converged=converge(u_previous,u_current,i_min , i_max - 1 , j_min , j_max - 1);

            MPI_Allreduce(&converged, &global_converged, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	    gettimeofday(&tcvf, NULL);
            tconv += (tcvf.tv_sec-tcvs.tv_sec)+(tcvf.tv_usec-tcvs.tv_usec)*0.000001;
        }
#endif

    }
    gettimeofday(&ttf,NULL);

    ttotal=(ttf.tv_sec-tts.tv_sec)+(ttf.tv_usec-tts.tv_usec)*0.000001;

    MPI_Reduce(&ttotal,&total_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&tcomp,&comp_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&tconv,&conv_time,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

    //----Rank 0 gathers local matrices back to the global matrix----//

    if (rank==0) {
        U=allocate2d(global_padded[0],global_padded[1]);
    }

    //*************TODO*******************//


    /*Fill your code here*/

    MPI_Gatherv(&u_current[1][1], 1, local_block, &(U[0][0]), scattercounts, scatteroffset, global_block, 0, MPI_COMM_WORLD);

   //************************************//


   //----Printing results----//

        //**************TODO: Change "Jacobi" to "GaussSeidelSOR" or "RedBlackSOR" for appropriate printing****************//


    if (rank==0) {
        printf("RedBlackSOR X %d Y %d Px %d Py %d Iter %d ComputationTime %lf TotalTime %lf Convergence Time %lf  midpoint %lf\n",global[0],global[1],grid[0],grid[1],t,comp_time,total_time,conv_time, U[global[0]/2][global[1]/2]);

#ifdef PRINT_RESULTS
        char * s=malloc(50*sizeof(char));
        sprintf(s,"resJacobiMPI_%dx%d_%dx%d",global[0],global[1],grid[0],grid[1]);
        fprint2d(s,U,global[0],global[1]);
        free(s);
#endif

    }
    MPI_Finalize();
    return 0;
}
