#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
      
    //Variables for MPI initialization
    int rank, size, i, provided;
    
    //Initialize mpi
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    //Variables for the timer
    double start_time, stop_time, elapsed_time;
    start_time = MPI_Wtime();

    //Call useful functions
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //let each mpi process do their computations
    int count_loc = 0;
    int number_reductions=log(size)/log(2); 
    //printf("log is: %d\n",number_reductions);
    
    double x, y, z, pi;
    int NUM_ITER_LOC=NUM_ITER/size;
    
    srand(SEED*rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    
    // Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < NUM_ITER_LOC; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            count_loc++;
        }
    }
    
    
    // Now do the comunication
    
    // start by the reductions
    /*The number of reductions was calculated before. It is basically
    log_2(size). So for size=8=2^3 we have 3 reduction levels.
    */
    for (int j=1;j<=number_reductions;j++){
        /* Rank 0 will always be in the recieving side. for each reduction level j
        the ranks that will recieve information are those that meet the criteria that
        the modulus(rank,2^j)==0
        */ 
    	if (rank==0 || rank % (int)pow(2,j) ==0){
   		int count_temp=0;
    	
   	 	//Recieve data
   	 	/* The first reduction allways recieve data from the following neighbour
   	 	subsequent reductions recieve data from the closer following neighbour that
   	 	recieved data in the previous reduction.
   	 	*/ 
   	 	if (j==1){
    		MPI_Recv(&count_temp, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    		} else {
    		MPI_Recv(&count_temp, 1, MPI_INT, rank+pow(2,(j-1)), 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    			}
    		
       	//add the value of the current reduction
       	count_loc=count_loc+count_temp; 
       	
    	} else {
    	
    		//Send data
    		if (j==1){
    		MPI_Send(&count_loc, 1, MPI_INT, rank-1, 0, MPI_COMM_WORLD);
    		/* for reductions > 1, the ranks that send data are the ones
    		that could enter the if statement in the previous reduction
    		but couldn't enter in the current one.
    		*/
    		} else if (rank % (int)pow(2,j-1)==0) {
    		MPI_Send(&count_loc, 1, MPI_INT, rank-pow(2,(j-1)), 0, MPI_COMM_WORLD);
    		}
    	}
    }
    
    /* At this stage, all reductions are done. Here the local count of rank 0
    has already recieved the sums of ALL other ranks counts. So the value of pi
    can be calculated directly with count_loc
    */
    if (rank==0){
    // Estimate Pi and display the result
    	pi = ((double)count_loc / (double)NUM_ITER) * 4.0;
    	
    	stop_time = MPI_Wtime();
	elapsed_time = stop_time - start_time;
    	
    	printf("The result is %f. It took %f seconds\n", pi, elapsed_time);	
    }
    
    MPI_Finalize();
    
    return 0;
}

