/*Calculate Pi using MPI Blocking Communication & Linear Reduction Algorithm */
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

    //Call useful functions
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //let each mpi process do their computations
    int count_loc = 0;
    
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
    start_time = MPI_Wtime();

    if (rank==0){
    	int count_glob=0;
    	int count[size-1];
    	
    	//Recieve data from all other ranks
    	for (int i=1;i<size;i++){
    	MPI_Recv(&count[i-1], 1, MPI_INT, i, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    	}
    	
       // Perform the calculation of pi
       //add the value of rank 0
       count_glob=count_loc; 
       // add the values of the other ranks	
    	for (int i=1;i<size;i++){
    	count_glob=count_glob+count[i-1];
    	}
    	
    	// Estimate Pi and display the result
    	pi = ((double)count_glob / (double)NUM_ITER) * 4.0;
    	
    	stop_time = MPI_Wtime();
	elapsed_time = stop_time - start_time;
    	
    	printf("The result is %f. It took %f seconds\n", pi, elapsed_time);	
    	
    } else { 
    
    //Send data to rank 0
    MPI_Send(&count_loc, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    
    
    MPI_Finalize();
    
    return 0;
}
