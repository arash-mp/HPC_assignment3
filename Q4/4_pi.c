
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
    
    // Do the communications
    // Each MPI process sends its rank to reduction, root MPI process collects the result
    int reduction_result = 0;
    int root_rank=0;
    MPI_Reduce(&count_loc, &reduction_result, 1, MPI_INT, MPI_SUM, root_rank, MPI_COMM_WORLD);
    
    if (rank==root_rank){
    	
    	// Estimate Pi and display the result
    	pi = ((double)reduction_result / (double)NUM_ITER) * 4.0;
    	
    	stop_time = MPI_Wtime();
	elapsed_time = stop_time - start_time;
    	
    	printf("The result is %f. It took %f seconds\n", pi, elapsed_time);	
    	}
    
    
    MPI_Finalize();
    
    return 0;
}


