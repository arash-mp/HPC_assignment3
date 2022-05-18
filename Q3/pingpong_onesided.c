#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	/* -------------------------------------------------------------------------------------------
		MPI Initialization 
	--------------------------------------------------------------------------------------------*/
	MPI_Init(&argc, &argv);

	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Status stat;

	if(size != 2){
		if(rank == 0){
			printf("This program requires exactly 2 MPI ranks, but you are attempting to use %d! Exiting...\n", size);
		}
		MPI_Finalize();
		exit(0);
	}

	/* -------------------------------------------------------------------------------------------
		Loop from 8 B to 1 GB
	--------------------------------------------------------------------------------------------*/

	for(int i=0; i<=27; i++){

		long int N = 1 << i;
	
   	 	// Allocate memory for A on CPU
		double *A = (double*)malloc(N*sizeof(double));
		
		// Allocate memory for R on CPU to recieve the data
		double *R = (double*)malloc(N*sizeof(double));
		

		// Initialize all elements of A to 0.0
		for(int i=0; i<N; i++){
			A[i] = 0.0;
		}
	
		// create the window to be exposed (each rank will expose A)
		MPI_Win win;
		MPI_Win_create(A,N*sizeof(double),sizeof(double),MPI_INFO_NULL,MPI_COMM_WORLD, &win);
	
		int tag1 = 10;
		int tag2 = 20;
		
		// Set up from which rank every other rank needs to recieve info
		int targetrank;
		if (rank==0){ // This way each origin will recieve from another target
			targetrank=1;
		} else {
			targetrank=0;
		}
		
			
		int loop_count = 50;

		// Warm-up loop
		for(int i=1; i<=5; i++){
			
			MPI_Win_fence(0, win);
    			MPI_Get(R, N, MPI_DOUBLE, targetrank,0, N, MPI_DOUBLE, win);
    			//Sync to make sure the get is complete
    			MPI_Win_fence(0, win);
		}

		// Time ping-pong for loop_count iterations of data transfer size 8*N bytes
		double start_time, stop_time, elapsed_time;
		start_time = MPI_Wtime();
	
		for(int i=1; i<=loop_count; i++){
			
			//Open the communication epoch in window win
			MPI_Win_fence(0, win);
			
			// here each rank recieve into array R from target rank set up before
			// rank 0 is recieving A into R from rank 1
			// rank 1 is recieving A into R from rank 0
			// These trasnfers will happen simultaneously
    			MPI_Get(R, N, MPI_DOUBLE, targetrank,0, N, MPI_DOUBLE, win);
    			
    			//Sync to make sure the get is complete
    			MPI_Win_fence(0, win);
		}

		stop_time = MPI_Wtime();
		elapsed_time = stop_time - start_time;

		long int num_B = 8*N;
		long int B_in_GB = 1 << 30;
		double num_GB = (double)num_B / (double)B_in_GB;
		//to calculate the time, do not devide by 2 as both messages happen at the same time. 
		double avg_time_per_transfer = elapsed_time / ((double)loop_count);

		if(rank == 0) printf("%10li\t%15.9f\n", num_B, avg_time_per_transfer);

		free(A);
		free(R);
		MPI_Win_free(&win);
		
		
		/* A conclussion from this is: when we use one sided communication,
		we can recieve data at the same time from both ranks. This allows to have
		a total communication time that is lower. But apparently, the communication 
		of a single mesagge is less efficient, as the system needs to share the bandwith 
		between all the messages happening at the same time.
		*/
		
	}

	MPI_Finalize();

	return 0;
}
