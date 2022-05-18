#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char *argv[])
{

	// Define the variables
	int rank, size, provided;

	MPI_Init_thread(&argc,&argv,MPI_THREAD_SINGLE,&provided);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	printf("Hello World from rank %d out of %d \n", rank,size);
	
	MPI_Finalize();
	
	
		
   return 0; 
}
