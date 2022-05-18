#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdbool.h>
 

// Define function to check the number of MPI ranks is perfect square or not.
int isPerfectSquare(int number)
{
    int iVar;
    float fVar;
 
    fVar=sqrt((double)number);
    iVar=fVar;
 
    if(iVar==fVar)
        return 1;
    else
        return 0;
}



#define ITER_OUTER_LOOP 5 
#define Matrix_Size 1024

// Define original matrix a*b=c
double matrix_a[Matrix_Size][Matrix_Size];
double matrix_b[Matrix_Size][Matrix_Size];
double matrix_c[Matrix_Size][Matrix_Size];

// Create a function to allocate 2d arrays
double **alloc_2d_array(int rows, int cols) {
  int i;
  double *data = (double *)malloc(rows*cols*sizeof(double));
  double **array= (double **)malloc(rows*sizeof(double*));
  for (i=0; i<rows; i++)
    array[i] = &(data[cols*i]);

  return array;
}

int main(int argc, char* argv[]){
  
    //Variables for MPI initialization
    int rank, size, provided;
    
    int i, j, k ;
    int rankprint=0;
    // Fill the array swith rand variables
    for (i = 0 ; i < Matrix_Size ; i++) {
      for (j = 0 ; j < Matrix_Size ; j++) {
        matrix_a[i][j] = (double) rand() / RAND_MAX;
        matrix_b[i][j] = (double) rand() / RAND_MAX;
        matrix_c[i][j] = 0.0;
      }
    } 
    
    // =====================Set Up Communications============================
        
    //Initialize mpi
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    //Call useful functions
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Check if number of MPI ranks is perfect square or not
    if (!isPerfectSquare(size)){
       if (rank==0){
       printf("Number of processes must be a perfect square \n");

       exit(0);
       }
        }

    
    // Ask MPI to decompose our processes in a 2D cartesian grid
    int dims[2] = {0, 0}; 
    MPI_Dims_create(size, 2, dims);
 
    // Make both dimensions periodic
    int periods[2] = {true, true}; // Periodic dimensions are neighbours.
    
    // Let MPI assign arbitrary ranks if it deems it necessary
    int reorder = false; 
 
    // Create a communicator given the 2D torus topology.
    MPI_Comm cartesian_communicator;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartesian_communicator);
 
    // My rank in the cartesian communicator communicator
    int cart_rank;
    MPI_Comm_rank(cartesian_communicator, &cart_rank);
 
    // Get my coordinates in the cartesian communicator
    int cart_coords[2];
    MPI_Cart_coords(cartesian_communicator, cart_rank, 2, cart_coords);
 
 
    
    // Find the neighbours of each rank
    // Declare neighbours 
    // From rookiehpc.com	
    enum DIRECTIONS {DOWN, UP, LEFT, RIGHT};
    char* neighbours_names[4] = {"down", "up", "left", "right"};
    int neighbours_ranks[4];
 
    // shift tells us our up and down neighbours
    MPI_Cart_shift(cartesian_communicator, 0, 1, &neighbours_ranks[UP], &neighbours_ranks[DOWN]);
 
    // shift tells us our left and right neighbours
    MPI_Cart_shift(cartesian_communicator, 1, 1, &neighbours_ranks[LEFT], &neighbours_ranks[RIGHT]);

 
    // From rookiehpc.com	
    // Partition the 2D cartesian topology along the first dimension, by preserving the second dimension
    int remain_dims[2] = {false, true};
    MPI_Comm subgrid_communicator;
    MPI_Cart_sub(cartesian_communicator, remain_dims, &subgrid_communicator);
  
    //rank in the subgrid communicator
    int sub_rank;
    MPI_Comm_rank(subgrid_communicator, &sub_rank);
    
	    
 
 
    //Variables for the timer
    double start_time, stop_time, elapsed_time;
    elapsed_time=0;
    
    if (rank==0){
      printf("============================================\n");
      printf("Number of MPI Processes: %d\n", size);
      printf("Number of entries: %d\n",Matrix_Size*Matrix_Size);
      printf("============================================\n");
    }
    
    // Iterate for verification
    for (int u=1;u<=ITER_OUTER_LOOP;u++){
    
    start_time = MPI_Wtime(); 
    
    // ===================== Implement the Algorithm ========================
 
    //Create local matrices to divide
    int Lsize=Matrix_Size/sqrt(size);
    // How many rows and columns in each tile
    int r = Lsize;
    int c = Lsize;
    // Allocate space for the tiles and comm buffers 
    double **A, **Atemp, **B, **Btemp, **C;
    A= alloc_2d_array(Lsize, Lsize);
    Atemp= alloc_2d_array(Lsize, Lsize);
    B= alloc_2d_array(Lsize, Lsize);
    Btemp= alloc_2d_array(Lsize, Lsize);
    C= alloc_2d_array(Lsize, Lsize);
    
    // Perform the domain decomposition - Copy entries from the global to local. 
    for (i = 0; i <  r; i++)
      for (j = 0; j < c; j++)
        A[i][j] = matrix_a[i+cart_coords[0]*Lsize][j+cart_coords[1]*Lsize]; //
        
    for (i = 0; i <  r; i++)
      for (j = 0; j < c; j++)
        B[i][j] = matrix_b[i+cart_coords[0]*Lsize][j+cart_coords[1]*Lsize]; // 
        
    for (i = 0; i <  r; i++)
      for (j = 0; j < c; j++)
        C[i][j] = matrix_c[i+cart_coords[0]*Lsize][j+cart_coords[1]*Lsize]; // 
        
  	
  	//Follow the Fox algorithm
       int totaliterations=sqrt(size);
       for (int it=0;it<totaliterations;it++){
      
       // First copy A before broadcasting. So you dont overwrite it. 
       for (i = 0; i <  r; i++)
         for (j = 0; j < c; j++)
           Atemp[i][j] = A[i][j]; //
        
       // broadcast the diagonal+it      	
	int broadcast_root=cart_coords[0]+it;
  	//If the entry is over number of horizontal tiles (squareroot(size))-1 since the rank start at 0. Then start from 0.
  	if (broadcast_root>(int)sqrt(size)-1){
  	broadcast_root=cart_coords[0]+it-(int)sqrt(size);
  	}
  	//Perform the broadcast
  	MPI_Bcast(&Atemp[0][0], Lsize*Lsize, MPI_DOUBLE, broadcast_root, subgrid_communicator);

	//Do the matmul for the tile
  	for (i = 0 ; i < Lsize ; i++) {
    	  for (j = 0 ; j < Lsize ; j++) {
      	    for (k = 0 ; k < Lsize ; k++) {
              C[i][j] += Atemp[i][k] * B[k][j];
      	    }
    	  }
  	}
  	
  	//send B up recieve down
    	MPI_Sendrecv(&B[0][0], Lsize*Lsize, MPI_DOUBLE, neighbours_ranks[UP], 0, &Btemp[0][0], Lsize*Lsize, MPI_DOUBLE, neighbours_ranks[DOWN], MPI_ANY_TAG, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
   	
  	//Update the recieved value in preparation from the next iteration
	for (i = 0; i <  r; i++)
         for (j = 0; j < c; j++)
           B[i][j] = Btemp[i][j]; 
  	}
  	
  	
  	
  // ===================== Print and validate results ========================

  
 //Calculate a global average of each tile
  double  ave = 0.0;
  double  maindiag=0.0;
  for (i = 0 ; i < Lsize ; i++) {
    for (j = 0 ; j < Lsize ; j++) {
      ave += C[i][j]/(double)(Matrix_Size*Matrix_Size);
      if (i==j && cart_coords[0]==cart_coords[1]){
        maindiag += C[i][j]/(double)(Matrix_Size*Matrix_Size);
      }
    }
  }
  
  //Communicate the results and perform reduction
  double reduction_result = 0;
  double reduction_diag = 0;
  int root_rank=rankprint;
  MPI_Reduce(&ave, &reduction_result, 1, MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
  MPI_Reduce(&maindiag, &reduction_diag, 1, MPI_DOUBLE, MPI_SUM, root_rank, MPI_COMM_WORLD);
  
  
  int row, columns;
  //Pint the matrix
  if (rank==0){
  	
  	printf("Average of elements of matrix for run %d = %8.6f \n", u, reduction_result);
  	printf("Trace of matrix: %d = %8.6f \n", u, reduction_diag);
 
	stop_time = MPI_Wtime();
	elapsed_time = elapsed_time + stop_time - start_time;
    	
    	printf("Execution time for run %d: %f seconds\n", u, stop_time - start_time);
        printf("------------------------------------------------------------------\n");
  }
  
  }
  
  if (rank==rankprint){
  	printf("============================================\n");
    	printf("Average execution time: %f seconds\n", elapsed_time/ITER_OUTER_LOOP);
        printf("============================================\n");
  }
  
 
      MPI_Finalize();

  return 0;
}



