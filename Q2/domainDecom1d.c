#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(int argc, char *argv[]){

    int rank, size, i, provided;
    
    // number of cells (global)
    int nxc = 128; // make sure nxc is divisible by size
    double L = 2*3.1415; // Length of the domain
    int templ,tempr;
    

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of nodes (local to the process): 0 and nxn_loc-1 are ghost cells
    // node 1 and nxn_loc-2 are shared by contigous cells
    // so the ones that need to be sent over are 2 and nxn_loc-3
    int nxn_loc = nxc/size + 3; // number of nodes is number cells + 1; we add also 2 ghost cells
    double L_loc = L/((double) size);
    double dx = L / ((double) nxc);
    
    // define out function
    double *f = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double *dfdx = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double *exact_df = calloc(nxn_loc, sizeof(double)); // allocate and fill with z

    for (i=1; i<(nxn_loc-1); i++)
      f[i] = sin(L_loc*rank + (i-1) * dx);
    
    // need to communicate and fill ghost cells f[0] and f[nxn_loc-1]
    // communicate ghost cells
     
    
     
    // Keep into consideration periodicity for rank 0 and size-1
    if (rank==0)
    {
         MPI_Sendrecv(&f[2],1,MPI_DOUBLE,size-1,0, &f[nxn_loc-1],1, MPI_DOUBLE,rank+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
         MPI_Sendrecv(&f[nxn_loc-3],
                     1,
                     MPI_DOUBLE,
                     rank+1,
                     0,
                     &f[0],
                     1,
                     MPI_DOUBLE,
                     size-1,
                     MPI_ANY_TAG,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
     }  
     else if (rank==size-1)
     {
      MPI_Sendrecv(&f[2],
                     1,
                     MPI_DOUBLE,
                     rank-1,
                     0,
                     &f[nxn_loc-1],
                     1,
                     MPI_DOUBLE,
                     0,
                     MPI_ANY_TAG,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
      MPI_Sendrecv(&f[nxn_loc-3],
                     1,
                     MPI_DOUBLE,
                     0,
                     0,
                     &f[0],
                     1,
                     MPI_DOUBLE,
                     rank-1,
                     MPI_ANY_TAG,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

     }   
     else
     {
      MPI_Sendrecv(&f[2],
                     1,
                     MPI_DOUBLE,
                     rank-1,
                     0,
                     &f[nxn_loc-1],
                     1,
                     MPI_DOUBLE,
                     rank+1,
                     MPI_ANY_TAG,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
       MPI_Sendrecv(&f[nxn_loc-3],
                     1,
                     MPI_DOUBLE,
                     rank+1,
                     0,
                     &f[0],
                     1,
                     MPI_DOUBLE,
                     rank-1,
                     MPI_ANY_TAG,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

     }   

    

    // here we finish the calculations

    // calculate first order derivative using central difference
    // here we need to correct value of the ghost cells!
    for (i=1; i<(nxn_loc-1); i++)
      dfdx[i] = (f[i+1] - f[i-1])/(2*dx);

    for (i=1; i<(nxn_loc-1); i++){
    // Real value of derivative
      exact_df[i] = cos(L_loc*rank + (i-1) * dx);
    }
    
    
      // Print f values
      if (rank==0){ // print only rank 0 for convenience
    
          printf("My rank %d of %d\n", rank, size );
          printf("Here are my values for dfdx-cos() including ghost cells\n");
          for (i=0; i<nxn_loc; i++)
	       printf("sin(x)=%f, numerical_dfdx=%f, cos(x)=%f\n", f[i], dfdx[i], exact_df[i]);
          printf("\n");
    }

    MPI_Finalize();
}





