// toward a solution to Ax=b
// factor A = L*U
// then LUx=b
// forward substitute to solve Ly = b
// back substitute to solve Ux = y

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/usr/local/Cellar/openblas/0.2.12/include/cblas.h"
#include "blas_fortran_double.h"
#include "blas_utilities.h"


int main(int argc, char* argv[]){
  printf("Factor a Matrix into its upper triangular portion\n");

  int n = 3; 

  double A[n*n];                // initial matrix
 
  double U[n*n];                // to hold factored matrix

  print_matrix(A, n, n);

  for (int col = 0; col < n; col++){
    *(U + col*n) = *(A + col*n);
  }

  for (int i = 1; i < n; i++){

  }

  print_matrix(U, n, n);

  return 0;

}




