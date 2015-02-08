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

void print_matrix(double * matrix);

int main(int argc, char* argv[]){
  printf("Factor a Matrix into its upper triangular portion\n");

  int n = 3; 

  double A[n*n] = { 2, 1, 3, -3, 4, 7, -12, -5, 3};                
  // matrix to be examined
  // note that this matrix is in column major order and represents the following matrix
  //     2  -3  -12
  // A=  1   4   -5
  //     3   7    3

  double U[n*n] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};                // to hold factored matrix

  print_matrix(A);

  for (int i = 0; i < n; i++){
    *(U + i*n) = *(A + i*n);
  }

  print_matrix(U);

  return 0;

}

void print_matrix(double * matrix, n){
  printf("\n");
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      printf("\t%f",*(matrix+j*n+i));
    }
    printf("\n");
  }
  printf("\n");
}
