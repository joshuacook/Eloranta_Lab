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

void print_matrix(double * matrix, int n);

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

  for (int col = 0; col < n; col++){
    *(U + col*n) = *(A + col*n);
  }

  for (int i = 1; i < n; i++){

  }

  print_matrix(U);

  return 0;

}

void print_matrix(double * matrix, int n){
  printf("\n");
  for (int col = 0; col < n; col++){
    for (int row = 0; row < n; row++){
      printf("\t%f",*(matrix+row*3+col)); // *(matrix+row*3+col) used to travers column major matrix
    }
    printf("\n");
  }
  printf("\n");
}
