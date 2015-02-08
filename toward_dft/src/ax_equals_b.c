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

  double A[3*3] = { 2, 1, 3, -3, 4, 7, -12, -5, 3};                // matrix to be examined
  // note that this matrix is in column major order and represents the following matrix
  //     2  -3  -12
  // A=  1   4   -5
  //     3   7    3
  
  double U[3*3] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};                // too hold factored matrix

  print_matrix(A);

  return 0;

}

void print_matrix(double * matrix){
  for (int i = 0; i < 3; i++){
    for (int j = 0; j < 3; j++){
      printf("\t%f",*(matrix+j*3+i));
    }
    printf("\n");
  }
}
