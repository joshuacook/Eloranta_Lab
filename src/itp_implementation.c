#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "blas_fortran_double.h"
#include "lapack_fortran_double.h"
#include "blas_utilities.h"
#include "itp_method.h"

int itp_method_test(double * H, int n, int print_mat, int DEBUG);


int main(int argc, char* argv[]){

  // Performance Testing
  double time_elapsed_in_seconds;

  // Matrix parameters
  int n = 100;
  double A[n*n];

  // create matrix to be used for all
  random_matrix(A,n,n);

  //--------------------- ITP Method Test  
  // start timer
  clock_t start = clock();
  // call method
  // get first eigenvector
  itp_method_test(A,2,1,0);  

  // end timer
  clock_t end = clock();  
  // store time
  printf("Found eigenvector in %.2g seconds\n", 
    time_elapsed_in_seconds = (end - start)/(double)CLOCKS_PER_SEC);
  
  //--------------------- ARPACK Method
  // start timer
  // call method
  // get first eigenvector
  // end timer
  // store time
  
  //--------------------- LAPACK Method
  // start timer
  // call method
  // get first eigenvector
  // end timer
  // store time
  
  return 0;
}
