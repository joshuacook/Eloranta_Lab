#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/usr/local/Cellar/openblas/0.2.12/include/cblas.h"
#include "blas_fortran_double.h"
#include "blas_utilities.h"

// FORTRAN Routines used
// dgemv_dgemv_(&no_trans, &n, &n, &alpha, B, &n, x, &one, &beta, y, &one);
// y := alpha*A*x + beta*y

// dscal_(&n, &mu_inv, y, &one);
// x = a*x

// dcopy_(&n, y,&one, x,&one);


int main(int argc, char* argv[])
{
  float eps = 10e-6;
  int i;
  printf("A simple power method iteration to find largest eigenvector\n");

  int n = 4;

  double A[n*n];                 // matrix to be examined
  double A_2[n*n];

  random_matrix(A, n, n, 1234);
  print_matrix(A, n, n);


  double alpha = 1.0;           // scaling factor
  double beta = 1.0;            // scaling factor 
  double mu;                    // norm
  double mu_inv;                // inverse norm

  int one = 1;                  // fortran requires pointers
  int two = 2;
  char no_trans='N';                    
  float conv = 1.0;

  // Calculate A^2
  dgemm_(&no_trans, &no_trans, &n, &n, &n, &alpha, A, &n, A, &n, &alpha, A_2, &n);

  // Initialize and display a random vector 
  double x[n] ; 
  random_vector(x,n);
  print_vector(x,n);
  
  // Will hold resulting vector
  double y[n] ;
  double temp[n];
  float var_1, var_2;

  while (conv > eps){
    // Perform matrix multiplication
    // row_order  transform lenY lenX alpha  a  lda  X  incX  beta  Y, incY 
    // y := alpha*A*x + beta*y
    printf("\tMultiply A*x to get y\n");
    dgemv_(&no_trans, &n, &n, &alpha, A, &n, x, &one, &beta, y, &one);
    
    // Normalize Vector
    // elements X incX Y incY 
    // x = a*x
    printf("\tCalculate mu, the magnitude of the vector, y\n");
    mu = sqrt(ddot_ (&n, y, &one, y, &one));    
    mu_inv = 1/mu;
    printf("\tmu: %f\n",mu);

    // elements alpha X intX Y intY(y:= a*x+y)
    print_vector(y,n);
    printf("\tNormalize the vector y\n");
    dscal_(&n, &mu_inv, y, &one);
    
    // Display result
    printf("\tDisplay new y\n");
    print_vector(y,n);
    
    // Copy into x for next iteration  

    dcopy_(&n, y,&one, x,&one);
    
    // stopping criteria
    dgemv_(&no_trans, &n, &n, &alpha, A_2, &n, y, &one, &beta, temp, &one);
    var_1 = ddot_ (&n, temp, &one, y, &one);
    printf("var_1: %f\n", var_1 );
    dgemv_(&no_trans, &n, &n, &alpha, A, &n, y, &one, &beta, temp, &one);
    var_2 = ddot_ (&n, temp, &one, y, &one);
    var_2 = var_2*var_2;
    printf("var_2: %f\n", var_2 );
    conv = sqrt(2)*sqrt(abs(var_1-var_2));
    printf("%f\n", conv );
  }

  return 0;
}
