#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/usr/local/Cellar/openblas/0.2.12/include/cblas.h"
#include "/Users/joshuacook/src/Eloranta_Lab/my_blas/blas_fortran_double.h"

// FORTRAN Routines used
// dgemv_dgemv_(&no_trans, &two, &two, &alpha, B, &two, x, &one, &beta, y, &one);
// y := alpha*A*x + beta*y

// dscal_(&two, &mu_inv, y, &one);
// x = a*x

// dcopy_(&two, y,&one, x,&one);


int main(int argc, char* argv[])
{
  int i;
  printf("A simple power method iteration to find largest eigenvector\n");

  double B[2*2] = 
    {  2,  1,
     -12,  -5};                 // matrix to be examined
     
  double alpha = 1.0;           // scaling factor
  double beta = 1.0;            // scaling factor 
  double mu;                    // norm
  double mu_inv;                // inverse norm

  int one = 1;                  // fortran requires pointers
  int two = 2;
  char no_trans='N';                    

  // Initialize and display a random vector 
  double x[2] ; 
  x[0] = rand()%10;
  x[1] = rand()%10;
  printf("x: %f, %f\n",x[0], x[1]);
  
  // Will hold resulting vector
  double y[2] ;

  for (i = 0; i < 4; i++){
    // Perform matrix multiplication
    // row_order  transform lenY lenX alpha  a  lda  X  incX  beta  Y, incY 
    // y := alpha*A*x + beta*y
    printf("Multiply A*x to get y\n");
    dgemv_(&no_trans, &two, &two, &alpha, B, &two, x, &one, &beta, y, &one);
    
    // Normalize Vector
    // elements X incX Y incY 
    // x = a*x
    printf("Calculate mu, the magnitude of the vector, y\n");
    mu = sqrt(ddot_ (&two, y, &one, y, &one));    
    mu_inv = 1/mu;
    printf("mu: %f\n",mu);

    // elements alpha X intX Y intY(y:= a*x+y)
    printf("y: %f, %f\n",y[0], y[1]);
    printf("Normalize the vector y\n");
    dscal_(&two, &mu_inv, y, &one);
    
    // Display result
    printf("Display new y\n");
    printf("y: %f, %f\n",y[0], y[1]);   
    
    // Copy into x for next iteration  

    dcopy_(&two, y,&one, x,&one);
  }


  printf("%f, %f\n", y[0],y[1] );

  return 0;
}
