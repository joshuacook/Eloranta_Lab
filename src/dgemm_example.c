#define min(x,y) (((x) < (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

int main()
{
  double *A, *B, *C;
  int m,n,k,i,j;
  double alpha, beta;

  printf ("\n This example computes real matrix C=alpha*A*B+beta*C using \n"
          " dgemm, where A, B, and C are matrices and alpha and beta are \n"
          " double precision scalars\n\n");

  m = 2000; k = 200; n = 1000;
  printf (" Initializing data for matrix multiplication C=A*B for matrix \n"
          " A(%ix%i) and matrix B(%ix%i)\n\n", m, k, k, n);
  alpha = 1.0; beta = 0.0;

  printf (" Allocating memory for matrices aligned on 64-byte boundary for \n"
          " performance \n\n");

  A = (double *)mkl_malloc: 
}
