
#include <stdio.h>
#include <stdlib.h>
#include "blas_fortran_double.h"

#define NEWLINE printf("\n");

void print_vector(double * A, int n);


int main()
{
  int n=4;
  double mu;
  int one = 1;

  double A[4] = {1,1,1,1};

  print_vector(A,n);

  mu = dnrm2_(&n,A,&one);

  printf("%f\n", mu);

}

void print_vector(double * A, int n){
  int i;
  printf("<");
  for (i = 0; i < n-1; i++){
    printf("%f, ", (float)A[i]);
  }
  printf("%f", (float)A[i]);
  printf(">");
  NEWLINE;
  NEWLINE;
}
