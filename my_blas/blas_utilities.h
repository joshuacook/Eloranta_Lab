/* This header contains multiple utility functions used in working within the BLAS framework */ 
/* Stylistically, BLAS functions typically take pointers as input. These functions are taking the variables themselves. */

#include <time.h>

#define NEWLINE printf("\n");

// Function Declarations
void random_matrix(double * A, int n, int m);
// void print_matrix(double * A, int n, int m);
// void print_vector(double * A, int n);
// void random_vector(double *A, int n);
// void identity_matrix(double * I, int n);

// Matrix and Vector Generators for testing
void random_matrix(double * A, int m, int n){
  time_t timer;
  int i;

  srand48(time(&timer));

  for (i = 0; i < m*n; i++){
    *(A + i) = drand48();
  }
}

