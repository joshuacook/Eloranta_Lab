/* This header contains multiple utility functions used in working within the BLAS framework */ 
/* Stylistically, BLAS functions typically take pointers as input. These functions are taking the variables themselves. */

#include <time.h>

#define NEWLINE printf("\n");

// Function Declarations
void print_matrix(double * A, int n, int m);
void print_vector(double * A, int n);
void random_matrix(double * A, int n, int m);
void random_vector(double *A, int n);
void identity_matrix(double * I, int n);

// Matrix and Vector Generators for testing
void random_matrix(double * A, int m, int n){
  time_t timer;
  int i;

  srand48(time(&timer));

  for (i = 0; i < m*n; i++){
    A[i] = drand48();
  }
}
void random_vector(double * A, int n){
  time_t timer;
  int i; 
  
  srand48(time(&timer));
  
  for (i = 0; i < n; i++){
    A[i] = drand48();
  } 
} 
void identity_matrix(double * I, int n){
  int i;
  
  for(i = 0; i < n + 1; i++){
    I[i*(n+1)] = 1;
  } 
} 

// Printing Functions
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
void print_matrix(double * A, int n, int m){
  int i,j;
  for (i = 0; i < m; i++){
    for (j = 0; j < n; j++){
      printf("%f ", (float)A[j*m+i]);
    } 
    NEWLINE;
  } 
  NEWLINE;
} 

