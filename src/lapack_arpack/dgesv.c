/*
  This is an example program using the header file dgesv.h, my
  implementation of the LAPACK linear equation solver.  See that
  header file for more documentation on its use.
  
  The matrix equation that this example solves has

    A = 2 5 9 13 
        2 6 10 14 
	2 7 13 21 
	2 8 18 40 

    b = 2 4 8 16

  The output should be

    -5.33333 1.66667 8.88178e-16 0.333333

  Scot Shaw
  1 February 2000 */

// Begin with some standard include files.

#include <math.h>
#include <stdio.h>
#include <fstream.h>
#include "dgesv.h"

void main()
{
  double **A, *b, *x;
  int n=4;
  
  // Generate the matrices for my equation

  A = new double*[n];
  for (int i=0; i<n; i++) A[i] = new double[n];
  for (int i=0; i<n; i++) for (int j=0; j<n; j++)
    A[i][j] = pow(i,j)+n*j+1;
  
  b = new double[n];
  for (int i=0; i<n; i++) b[i] = pow(2,i+1);

  // Call the LAPACK solver

  x = new double[n];
  
  dgesv(A, b, n, x);
  
  // Output the results

  for (int i=0; i<n; i++) cout << x[i] << " "; cout << "\n";
}






