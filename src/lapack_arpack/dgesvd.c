/*
  This is an example program using the header file dgesvd.h, my
  implementation of the LAPACK singular value decomposition routine.
  See that header file for more information about the routine.

  For this example, we work with the four by three matrix A[i,j] =
  n*i+j+1.  The output of the program should be

    1	2	3	
    4	5	6	
    7	8	9	
    10	11	12	
    
    25.4624	1.29066	2.21303e-15	
    
    -0.140877	0.824714	0.540521	
    -0.343946	0.426264	-0.654711	
    -0.547016	0.0278135	-0.312142	
    -0.750086	-0.370637	0.426331	
    
    -0.504533	-0.574516	-0.644498	
    -0.760776	-0.0571405	0.646495	
    -0.408248	0.816497	-0.408248	

  First we have the matrix being decomposed.  Next are the three
  singular values.  Then the matrix U, and finally the matrix VT.

  Scot Shaw
  24 January 2000 */

// Begin with some standard include files.

#include <math.h>
#include <stdio.h>
#include <fstream.h>
#include "dgesvd.h"

void main(int nargs, char *args[])
{
  double **A, *S, **U, **VT;
  int i, j, k;

  int m=4; // rows
  int n=3; // columns
  int minmn = 3;

  // Generate the matrix for my equation

  A = new double*[m];
  for (int i=0; i<m; i++) A[i] = new double[n];
  for (int i=0; i<m; i++) for (int j=0; j<n; j++) A[i][j] = n*i+j+1;

  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) cout << A[i][j] << "\t";
    cout << "\n";
  } cout << "\n";

  // Allocate the storage space for the output arrays and call the solver

  S = new double[minmn];
  U = new double*[m]; for (int i=0; i<m; i++) U[i] = new double[minmn];
  VT = new double*[minmn]; for (int i=0; i<minmn; i++) VT[i] = new double[n];

  dgesvd(A, m, n, S, U, VT);

  // Output the results

  for (int i=0; i<minmn; i++) cout << S[i] << "\t"; cout << "\n\n";

  for (int i=0; i<m; i++) {
    for (int j=0; j<minmn; j++) cout << U[i][j] << "\t"; cout << "\n";
  } cout << "\n";

  for (int i=0; i<minmn; i++) {
    for (int j=0; j<n; j++) cout << VT[i][j] << "\t"; cout << "\n";
  } cout << "\n";
}









