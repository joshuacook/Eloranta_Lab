/*
  This is an example program using the header file
  dgeev.h, my implementation of the LAPACK matrix
  solver routine.  See that header file for more
  documentation on its use.
  
  For this example, we will work with the four by four
  matrix H[i,j] = i^j + j^i.  This is symmetric, so the
  eigenvalues are real.  The output of this program
  should be:
  
  Eigenvalue  0: 0.151995 + i 0
  Eigenvector 0: (0.610828, -0.72048, 0.324065, -0.0527271)
  Eigenvalue  1: 1.80498 + i 0
  Eigenvector 1: (-0.7614, -0.42123, 0.481884, -0.103072)
  Eigenvalue  2: 17.6352 + i 0
  Eigenvector 2: (-0.216885, -0.547084, -0.764881, 0.26195)
  Eigenvalue  3: 556.408 + i 0
  Eigenvector 3: (0.0110019, 0.064609, 0.278795, 0.958112)

  Scot Shaw
  7 September 1999
*/

// Begin with some standard include files.

#include <math.h>
#include <stdio.h>
#include <fstream>
#include "dgeev.h"

void main(int nargs, char *args[])
{
  double **H, *Er, *Ei, **Evecs;
  int i, j, k;
  
  int n=4;
  
  // Generate the matrix for my equation
  
  H = new double*[n];
  for (i=0; i<n; i++) H[i] = new double[n];
  for (i=0; i<n; i++) for (j=0; j<n; j++)
    H[i][j] = pow(i+1, j+1) + pow(j+1, i+1);
  
  // Call the LAPACK solver
  
  Er = new double[n]; 
  Ei = new double[n]; 
  Evecs = new double*[n];
  for (i=0; i<n; i++) Evecs[i] = new double[n];

  dgeev(H, n, Er, Ei, Evecs);
  
  // Output the results

  for (i=0; i<n; i++) {
    cout << "Eigenvalue  " << i << ": " << Er[i] << " + i " << Ei[i] << "\n";
      cout << "Eigenvector " << i << ": ("
      << Evecs[0][i] << ", "
      << Evecs[1][i] << ", "
      << Evecs[2][i] << ", "
      << Evecs[3][i] << ")\n";
  }
}




