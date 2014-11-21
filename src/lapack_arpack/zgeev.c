/*
  This is an example program using the header file zgeev.h, my
  implementation of the LAPACK complex matrix solver routine.  See
  that header file for more documentation on its use.
  
  For this example, we will work with the four by four matrix

    1+i 2    3     4
    5   6+2i 7     8
    9   10   11+3i 12
    14  14   15    16+4i

  The output of this program should be:

    Eigenvalue  0: (36.279,3.14823)
    Eigenvector 0: (0.150123,-0.0114246) (0.347632,-0.0184352) (0.545571,-0.0145593) (0.747185,0) 
    Eigenvalue  1: (-2.15437,1.81526)
    Eigenvector 1: (0.673556,0) (0.247071,-0.251264) (-0.191337,-0.100045) (-0.511194,0.337946) 
    Eigenvalue  2: (-0.0550235,3.19742)
    Eigenvector 2: (-0.12827,0.0978439) (-0.0509522,0.10703) (0.77725,0) (-0.57738,-0.149788) 
    Eigenvalue  3: (-0.0695699,1.83909)
    Eigenvector 3: (-0.435554,0.203597) (0.816338,0) (-0.101738,0.0171173) (-0.258111,-0.158645) 

  Scot Shaw
  28 January 2000 */

#include <iostream.h>
#include "zgeev.h"

const complex<double> I(0,1);

void main()
{
  complex<double> **H, *E, **Evecs;
  int n=4;

  // Allocate memory

  H = new complex<double>*[n];
  for (int i=0; i<n; i++) H[i] = new complex<double>[n];
  Evecs = new complex<double>*[n];
  for (int i=0; i<n; i++) Evecs[i] = new complex<double>[n];

  // Generate the matrix

  for (int i=0; i<n; i++) for (int j=0; j<n; j++) H[i][j] = i*n+j+1;
  for (int i=0; i<n; i++) H[i][i] += (i+1)*I;
  H[3][0] += 1;

  E = new complex<double>[n];

  // Call the solver

  zgeev(H, n, E, Evecs);

  // Print the results

  for (int i=0; i<n; i++) {
    cout << "Eigenvalue  " << i << ": " << E[i] << "\n";
    cout << "Eigenvector " << i << ": ";
    for (int j=0; j<n; j++) cout << Evecs[j][i] << " ";
    cout << "\n";
  }
}




