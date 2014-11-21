/*
  This is an example program using the header file zgeev.h, my
  implementation of the LAPACK complex matrix solver routine.  See
  that header file for more documentation on its use.
  
  For this example, we will work with the four by four matrix

    1+i 2    3     4
    5   6+2i 7     8
    9   10   11+3i 12
    14  14   15    16+4i

  on the left hand side.  To solve an Ax=b problem, we take the vector

    1 2 3 4

  on the right hand side.  We also fully invert the same matrix.  The
  output of this program should be:

    (-0.0351124,-0.0688202)
    (0.0660112,-0.00561798)
    (0.0997191,0.0154494)
    (0.133778,0.0172051)

    (-0.275281,-0.359551)	(-0.0575843,0.167135)	(0.0266854,0.0323034)	(0.0688202,-0.0351124)	
    (-0.0224719,0.235955)	(0.0182584,-0.379213)	(0.00983146,0.0842697)	(0.00561798,0.0660112)	
    (0.0617978,0.101124)	(0.0435393,0.105337)	(0.00421348,-0.231742)	(-0.0154494,0.0997191)	
    (0.19382,-0.0351124)	(0.0143961,0.0832163)	(-0.00667135,0.116924)	(-0.0172051,-0.116222)	

  Scot Shaw
  28 January 2000 */

#include <iostream.h>
#include "zgesv.h"

const complex<double> I(0,1);

void main()
{
  complex<double> **A, *b, *x;
  int n=4;

  // First we solve for a single right hand side

  // Allocate memory

  A = new complex<double>*[n];
  for (int i=0; i<n; i++) A[i] = new complex<double>[n];
  b = new complex<double>[n];
  x = new complex<double>[n];

  // Generate the matrices

  for (int i=0; i<n; i++) for (int j=0; j<n; j++) A[i][j] = i*n+j+1;
  for (int i=0; i<n; i++) A[i][i] += (i+1)*I;
  A[3][0] += 1;

  for (int i=0; i<n; i++) b[i] = i+1;

  // Call the solver

  zgesv(A, b, n, x);

  // Print the results

  for (int i=0; i<n; i++) cout << x[i] << "\n";

  // Now we completely invert the same matrix

  // Allocate memory

  complex<double> **AI;
  AI = new complex<double>*[n];
  for (int i=0; i<n; i++) AI[i] = new complex<double>[n];

  // Generate the matrix

  for (int i=0; i<n; i++) for (int j=0; j<n; j++) A[i][j] = i*n+j+1;
  for (int i=0; i<n; i++) A[i][i] += (i+1)*I;
  A[3][0] += 1;

  // Call the solver

  zgesv(A, AI, n);

  // Print the results

  cout << endl;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) cout << AI[i][j] << "\t";
    cout << endl;
  }
}
