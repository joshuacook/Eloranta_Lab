/*
  This is an example program using the header file zgtsv.h, my
  implementation of the LAPACK complex matrix solver routine.  See
  that header file for more documentation on its use.
  
  For this example, we will work with the four by four matrix

    2+I  -3    0   0
    -1    2  -2+I  0
     0  -2+I   2   -1
     0    0   -3   2

  on the left hand side, and the vector

    1 2 3 4+I

  on the right hand side.  We find first the solution to this set of
  linear equations, and then find the inverse of the matrix.  The
  output of this program should be:

    (-2.82574,-0.942574)
    (-1.90297,-1.5703)
    (-0.752475,-1.47525)
    (0.871287,-1.71287)

    (0.340594,-0.205941) (-0.112871,-0.0712871) (-0.594059,-0.0594059)
      (-0.29703,-0.029703)	
    (-0.0376238,-0.0237624) (-0.0514851,-0.0851485) (-0.376238,-0.237624)
      (-0.188119,-0.118812)	
    (-0.19802,-0.019802) (-0.376238,-0.237624) (0.019802,-0.19802)
      (0.00990099,-0.0990099)	
    (-0.29703,-0.029703) (-0.564356,-0.356436) (0.029703,-0.29703)
      (0.514851,-0.148515)	

  Scot Shaw
  7 March 2000 */

#include <iostream.h>
#include "zgtsv.h"

const complex<double> I(0,1);

void main()
{
  complex<double> *Al, *Am, *Au, *b, *x;
  int n=4;

  // Allocate memory

  Al = new complex<double>[n-1];
  Am = new complex<double>[n];
  Au = new complex<double>[n-1];
  b = new complex<double>[n];
  x = new complex<double>[n];

  // Generate the matrices

  for (int i=0; i<n; i++) Am[i] = 2;
  for (int i=0; i<n-1; i++) {
    Al[i] = -(i+1);
    Au[i] = -(3-i);
  }
  Al[1] += I;
  Am[0] += I;
  Au[1] += I;

  for (int i=0; i<n; i++) b[i] = i+1;
  b[3] += I;

  // Call the solver

  zgtsv(Al, Am, Au, b, n, x);

  // Print the results

  for (int i=0; i<n; i++) cout << x[i] << "\n";


  // Now try an inversion of the same matrix

  // Generate the matrices

  for (int i=0; i<n; i++) Am[i] = 2;
  for (int i=0; i<n-1; i++) {
    Al[i] = -(i+1);
    Au[i] = -(3-i);
  }
  Al[1] += I;
  Am[0] += I;
  Au[1] += I;

  // Allocate space for the return matrix

  complex<double> **Ai;
  Ai = new complex<double>*[n];
  for (int i=0; i<n; i++) Ai[i] = new complex<double>[n];

  // Call the solver

  zgtsv(Al, Am, Au, n, Ai);
  
  // Print the results

  cout << endl;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) cout << Ai[i][j] << "\t";
    cout << endl;
  }
}

