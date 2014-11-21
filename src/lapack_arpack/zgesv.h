/*
  This file has my implementation of the LAPACK routine dgesv for C++.
  This program solves the general set of complex linear equations of
  the form Ax = b, where A is a given n by n complex matrix and b is a
  given n element vector.  The n element vector x is returned.  The
  function call is of the form

    void zgesv(complex<double> **A, complex<double> *b, int n,
               complex<double> *x)

    A: the left hand size n by n matrix
    b: the right hand side n element vector
    n: the dimension of A and b
    x: the n element vector for returned values

  There is also a specialized version of the call, to be used when we
  wish to fully invert the matrix A.  In this case, one uses

    void zgesv(complex<double> **A, complex<double> **AI, int n)

    A: the matrix to invert
    AI: the matrix for returning the inverse
    n: the dimension of A and AI.

  This function uses the C++ object type complex<double> found in the
  header file complex.h to interface with the Fortran type complex16.
  Briefly, here is how to work with these complex numbers.  Declare
  them as usual.  To add to the real component of the number, use +=
  with a real right hand side.  To add to the imaginary part of the
  number, I have had to declare a number complex<double> I(0,1) which
  makes I be the imaginary number i.  Then one can use a command like
  "x += .3*I" to add to the imaginary component of a complex variable
  x.  (There may be a better way of doing this, but I haven't found it
  yet.)  For more information about using these complex numbers, take
  a look at the header file complex.h.

  Scot Shaw
  1 December 2000 */


#include <math.h>
#include <complex.h>

void zgesv(complex<double> **A, complex<double> *b, int n,
	   complex<double> *x);
void zgesv(complex<double> **A, complex<double> **AI, int n);

complex<double> *zgesv_ctof(complex<double> **in, int rows, int cols);
void zgesv_ftoc(complex<double> *in, complex<double> **out,
                int rows, int cols);

extern "C" void zgesv_(int *n, int *nrhs, complex<double> *a, int *lda,
		       int *ipiv, complex<double> *b, int *ldb, int *info);

void zgesv(complex<double> **A, complex<double> *b, int n, complex<double> *x)
{
  int nrhs, lda, ldb, info;
  complex<double> *a;
  int *ipiv;

  nrhs = 1; /* The Fortran routine that I am calling will solve a
	       system of equations with multiple right hand sides.
	       nrhs gives the number of right hand sides, and then b
	       must be a matrix of dimension n by nrhs. */

  lda = n; // The leading dimension of A
  a = zgesv_ctof(A, lda, n); // Convert A to a Fortran style matrix

  ipiv = new int[n]; // Allocate memory to store pivot information

  ldb = n;

  /* The Fortran routine replaces the input information in b with the
     results.  Since we are interested in the results, we put the
     initial conditions in the output array and use that in the
     call. */

  for (int i=0; i<n; i++) x[i] = b[i];

  // Now the function call

  zgesv_(&n, &nrhs, a, &lda, ipiv, x, &ldb, &info);

  // Clean up the memory before returning to the calling program

  delete ipiv;
  delete a;
}


void zgesv(complex<double> **A, complex<double> **AI, int n)
{
  int nrhs, lda, ldb, info;
  complex<double> *a, *b;
  int *ipiv;

  nrhs = n; // The number of right hand sides
  lda = n; // The leading dimension of A
  a = zgesv_ctof(A, lda, n); // Convert A to a Fortran style matrix

  ipiv = new int[n]; // Allocate memory to store pivot information

  ldb = n;

  /* To invert the matrix we solve an AX=B problem with X = AI and B
     the identity matrix.  We thus need to generate the identity
     matrix as an initial right hand side. */

  b = new complex<double>[n*n];
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
    b[i*n+j] = 0.0;
    if (i==j) b[i*n+j] = 1.0;
  }

  // Now the function call

  zgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

  // Convert the returned inverse to the C-style matrix AI.

  zgesv_ftoc(b, AI, n, n);

  // Clean up the memory before returning to the calling program

  delete ipiv;
  delete b;
  delete a;
}


complex<double>* zgesv_ctof(complex<double> **in, int rows, int cols)
{
  complex<double> *out;
  int i, j;

  out = new complex<double>[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*cols] = in[i][j];
  return(out);
}


void zgesv_ftoc(complex<double> *in, complex<double> **out,
                 int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*cols];
}


