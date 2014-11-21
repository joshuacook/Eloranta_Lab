/*
  This file has my implementation of the LAPACK routine dgesv for C++.
  This program solves the set of complex linear equations of the form
  Ax = b, where A is a given n by n tridiagonal complex matrix and b
  is a given n element vector.  The n element vector x is returned.
  The function call is of the form

    void zgtsv(complex<double> *Al, complex<double> *Am,
               complex<double> *Au, complex<double> *b, int n,
               complex<double> *x)

    Al: the n-1 sub-diagonal elements of A
    Am: the n diagonal elements of A
    Au: the n-1 super-diagonal elements of A
    b: the right hand side n element vector
    n: the dimension of A and b
    x: the n element vector for returned values

  There is a second form of the call, wherein we assume that the right
  hand side is the identity matrix; hence we are finding the inverse
  of the matrix A.
  
    void zgtsv(complex<double> *Al, complex<double> *Am,
               complex<double> *Au, int n, complex<double> **X)

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
  7 March 2000 */


#include <math.h>
#include <complex.h>

void zgtsv(complex<double> *Al, complex<double> *Am,
           complex<double> *Au, complex<double> *b, int n,
	   complex<double> *x);

void zgtsv_ftoc(complex<double> *in, complex<double> **out,
		int rows, int cols);

extern "C" void zgtsv_(int *n, int *nrhs, complex<double> *al, 
		       complex<double> *am, complex<double> *au,
		       complex<double> *b, int *ldb, int *info);

void zgtsv(complex<double> *Al, complex<double> *Am,
           complex<double> *Au, complex<double> *b, int n,
	   complex<double> *x)
{
  int nrhs, ldb, info;
  int *ipiv;

  nrhs = 1; /* The Fortran routine that I am calling will solve a
	       system of equations with multiple right hand sides.
	       nrhs gives the number of right hand sides, and then b
	       must be a matrix of dimension n by nrhs. */
  ldb = n;

  /* The Fortran routine replaces the input information in b with the
     results.  Since we are interested in the results, we put the
     initial conditions in the output array and use that in the
     call. */

  for (int i=0; i<n; i++) x[i] = b[i];

  /* The Fortran routine also screws up the elements of the
     tridiagonal matrix that is fed in.  On the chance that I will
     need this information again in the calling program, I'll make a
     duplicate of the information to send to Fortran. */

  complex<double> *cAl, *cAm, *cAu;

  cAl = new complex<double>[n-1];
  cAm = new complex<double>[n];
  cAu = new complex<double>[n-1];

  for (int i=0; i<n; i++) cAm[i] = Am[i];
  for (int i=0; i<n-1; i++) { cAl[i] = Al[i]; cAu[i] = Au[i]; }

  // Now the function call

  zgtsv_(&n, &nrhs, cAl, cAm, cAu, x, &ldb, &info);

  delete(cAl);
  delete(cAm);
  delete(cAu);
}


void zgtsv(complex<double> *Al, complex<double> *Am,
           complex<double> *Au, int n, complex<double> **X)
{
  int nrhs, ldb, info;
  int *ipiv;
  complex<double> *x;

  nrhs = n;
  ldb = n;
  x = new complex<double>[n*n];

  for (int i=0; i<n; i++) x[i*(n+1)] = 1;
  
  complex<double> *cAl, *cAm, *cAu;

  cAl = new complex<double>[n-1];
  cAm = new complex<double>[n];
  cAu = new complex<double>[n-1];

  for (int i=0; i<n; i++) cAm[i] = Am[i];
  for (int i=0; i<n-1; i++) { cAl[i] = Al[i]; cAu[i] = Au[i]; }

  zgtsv_(&n, &nrhs, cAl, cAm, cAu, x, &ldb, &info);

  // Convert to the C form of the matrix
  
  zgtsv_ftoc(x, X, ldb, nrhs);

  delete(cAl);
  delete(cAm);
  delete(cAu);
}


void zgtsv_ftoc(complex<double> *in, complex<double> **out,
                 int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*cols];
}



