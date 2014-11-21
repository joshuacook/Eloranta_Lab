/*
  This file has my implementation of the LAPACK routine zgeev for C++.
  This program solves for the eigenvalues and, if desired, the
  eigenvectors for a complex Hamiltonian matrix H.  Both real and
  imaginary components of the eigenvalues are returned.

  There are two function calls defined in this header, of the
  forms

    void zgeev(complex<double> **H, int n, complex<double> *E)
    void zgeev(complex<double> **H, int n, complex<double> *E,
               complex<double> **Evecs)

    H: the n by n Hamiltonian matrix that we are solving
    n: the order of the square matrix H
    E: an n-element array to hold the real parts of the eigenvalues
       of H.
    Evecs: an n by n matrix to hold the eigenvectors of H, if
           they are requested.

  The function call is defined twice, so that whether or not
  eigenvectors are called for the proper version is called.

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
  28 January 1999 */


#include <math.h>
#include <complex.h>

void zgeev(complex<double> **H, int n, complex<double> *E);
void zgeev(complex<double> **H, int n, complex<double> *E,
	   complex<double> **Evecs);

complex<double> *zgeev_ctof(complex<double> **in, int rows, int cols);
void zgeev_ftoc(complex<double> *in, complex<double> **out,
		int rows, int cols);

void zgeev_sort(int n, complex<double> *E);
void zgeev_sort(int n, complex<double> *E, complex<double> **Evecs);


extern "C" void zgeev_(char *jobvl, char *jobvr, int *n, complex<double> *a,
		       int *lda, complex<double> *w, complex<double> *vl,
		       int *ldvl, complex<double> *vr, int *ldvr,
		       complex<double> *work, int *lwork, double *rwork,
		       int *info);

void zgeev(complex<double> **H, int n, complex<double> *E)
{
  char jobvl, jobvr;
  int lda, ldvl, ldvr, lwork, info;
  double *rwork;
  complex<double> *a, *w, *vl, *vr, *work;
  
  jobvl = 'N'; /* V/N to calculate/not calculate the left eigenvectors
		  of the matrix H.*/
  
  jobvr = 'N'; // As above, but for the right eigenvectors.
  
  lda = n; // The leading dimension of the matrix a.
  a = zgeev_ctof(H, n, lda); /* Convert the matrix H from double pointer
				C form to single pointer Fortran form. */

  /* Whether we want them or not, we need to define the matrices
     for the eigenvectors, and give their leading dimensions.
     We also create a vector for work space. */

  ldvl = n;
  vl = new complex<double>[ldvl*n];
  ldvr = n;
  vr = new complex<double>[ldvr*n];
  lwork = 4*n;
  work = new complex<double>[lwork];
  rwork = new double[2*n];

  zgeev_(&jobvl, &jobvr, &n, a, &lda, E, vl, &ldvl, vr,
	 &ldvr, work, &lwork, rwork, &info);

  zgeev_sort(n, E);

  // Clean up the memory usage

  delete a;
  delete vl;
  delete vr;
  delete work;
  delete rwork;
}


void zgeev(complex<double> **H, int n, complex<double> *E,
	   complex<double> **Evecs)
{
  char jobvl, jobvr;
  int lda, ldvl, ldvr, lwork, info;
  double *rwork;
  complex<double> *a, *w, *vl, *vr, *work;
  
  jobvl = 'N'; /* V/N to calculate/not calculate the left eigenvectors
		  of the matrix H.*/
  
  jobvr = 'V'; // As above, but for the right eigenvectors.
  
  lda = n; // The leading dimension of the matrix a.
  a = zgeev_ctof(H, n, lda); /* Convert the matrix H from double pointer
				C form to single pointer Fortran form. */

  /* Whether we want them or not, we need to define the matrices
     for the eigenvectors, and give their leading dimensions.
     We also create a vector for work space. */

  ldvl = n;
  vl = new complex<double>[ldvl*n];
  ldvr = n;
  vr = new complex<double>[ldvr*n];
  lwork = 4*n;
  work = new complex<double>[lwork];
  rwork = new double[2*n];

  zgeev_(&jobvl, &jobvr, &n, a, &lda, E, vl, &ldvl, vr,
	 &ldvr, work, &lwork, rwork, &info);

  // Now I need to get the eigenvectors into the output array

  zgeev_ftoc(vr, Evecs, ldvr, n);
  zgeev_sort(n, E, Evecs);

  // Clean up the memory usage

  delete a;
  delete vl;
  delete vr;
  delete work;
  delete rwork;
}


complex<double>* zgeev_ctof(complex<double> **in, int rows, int cols)
{
  complex<double> *out;
  int i, j;

  out = new complex<double>[rows*cols];
  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i+j*cols] = in[i][j];
  return(out);
}


void zgeev_ftoc(complex<double> *in, complex<double> **out,
		 int rows, int cols)
{
  int i, j;

  for (i=0; i<rows; i++) for (j=0; j<cols; j++) out[i][j] = in[i+j*cols];
}


// Sort the eigenenergies by their real components.

void zgeev_sort(int n, complex<double> *E)
{
  complex<double> temp;
  double min, v;
  int imin, ct;

  for (int j=0; j<n-1; j++) {
    min = fabs(E[j].real());
    imin = j;
    for (int i=j+1; i<n; i++) {
      if (fabs(E[i].real())<min) {
	min = fabs(E[i].real());
	imin = i;
      }
    }
    if (imin!=j) {
      temp = E[imin];
      E[imin] = E[j];
      E[j] = temp;
    }
  }
}


void zgeev_sort(int n, complex<double> *E, complex<double> **Evecs)
{
  complex<double> temp;
  double min, v;
  int imin, ct;

  for (int j=0; j<n-1; j++) {
    min = fabs(E[j].real());
    imin = j;
    for (int i=j+1; i<n; i++) {
      if (fabs(E[i].real())<min) {
	min = fabs(E[i].real());
	imin = i;
      }
    }
    if (imin!=j) {
      temp = E[imin];
      E[imin] = E[j];
      E[j] = temp;

      for (int i=0; i<n; i++) {
	temp = Evecs[i][imin];
	Evecs[i][imin] = Evecs[i][j];
	Evecs[i][j] = temp;
      }
    }
  }
}
