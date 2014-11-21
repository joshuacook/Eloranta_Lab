/*
  This is an example program using the header file
  znaupd.h, my implementation of the ARPACK sparse matrix
  solver routine.  See that header file for more
  documentation on its use.
  
  For this example, we will work with the four by four matrix H[i,j] =
  i^j + j^i + I delta(i) delta(j-3).  The output of this program
  should be:

  Eigenvalue  0: (0.150577,-0.0321152)
  Eigenvector 0: ((-0.253177,-0.557293), (0.326799,0.641205),
                 (-0.152555,-0.285208), (0.025262,0.0461506))
  Eigenvalue  1: (1.80622,0.0783858)
  Eigenvector 1: ((0.573869,-0.502308), (0.287917,-0.305961),
                 (-0.344501,0.335538), (0.0742376,-0.0712374))
  
  Scot Shaw
  7 September 1999 */

// Begin with some standard include files.

#include <math.h>
#include <stdio.h>
#include <fstream.h>
#include <complex.h>

/*
  To use the routines from znaupd.h, we need to define our own matrix
  multiplication routine.  Here I show one method that would exploit
  sparsity by storing only the non-zero matrix elements (though for
  this example we won't actually have a sparse matrix to work with).

  There is a routine for the multiplication, and a global variable T
  that will hold the non-zero matrix elements and their indicies.  The
  global integer tlen will hold the number of non-zero elements.  The
  matrix multiplication function needs to have the form below, and
  must be defined before the file znaupd.h is included.  */

void av(int n, complex<double> *in, complex<double> *out);
#include "znaupd.h"

const complex<double> I(0,1);

typedef struct {
  int i, j;
  complex<double> v;
} Telt;

Telt *T;
int tlen;

void main()
{
  int n, nev, i, j;
  complex<double> **Evecs, *Evals;
  
  n = 4; // The order of the matrix

  /*
    Here we generate the matrix T for the multiplication.  It helps if
    we know the number of non-zero matrix elements before we find
    them, so that we can save on storage space.  If not, generate
    enough storage space for T to cover an upper limit.  */

  tlen = n*n;
  T = new Telt[tlen];

  tlen = 0;
  for (i=0; i<n; i++) for (j=0; j<n; j++) {
    T[tlen].i = i;
    T[tlen].j = j;
    T[tlen].v = pow(i+1, j+1) + pow(j+1, i+1);
    if ((i==0)&&(j==3)) T[tlen].v += I;
    tlen++;
  }

  /*
    We will calculate both the eigenvalues and the eigenvectors in
    this example.  To not calculate the vectors, just leave that
    argument out of the call to dsaupd.

    We will find only three of the eigenvalues, to make it more clear
    how the eigenvectors come out of the routine.  Besides that, the
    ARPACK routine won't calculate all of the eigenvalues.  I'm not
    sure if that is intentional, or if I am missing something.  */

  nev = 2; // The number of values to calculate

  Evals = new complex<double>[nev];
  Evecs = new complex<double>*[n];
  for (i=0; i<n; i++) Evecs[i] = new complex<double>[nev];

  znaupd(n, nev, Evals, Evecs);

  // Now we output the results.
  
  for (i=0; i<nev; i++) {
    cout << "Eigenvalue  " << i << ": " << Evals[i] << "\n";
      cout << "Eigenvector " << i << ": ("
      << Evecs[0][i] << ", "
      << Evecs[1][i] << ", "
      << Evecs[2][i] << ", "
      << Evecs[3][i] << ")\n";
  }
}


void av(int n, complex<double> *in, complex<double> *out)
{
  int i, j;

  for (i=0; i<n; i++) out[i] = 0;
  for (i=0; i<tlen; i++) out[T[i].i] += in[T[i].j] * T[i].v;
}
