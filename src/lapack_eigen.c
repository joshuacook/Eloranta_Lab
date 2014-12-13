/* finding the eigenvalues of a complex matrix */

#include <stdio.h>

// #include 
#define size 3				/* dimension of matrix */

extern void dgeev_(char *jobvl, char *jobvr, int *n, double * a,
	int *lda, double * w, double * vl,
	int *ldvl, double * vr, int *ldvr,
	double * work, int *lwork, double *rwork, int *info);

int main()
{
double A[3][3], b[3], DUMMY, WORK[6];
double AT[2*size*size];			/* for transformed matrix */
int i, j, ok, c1, c2, c3;
char c4;

A[0][0]=3.1;			/* the input matrix */
A[0][1]=1.3;	
A[0][2]=-5.7;	
A[1][0]=1.0;	
A[1][1]=-6.9;	
A[1][2]=5.8;	
A[2][0]=3.4;	
A[2][1]=7.2;	
A[2][2]=-8.8;	

for (i=0; i<size; i++)		/* to call a Fortran routine from C we */
{				/* have to transform the matrix */
  for(j=0; j<2*size; j++) 
  {
     AT[2*(j+size*i)]=A[j][i];
  }		
}

c1=size;			/* and put all numbers and characters */ 
c2=2*size;    			/* we want to pass */
c3=1;				/* to the routine in variables */
c4='N';

/* find solution using LAPACK routine ZGEEV, all the arguments have to */
/* be pointers and you have to add an underscore to the routine name */
dgeev_(&c4, &c4,&c1, AT, &c1, b, &DUMMY, &c3, &DUMMY, &c3, WORK, &c2, WORK, &ok);      

/*
 parameters in the order as they appear in the function call
    no left eigenvectors, no right eigenvectors, order of input matrix A,
    input matrix A, leading dimension of A, array for eigenvalues, 
    array for left eigenvalue, leading dimension of DUMMY, 
    array for right eigenvalues, leading dimension of DUMMY,
    workspace array dim>=2*order of A, dimension of WORK
    workspace array dim=2*order of A, return value */

if (ok==0)				/* output of eigenvalues */
{
   for (i=0; i<size; i++)
   {
      printf("%f\t%f\n", b[i], b[i]);
   }
}
else printf("An error occured");
}
