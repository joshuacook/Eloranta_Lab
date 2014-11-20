#include <stdio.h>

extern void dgemv_(char* trans,int*  m,int* n,double* alpha,double*  a,int* lda,double*  x,int* incx, double* beta,double* y,int* incy);

int main()
{
  int i;
  int m = 2;
  int n = 2;
  int incX = 1;
  int incY = 1;
  double alpha = 1;
  double beta = 0;
  double a[2*2] = {  1,  0,
                     0,  1};
  double x[2]  = {2,3};
  double y[2];
  char trans = 'N';
  dgemv_(&trans,&m,&n,&alpha,a,&m,x,&incX,&beta,y,&incY);

  for(i=0; i<2; i++)
		printf("%lf ", y[i]);
  printf("\n");

  return 0;
}


// gcc -o test_cblas_open test_cblas_open.c -I ~/bin/OpenBLAS/ -L ~/bin/OpenBLAS/ -lopenblas -lpthread -lgfortran
