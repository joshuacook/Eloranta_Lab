#include <cblas.h>
#include <stdio.h>

int main()
{
  int i;
  double a[2*2] = {  1,  0,
                     0,  1};
  double x[2]  = {2,3};
  double y[2];
  
        	/* row_order      transform     lenY lenX alpha  a  lda  X  incX  beta  Y, incY */
   cblas_dgemv(102, 111, 2,   2,   1,     a,   2, x, 1,    0,    y, 1);

  for(i=0; i<2; i++)
		printf("%lf ", y[i]);
  printf("\n");

  return 0;
}


// gcc -o test_cblas_open test_cblas_open.c -I ~/bin/OpenBLAS/ -L ~/bin/OpenBLAS/ -lopenblas -lpthread -lgfortran