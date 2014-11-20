#include <stdio.h>
#include <stdlib.h>
#include "../../OpenBLAS/cblas.h"
#include <math.h>

extern void dgemv_(char* Trans, int* M, int* N, double* alpha, double* A, int* lda, double* X, int* incX, double* beta, double* Y, int* incY );
extern double ddot_(int* N,double* X,int* incX,double* Y,int* incY);
extern void dscal_(int* N,double* alpha,double* X, int* incX); 
extern void dcopy_(int* N,double* X, int* incX, double* Y, int* incY);
int main(int argc, char* argv[])
{
	int i;
	printf("test!\n");

	double B[2*2] = {  2,  1,
					 -12,  -5};
  double alpha = 1.0;
  double beta = 1.0;					 
	double mu;
    double mu_inv;
    int one = 1;
    int two = 2;
    char no_trans='N';					 
	double x[2] ;	
	x[0] = rand()%10;
	x[1] = rand()%10;
	printf("x: %f, %f\n",x[0], x[1]);
	

	double y[2] ;

	for (i = 0; i < 100; i++){
	      	// row_order  transform lenY lenX alpha  a  lda  X  incX  beta  Y, incY 
		dgemv_(&no_trans, &two, &two, &alpha, B, &two, x, &one, &beta, y, &one);
			// elements X incX Y incY 
		mu = sqrt(ddot_ (&two, y, &one, y, &one));
		mu_inv = 1/mu;
        printf("mu: %f\n",mu);
			// elements alpha X intX Y intY(y:= a*x+y)
		printf("y: %f, %f\n",y[0], y[1]);
		dscal_(&two, &mu_inv, y, &one);
		printf("y: %f, %f\n",y[0], y[1]);		
		dcopy_(&two, y,&one, x,&one);
		printf("y: %f, %f\n",y[0], y[1]);
	}


	printf("%f, %f\n", y[0],y[1] );

	return 0;
}
