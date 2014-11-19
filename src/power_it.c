#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>
#include <math.h>

int main(int argc, char* argv[])
{
	int i;
	printf("test!\n");

	double B[2*2] = {  2,  1,
					 -12,  -5};
					 
	double mu;
					 
	double x[2] ;	
	x[0] = rand()%10;
	x[1] = rand()%10;
	printf("x: %f, %f\n",x[0], x[1]);
	

	double y[2] ;

	for (i = 0; i < 100; i++){
	      	// row_order  transform lenY lenX alpha  a  lda  X  incX  beta  Y, incY 
		cblas_dgemv(102,111, 2, 2, 1, B, 2, x, 1, 1, y, 1);
			// elements X incX Y incY 
		mu = sqrt(cblas_ddot (2, y, 1, y, 1));
		printf("mu: %f\n",mu);
			// elements alpha X intX Y intY(y:= a*x+y)
		printf("y: %f, %f\n",y[0], y[1]);
		cblas_dscal (2, 1/mu, y, 1);
		printf("y: %f, %f\n",y[0], y[1]);		
		cblas_dcopy (2, y, 1, x, 1);
		printf("y: %f, %f\n",y[0], y[1]);
	}


	printf("%f, %f\n", y[0],y[1] );

	return 0;
}