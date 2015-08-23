#include <stdio.h>
#include <stdlib.h>
// #include "../../OpenBLAS/cblas.h"
#include <math.h>

int main(int argc, char* argv[])
{
	int i;
	printf("test!\n");

	double B[2*2] = {  2,  1,
					 -12,  -5};
    int one = 1;
    int two = 2;
    char no_trans='N';					 
	double mu;
					 
	double x[2] ;	
	x[0] = rand()%10;
	x[1] = rand()%10;
	printf("x: %f, %f\n",x[0], x[1]);
	
    double y[2] = {0,0};
    
    dgemv_(&no_trans, &two, &two, &one, &B, &two, &x, &one, &one, &y, &one);


	printf("%f, %f\n", y[0],y[1] );

	return 0;
}
