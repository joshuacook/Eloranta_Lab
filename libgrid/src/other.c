#include "grid.h"
#include "private.h"

/*
 * Polynomial interpolation routine.
 *
 * xa    = Array of x values (double *).
 * ya    = Array of y values (double *).
 * n     = Number of points for (x,y) (long).
 * x     = Point at which the approximation is obtained (double).
 * dy    = Error estimate (double *).
 *
 * Returns the interpolated value.
 *
 * Source: Numerical Recipes (indexing from zero).
 *
 */

EXPORT double grid_polynomial_interpolate(double *xa, double *ya, long n, double x, double *dy) {

  long i, m, ns = 1;
  double den, dif, dift, ho, hp, w, *c, *d, y;

  /* Indexing: preserve Num. Recp's messed up indexing and just shift xa and ya accordingly */
  dif = fabs(x - xa[0]);
  if(!(c = (double *) malloc(sizeof(double) * (n+1))) || !(d = (double *) malloc(sizeof(double) * (n+1)))) {
    fprintf(stderr, "libgrid: Out of memory in grid_polynomial_interpolate().\n");
    exit(1);
  }
  
  for (i = 1; i <= n; i++) {
    if((dift = fabs(x - xa[i-1])) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i-1];
    d[i] = ya[i-1];
  }
  y = ya[ns-1]; ns--;
  for (m = 1; m < n; m++) {
    for (i = 1; i <= n - m; i++) {
      ho = xa[i-1] - x;
      hp = xa[i + m - 1] - x;
      w = c[i+1] - d[i];
      if((den = ho - hp) == 0.0) {
	fprintf(stderr, "libgrid: Polynomial interpolation failed.\n");
	exit(1);
      }
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    y += (*dy = (2 * ns < (n - m) ? c[ns+1]:d[ns--]));
  }
  free(c);
  free(d);
  return y;
}

/* Second derivative for spline interpolation 
 * x    = Array of x values (double *).
 * y    = Array of y values (double *).
 * n    = Number of points for (x,y) (long).
 * yp1	= First derivative at x[1] (double)
 * ypn 	= First derivative at x[n+1] (double)
 * y2	= Second derivative computed. 
 *
 * Must be called only once before calling grid_spline_interpolate
 *
 * If yp1 and yp2 are larger than 1.e30 the second derivative at the boundaries is 0.
 *
 * Source: Numerical Recipes (indexing from zero).
 *
 */
EXPORT  void grid_spline_ypp( double *x, double *y, long n, double yp1, double ypn , double *y2 ){
	long i, k;
	double p, qn, sig, un, *u ;

  if(!(u = (double *) malloc(sizeof(double) * n )) )  {
    fprintf(stderr, "libgrid: Out of memory in grid_spline_ypp().\n");
    exit(1);
  }
  if(yp1 > 0.99e30)
	  y2[0]=u[0]=0.0 ;
  else{
	  y2[0]= -0.5 ;
	   u[0]= (3.0 / (x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1) ;
  }
  for(i=1 ; i <= n-2 ; i++) { 
	  sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]) ;
	  p = sig*y2[i-1] + 2.0 ;
	  y2[i] = (sig-1.0) / p ;
	  u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]) ;
	  u[i] = ( 6.0*u[i] / (x[i+1]-x[i-1] ) - sig*u[i-1]) / p ;
  }
  if(ypn > 0.99e30)
	  qn = un = 0.0 ;
  else{
	  qn = 0.5 ;
	  un = (3.0 / (x[n]-x[n-1])) * (ypn - (y[n]-y[n-1]) / (x[n]-x[n-1])) ;
  }
  y2[n-1] = ( un - qn * u[n-2] ) / ( qn*y2[n-2]+1.0 ) ;
  for(k=n-2 ; k>=0 ; k--)
	  y2[k] = y2[k]*y2[k+1] + u[k];
 free(u) ;
}


/* Same as before but allocating and returning the y2 */
EXPORT double *grid_spline_ypp_new( double *x, double *y, long n, double yp1, double ypn){
	double *y2 ;
  if(!(y2 = (double *) malloc(sizeof(double) * n )) )  {
    fprintf(stderr, "libgrid: Out of memory in grid_spline_ypp_new().\n");
    exit(1);
  }
  grid_spline_ypp(x, y, n, yp1, ypn, y2) ;
  return y2 ;
}

/* Spline interpolation 
 * xa    = Array of x values (double *).
 * ya    = Array of y values (double *).
 * y2a	 = Array of second derivatives of y values (double *). 
 *         Computed in grid_spline_ypp
 * n     = Number of points for (x,y) (long).
 * x     = Point at which the approximation is obtained (double).
 *
 * Returns the interpolated value.
 *
 * Source: Numerical Recipes (indexing from zero).
 *
 */
EXPORT double grid_spline_interpolate(double *xa, double *ya, double *y2a, long n, double x){
	long klo, khi, k;
	double h, b, a ;

	klo = 0;
	khi = n-1;
	while (khi-klo > 1) {
		k = (khi+klo) >> 1 ;
	       	if(xa[k] > x) khi = k ;
		else klo = k ;
	}
	h = xa[khi]-xa[klo] ;
	if(h==0.0){
    		fprintf(stderr, "libgrid: Bad x input values for grid_spline_interpolation() (repeated value).\n");
    		exit(1);
  	}
	a = ( xa[khi] - x ) / h ;
	b = ( x - xa[klo] ) / h ;
	return a * ya[klo] + b * ya[khi] + ( (a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi] ) * (h*h) / 6.0 ;
}

