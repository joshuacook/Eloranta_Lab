/*
 * Routines for 2D complex cylindrical grids.
 *
 */

#include "grid.h"
#include "private.h"
#include "private2d.h"

/*
 * Multiply a given grid by a function (cylindrical coordinates).
 *
 * grid = destination grid for the operation (cgrid2d *).
 * func = function providing the mapping (double complex (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and z,r are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_product_func_cyl(cgrid2d *grid, double complex (*func)(void *arg, double z, double r), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double z, r, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,step,func,value) private(i,j,ij,z,r) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    z = (i - nx/2.0) * step;
    r = j * step;
    value[ij] *= func(farg, z, r);
  }
}

/*
 * Map a given function onto 2D cylindrical grid.
 *
 * grid = destination grid for the operation (cgrid2d *).
 * func = function providing the mapping (double complex (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and z,r are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_map_cyl(cgrid2d *grid, double complex (*func)(void *arg, double z, double r), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double z, r, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,step,func,value) private(i,j,ij,z,r) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    z = (i - nx/2.0) * step;
    r = j * step;
    value[ij] = func(farg, z, r);
  }
}

/*
 * Map a given function onto 2D cylindrical grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid.
 * *
 * grid = destination grid for the operation (cgrid2d *).
 * func = function providing the mapping (double complex (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and (z,r) is the point (doubles) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 * ns   = number of intermediate points to be used in smoothing (int).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_smooth_map_cyl(cgrid2d *grid, double complex (*func)(void *arg, double z, double r), void *farg, int ns) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double zc, rc, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,ns,step,func,value) private(i,j,zc,rc) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;    
    zc = (i - nx/2.0) * step;
    rc = j * step;
    value[ij] = linearly_weighted_integralc2d(func, farg, zc, rc, step, ns);
  }
}

/*
 * Map a given function onto 2D cylindrical grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid. Limits for intermediate steps and
 * tolerance can be given.
 *
 * grid   = destination grid for the operation (cgrid2d *).
 * func   = function providing the mapping (double complex (*)(void *, double, double, double)).
 *          The first argument (void *) is for external user specified data
 *          and (z,r) is the point (doubles) where the function is evaluated.
 * farg   = pointer to user specified data (void *).
 * min_ns = minimum number of intermediate points to be used in smoothing (int).
 * max_ns = maximum number of intermediate points to be used in smoothing (int).
 * tol    = tolerance for weighing (double).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_adaptive_map_cyl(cgrid2d *grid, double complex (*func)(void *arg, double z, double r), void *farg, int min_ns, int max_ns, double tol) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny, ns;
  double zc, rc, step = grid->step;
  double tol2 = tol * tol;
  double sum, sump;
  double complex *value = grid->value;
  
  if (min_ns < 1) min_ns = 1;
  if (max_ns < min_ns) max_ns = min_ns;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,min_ns,max_ns,step,func,value,tol2) private(i,j,ns,zc,rc,sum,sump) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;    
    zc = (i - nx/2.0) * step;
    rc = j * step;          
    sum  = func(farg, zc, rc); sump = 0.0;
    for(ns = min_ns; ns <= max_ns; ns *= 2) {
      sum  = linearly_weighted_integralc2d(func, farg, zc, rc, step, ns);
      sump = linearly_weighted_integralc2d(func, farg, zc, rc, step, ns+1);
      if (sqnorm(sum - sump) < tol2) break;
    }
    
    /*
      if (ns >= max_ns)
      fprintf(stderr, "#");
      else if (ns > min_ns + 1)
      fprintf(stderr, "+");
      else
      fprintf(stderr, "-");
    */
    
    value[ij] = 0.5 * (sum + sump);
  }
  
  /*fprintf(stderr, "\n");*/
}

/*
 * Integrate over a grid (cylindrical coordinates).
 *
 * grid = grid to be integrated (cgrid2d *).
 *
 * Returns the integral value (double complex).
 *
 */

EXPORT double complex cgrid2d_integral_cyl(const cgrid2d *grid) {

  long ij, j, ny = grid->ny, nxy = grid->nx * ny;
  double complex sum = 0.0;
  double complex *value = grid->value;
  double r, step = grid->step;
  
#pragma omp parallel for firstprivate(nxy,value,step,ny) private(ij,r,j) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    j = ij % ny;
    r = step * (double) j;
    sum += value[ij] * r;
  }
  return sum * grid->step * grid->step * 2.0 * M_PI;
}

/* 
 * Integrate over the grid squared (int |grid|^2) (cylindrical coordinates).
 *
 * grid = grid to be integrated (cgrid2d *).
 *
 * Returns the integral (double complex).
 *
 */

EXPORT double cgrid2d_integral_of_square_cyl(const cgrid2d *grid) {

  long j, ij, ny = grid->ny, nxy = grid->nx * ny;
  double sum = 0.0, r, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,ny,step) private(ij,j,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    j = ij % ny;
    r = step * (double) j;
    sum += sqnorm(value[ij]) * r;
  }
  
  return sum * grid->step * grid->step * 2.0 * M_PI;
}

/*
 * Calculate expectation value of cylindrical grid over the probability density given by the other.
 * (int gridb grida^2).
 *
 * grida = grid giving the probability (grida^2) (cgrid2d *).
 * gridb = grid to be averaged (cgrid2d *).
 *
 * Returns the average value (double complex).
 *
 */

EXPORT double complex cgrid2d_grid_expectation_value_cyl(const cgrid2d *grida, const cgrid2d *gridb) {

  long ij, j, ny = grida->ny, nxy = grida->nx * ny;
  double step = grida->step;
  double complex *avalue = grida->value, r;
  double complex *bvalue = gridb->value, sum = 0.0;
  
#pragma omp parallel for firstprivate(ny,nxy,avalue,bvalue,step) private(ij,j,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    j = ij % ny;
    r = step * (double) j;
    sum += avalue[ij] * avalue[ij] * bvalue[ij] * r;
  }
  
  return sum * step * step * 2.0 * M_PI;
}

/*
 * Calculate the expectation value of a function over a grid.
 * (int grida func grida = int func |grida|^2).
 *
 * func  = function to be averaged (double complex (*)(void *, double, double, double)).
 *         The arguments are: optional arg, grida(z,r), z, r.
 * grida = grid giving the probability (grida^2) (cgrid2d *).
 *
 * Returns the average value (double complex).
 *
 */
 
EXPORT double complex cgrid2d_grid_expectation_value_func_cyl(void *arg, double complex (*func)(void *arg, double complex val, double z, double r), const cgrid2d *grida) {
   
  long ij, i, j, nx = grida->nx, ny = grida->ny, nxy = nx * ny;
  double complex sum = 0.0;
  double complex *avalue = grida->value;
  double z, r, step = grida->step;
  
#pragma omp parallel for firstprivate(func,arg,nxy,nx,ny,avalue,step) private(z,r,i,j,ij) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    z = (i - nx/2.0) * step;
    r = j * step;
    sum += sqnorm(avalue[ij]) * func(arg, avalue[ij], z, r) * r;
  }
  
  return sum * step * step * 2.0 * M_PI;
}

/* 
 * Integrate over cylindrical grid multiplied by weighting function (int grid w(z,r)).
 *
 * grid   = grid to be integrated over (cgrid2d *).
 * weight = function defining the weight (double complex (*)(double, double)). The arguments are (z,r) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double complex cgrid2d_weighted_integral_cyl(const cgrid2d *grid, double complex (*weight)(void *farg, double z, double r), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, ny = grid->ny, nx = grid->nx;
  double z, r, step = grid->step;
  double complex sum = 0.0;
  double complex w;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nx,ny,step,value,weight,farg) private(w,i,j,ij,z,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;    
    z = (i - nx/2.0) * step;
    r = step * (double) j;
    w = weight(farg, z, r);
    sum += w * value[ij] * r;
  }
  
  return sum * step * step * 2.0 * M_PI;
}

/* 
 * Integrate over square of the grid multiplied by weighting function (int grid^2 w(z,r)).
 *
 * grid   = grid to be integrated over (cgrid2d *).
 * weight = function defining the weight (double complex (*)(double, double)).
 *          The arguments are (z,r) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double cgrid2d_weighted_integral_of_square_cyl(const cgrid2d *grid, double complex (*weight)(void *farg, double z, double r), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double z, r, step = grid->step;
  double complex sum = 0.0, w;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nx,ny,step,value,weight,farg) private(w,i,j,ij,z,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {    
    i = ij / ny;
    j = ij % ny;    
    z = (i - nx/2.0) * step;
    r = j * step;
    w = weight(farg, z, r);
    sum += w * sqnorm(value[ij]) * r;
  }
  
  return sum * step * step * 2.0 * M_PI;
}

/* 
 * Differentiate a grid with respect to z.
 * This is used in dot product: grad . f along z.  (just a derivative with respect to z)
 *
 * grid     = grid to be differentiated (cgrid2d *).
 * gradient = differentiated grid output (cgrid2d *).
 * 
 * No return value.
 *
 */

EXPORT void cgrid2d_fd_gradient_z_cyl(const cgrid2d *grid, cgrid2d *gradient) {

  long i, j, ij;
  long ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double complex *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in cgrid2d_fd_gradient_cyl_z().\n");
    return;
  }
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta,grid) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = inv_delta * (cgrid2d_value_at_index(grid, i+1, j) - cgrid2d_value_at_index(grid, i-1, j));
  }
}

/* 
 * Differentiate a grid with respect to r, i.e. (1/r) * d(rf(r))/dr = df(r)/dr + (1/r)f(r).
 * This used in dot product: grad . f along r. Note that gradient in cylindrical coordinates
 * does NOT have the additional 1/r term.
 *
 * grid     = grid to be differentiated (cgrid2d *).
 * gradient = differentiated grid output (cgrid2d *).
 * 
 * No return value.
 *
 */

EXPORT void cgrid2d_fd_gradient_r_cyl(const cgrid2d *grid, cgrid2d *gradient) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double step = grid->step;
  double inv_delta = 1.0 / (2.0 * step);
  double complex *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in rgrid2d_fd_gradient_cyl_r().\n");
    return;
  }
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta,grid,step) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = inv_delta * (cgrid2d_value_at_index(grid, i, j+1) - cgrid2d_value_at_index(grid, i, j-1));
  }
}

/*
 * Calculate gradient of a cylindrical grid.
 *
 * grid       = grid to be differentiated twice (cgrid2d *).
 * gradient_z = z output grid for the operation (cgrid2d *).
 * gradient_r = r output grid for the operation (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_fd_gradient_cyl(const cgrid2d *grid, cgrid2d *gradient_z, cgrid2d *gradient_r) {

  cgrid2d_fd_gradient_z_cyl(grid, gradient_z);
  cgrid2d_fd_gradient_r_cyl(grid, gradient_r);
}

/*
 * Calculate laplacian of cylindrical grid.
 *
 * grid    = source grid (cgrid2d *).
 * laplace = output grid for the operation (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_fd_laplace_cyl(const cgrid2d *grid, cgrid2d *laplace) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * ny;
  double step = grid->step, inv_delta2 = 1.0 / (step * step);
  double inv_delta = 1.0 / (2.0 * step), r;
  double complex *lvalue = laplace->value;
  
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta2,inv_delta,grid,step) private(ij,i,j,r) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    r = step * (double) j;
    lvalue[ij] = inv_delta2 * (-4.0 * cgrid2d_value_at_index(grid, i, j) + cgrid2d_value_at_index(grid, i, j+1) + cgrid2d_value_at_index(grid, i, j-1) + cgrid2d_value_at_index(grid, i+1, j) + cgrid2d_value_at_index(grid, i-1, j));
    if(r != 0.0)
      lvalue[ij] += inv_delta * (cgrid2d_value_at_index(grid, i, j+1) - cgrid2d_value_at_index(grid, i, j-1)) / r;
  }
}

/*
 * Access grid point at given (z,r) point (with linear interpolation).
 *
 * grid = grid to be accessed (cgrid2d *).
 * x    = x value (double).
 * y    = y value (double).
 *
 * Returns grid value at (x,y).
 *
 */

EXPORT inline double complex cgrid2d_value_cyl(const cgrid2d *grid, double z, double r) {

  double complex f00, f10, f01, f11;
  double omr, omz, step = grid->step;
  long i, j;
  
  /* i to index and 0 <= z < 1 */
  i = (z /= step);
  if (z < 0) i--;
  z -= i;
  i += grid->nx / 2;
  
  /* j to index and 0 <= r < 1 */
  r = r / step;
  j = (long) r;
  r = r - j;
  
  /* linear extrapolation 
   *
   * f(z,r) = (1-z) (1-r) f(0,0) + x (1-r)f(1,0) + (1-z) r f(0,1)
   *          + z     r   f(1,1)
   */ 
  f00 = cgrid2d_value_at_index(grid, i, j);
  f10 = cgrid2d_value_at_index(grid, i+1, j);
  f01 = cgrid2d_value_at_index(grid, i, j+1);
  f11 = cgrid2d_value_at_index(grid, i+1, j+1);
  
  omz = 1.0 - z;
  omr = 1.0 - r;

  return omz * omr * f00 + z * omr * f10 + omz * r * f01 + z * r * f11;
}

/*
 * Calculate overlap between two grids (int grida^*gridb; cylindrical coords).
 *
 * grida = 1st grid (complex conjugated) (cgrid2d *).
 * gridb = 2nd grid (no complex conjugation) (cgrid2d *).
 *
 * Returns the value of the overlap integral (double complex).
 *
 */

EXPORT double complex cgrid2d_integral_of_conjugate_product_cyl(const cgrid2d *grida, const cgrid2d *gridb) {

  long j, ij, nx = grida->nx, ny = grida->ny, nxy = nx * ny;
  double step = grida->step, r;
  double complex sum = 0.0;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,ny,nxy,step,avalue,bvalue) private(ij,j,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    j = ij % ny;    
    r = step * (double) j;
    sum += conj(avalue[ij]) * bvalue[ij] * r;
  }

  return sum * step * step * 2.0 * M_PI;
}

/*
 * Extrapolate between two different grid sizes.
 *
 * dest = Destination grid (cgrid2d *; output).
 * src  = Source grid (cgrid2d *; input).
 *
 */

EXPORT void cgrid2d_extrapolate_cyl(cgrid2d *dest, cgrid2d *src) {

  long i, j, nx = dest->nx, ny = dest->ny;
  double step = dest->step, z, r, offset = step * src->ny/2.0;

  for (i = 0; i < nx; i++) {
    z = (i - nx/2.0) * step;
    for (j = 0; j < ny; j++) {
      r = j * step - offset;
      dest->value[i * ny + j] = cgrid2d_value(src, z, r);
    }
  }
}
