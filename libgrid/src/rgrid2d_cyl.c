/*
 * Routines for 2D real cylindrical grids.
 *
 * NOTE: Coordinates are (z,r).
 *
 */

#include "grid.h"
#include "private.h"
#include "private2d.h"

/*
 * Multiply a given grid by a function (cylindrical coordinates).
 *
 * grid = destination grid for the operation (rgrid2d *).
 * func = function providing the mapping (double (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and z,r are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_product_func_cyl(rgrid2d *grid, double (*func)(void *arg, double z, double r), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double z, r, step = grid->step;
  double *value = grid->value;
  
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
 * grid = destination grid for the operation (rgrid2d *).
 * func = function providing the mapping (double (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and z,r are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_map_cyl(rgrid2d *grid, double (*func)(void *arg, double z, double r), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double z, r, step = grid->step;
  double *value = grid->value;
  
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
 * grid = destination grid for the operation (rgrid2d *).
 * func = function providing the mapping (double (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and (z,r) is the point (doubles) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 * ns   = number of intermediate points to be used in smoothing (int).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_smooth_map_cyl(rgrid2d *grid, double (*func)(void *arg, double z, double r), void *farg, int ns) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double zc, rc, step = grid->step;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,ns,step,func,value) private(i,j,zc,rc) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;    
    zc = (i - nx/2.0) * step;
    rc = j * step;
    value[ij] = linearly_weighted_integralr2d(func, farg, zc, rc, step, ns);
  }
}

/*
 * Map a given function onto 2D cylindrical grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid. Limits for intermediate steps and
 * tolerance can be given.
 *
 * grid   = destination grid for the operation (rgrid2d *).
 * func   = function providing the mapping (double (*)(void *, double, double, double)).
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

EXPORT void rgrid2d_adaptive_map_cyl(rgrid2d *grid, double (*func)(void *arg, double z, double r), void *farg, int min_ns, int max_ns, double tol) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny, ns;
  double zc, rc, step = grid->step;
  double tol2 = tol * tol;
  double sum = 0.0, sump;
  double *value = grid->value;
  
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
      sum  = linearly_weighted_integralr2d(func, farg, zc, rc, step, ns);
      sump = linearly_weighted_integralr2d(func, farg, zc, rc, step, ns+1);
      if ((sum - sump)*(sum - sump) < tol2) break;
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
 * Integrate over a cylindrical grid.
 *
 * grid = grid to be integrated (rgrid2d *).
 *
 * Returns the integral value (double).
 *
 */

EXPORT double rgrid2d_integral_cyl(const rgrid2d *grid) {

  long ij, j, ny = grid->ny, nxy = grid->nx * ny;
  double sum = 0.0;
  double *value = grid->value, r, step = grid->step;
  
#pragma omp parallel for firstprivate(ny,nxy,value,step) private(ij,r,j) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    j = ij % ny;
    r = step * (double) j;
    sum += value[ij] * r;
  }

  return sum * step * step * 2.0 * M_PI;
}

/*
 * Integrate over cylindrical grid with limits.
 *
 * grid = grid to be integrated (rgrid2d *).
 * zl   = lower limit for z (double).
 * zu   = upper limit for z (double).
 * rl   = lower limit for r (double).
 * ru   = upper limit for r (double).
 *
 * Returns the integral value (double).
 *
 */

EXPORT double rgrid2d_integral_region_cyl(const rgrid2d *grid, double zl, double zu, double rl, double ru) {

  long iu, il, i, ju, jl, j, nx = grid->nx, ny = grid->ny;
  double *value = grid->value, sum;
  double step = grid->step, r;
   
  il = zl / step + nx/2;
  iu = zu / step + nx/2;
  jl = rl / step;
  ju = ru / step;
  
  sum = 0.0;
#pragma omp parallel for firstprivate(il,iu,jl,ju,nx,ny,value,step) private(i,j,r) reduction(+:sum) default(none) schedule(runtime)
  for (i = il; i < iu; i++)
    for (j = jl; j < ju; j++) {
      r = step * (double) j;
      sum += value[i * ny  + j] * r;
    }
  return sum * step * step * 2.0 * M_PI;
}

/* 
 * Integrate over cylindrical grid squared (int grid^2).
 *
 * grid = grid to be integrated (rgrid2d *).
 *
 * Returns the integral (double).
 *
 */

EXPORT double rgrid2d_integral_of_square_cyl(const rgrid2d *grid) {

  long j, ij, nr = grid->ny, nzr = grid->nx * nr;
  double sum = 0.0, r;
  double *value = grid->value, step = grid->step;;
  
#pragma omp parallel for firstprivate(nr,nzr,value,step) private(ij,j,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nzr; ij++) {
    j = ij % nr;
    r = step * (double) j;
    sum += value[ij] * value[ij] * r;
  }
  
  return sum * step * step * 2.0 * M_PI;
}

/*
 * Calculate overlap between two cylindrical 2D grids (int grida gridb).
 *
 * grida = 1st grid (rgrid3d *).
 * gridb = 2nd grid (rgrid3d *).
 *
 * Returns the value of the overlap integral (double).
 *
 */

EXPORT double rgrid2d_integral_of_product_cyl(const rgrid2d *grida, const rgrid2d *gridb) {

  long j, ij, nzr = grida->nx * grida->ny, nr = grida->ny;
  double sum = 0.0, r;
  double *avalue = grida->value, *bvalue = gridb->value, step = grida->step;
  
#pragma omp parallel for firstprivate(nzr,avalue,bvalue,step,nr) private(ij,j,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nzr; ij++) {
    j = ij % nr;
    r = step * (double) j;
    sum += r * avalue[ij] * bvalue[ij];
  }

  return sum * grida->step * grida->step * 2.0 * M_PI;
}

/*
 * Calculate expectation value of cylindrical grid over the probability density given by the other.
 * (int gridb grida^2).
 *
 * grida = grid giving the probability (grida^2) (rgrid2d *).
 * gridb = grid to be averaged (rgrid2d *).
 *
 * Returns the average value (double).
 *
 */

EXPORT double rgrid2d_grid_expectation_value_cyl(const rgrid2d *grida, const rgrid2d *gridb) {

  long ij, j, ny = grida->ny, nxy = grida->nx * ny;
  double sum = 0.0, step = grida->step;
  double *avalue = grida->value, r;
  double *bvalue = gridb->value;
  
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
 * func  = function to be averaged (double (*)(void *, double, double, double)).
 *         The arguments are: optional arg, grida(z,r), z, r.
 * grida = grid giving the probability (grida^2) (rgrid2d *).
 *
 * Returns the average value (double).
 *
 */
 
EXPORT double rgrid2d_grid_expectation_value_func_cyl(void *arg, double (*func)(void *arg, double val, double z, double r), const rgrid2d *grida) {
   
  long ij, i, j, nx = grida->nx, ny = grida->ny, nxy = nx * ny;
  double sum = 0.0;
  double *avalue = grida->value;
  double z, r, step = grida->step;
  
#pragma omp parallel for firstprivate(func,arg,nxy,nx,ny,avalue,step) private(z,r,i,j,ij) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    z = (i - nx/2.0) * step;
    r = j * step;
    sum += avalue[ij] * avalue[ij] * func(arg, avalue[ij], z, r) * r;
  }
  
  return sum * step * step * 2.0 * M_PI;
}


/* 
 * Integrate over cylindrical grid multiplied by weighting function (int grid w(z,r)).
 *
 * grid   = grid to be integrated over (rgrid2d *).
 * weight = function defining the weight (double (*)(double, double)). The arguments are (z,r) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double rgrid2d_weighted_integral_cyl(const rgrid2d *grid, double (*weight)(void *farg, double z, double r), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, ny = grid->ny, nx = grid->nx;
  double z, r, step = grid->step;
  double sum = 0.0;
  double w;
  double *value = grid->value;
  
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
 * grid   = grid to be integrated over (rgrid2d *).
 * weight = function defining the weight (double (*)(double, double)).
 *          The arguments are (z,r) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double rgrid2d_weighted_integral_of_square_cyl(const rgrid2d *grid, double (*weight)(void *farg, double z, double r), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double z, r, step = grid->step;
  double sum = 0.0, w;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nx,ny,step,value,weight,farg) private(w,i,j,ij,z,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {    
    i = ij / ny;
    j = ij % ny;    
    z = (i - nx/2.0) * step;
    r = j * step;
    w = weight(farg, z, r);
    sum += w * value[ij] * value[ij] * r;
  }
  
  return sum * step * step * 2.0 * M_PI;
}

/* 
 * Differentiate a grid with respect to z.
 * This is used in dot product: grad . f along z.  (just a derivative with respect to z)
 *
 * grid     = grid to be differentiated (rgrid2d *).
 * gradient = differentiated grid output (rgrid2d *).
 * 
 * No return value.
 *
 */

EXPORT void rgrid2d_fd_gradient_cyl_z(const rgrid2d *grid, rgrid2d *gradient) {

  rgrid2d_fd_gradient_x(grid, gradient);
}

/* 
 * Differentiate a grid with respect to r.
 * This is used in dot product: grad . f along r.  (just a derivative with respect to r)
 *
 * grid     = grid to be differentiated (rgrid2d *).
 * gradient = differentiated grid output (rgrid2d *).
 * 
 * No return value.
 *
 */

EXPORT void rgrid2d_fd_gradient_cyl_r(const rgrid2d *grid, rgrid2d *gradient) {

  rgrid2d_fd_gradient_y(grid, gradient);
}

#if 0
//THIS WAS USED for grad . A but we don't really need it anywhere??
// And this was certainly wrong for calculating the plain gradient...

/* 
 * Differentiate a grid with respect to r, i.e. (1/r) * d(rf(r))/dr = df(r)/dr + (1/r)f(r).
 * This used in dot product: grad . f along r. Note that gradient in cylindrical coordinates
 * does NOT have the additional 1/r term.
 *
 * grid     = grid to be differentiated (rgrid2d *).
 * gradient = differentiated grid output (rgrid2d *).
 * 
 * No return value.
 *
 */

EXPORT void rgrid2d_fd_gradient_cyl_r(const rgrid2d *grid, rgrid2d *gradient) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double step = grid->step;
  double inv_delta = 1.0 / (2.0 * step);
  double *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in rgrid2d_fd_gradient_cyl_r().\n");
    return;
  }
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta,grid,step) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = inv_delta * (rgrid2d_value_at_index(grid, i, j+1) - rgrid2d_value_at_index(grid, i, j-1));
  }
}
#endif


/*
 * Calculate gradient of a cylindrical grid.
 *
 * grid       = grid to be differentiated twice (rgrid2d *).
 * gradient_z = z output grid for the operation (rgrid2d *).
 * gradient_r = r output grid for the operation (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_fd_gradient_cyl(const rgrid2d *grid, rgrid2d *gradient_z, rgrid2d *gradient_r) {

  rgrid2d_fd_gradient_cyl_z(grid, gradient_z);
  rgrid2d_fd_gradient_cyl_r(grid, gradient_r);
}

/*
 * Calculate laplacian of cylindrical grid.
 *
 * grid    = source grid (rgrid2d *).
 * laplace = output grid for the operation (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_fd_laplace_cyl(const rgrid2d *grid, rgrid2d *laplace) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double step = grid->step, inv_delta2 = 1.0 / (step * step);
  double inv_delta = 1.0 / (2.0 * step);
  double *lvalue = laplace->value, r;
  
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta2,inv_delta,grid,step) private(ij,i,j,r) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    r = step * (double) j;
    lvalue[ij] = inv_delta2 * (-4.0 * rgrid2d_value_at_index(grid, i, j) + rgrid2d_value_at_index(grid, i, j+1) + rgrid2d_value_at_index(grid, i, j-1) + rgrid2d_value_at_index(grid,i+1,j) + rgrid2d_value_at_index(grid,i-1,j));
    if(r != 0.0)
      lvalue[ij] += inv_delta * (rgrid2d_value_at_index(grid, i, j+1) - rgrid2d_value_at_index(grid, i, j-1)) / r;
  }
}

/*
 * Access grid point at given (z,r) point (with linear interpolation).
 *
 * grid = grid to be accessed (rgrid2d *).
 * x    = x value (double).
 * y    = y value (double).
 *
 * Returns grid value at (x,y).
 *
 */

EXPORT inline double rgrid2d_value_cyl(const rgrid2d *grid, double z, double r) {

  double f00, f10, f01, f11;
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
  f00 = rgrid2d_value_at_index(grid, i, j);
  f10 = rgrid2d_value_at_index(grid, i+1, j);
  f01 = rgrid2d_value_at_index(grid, i, j+1);
  f11 = rgrid2d_value_at_index(grid, i+1, j+1);
  
  omz = 1 - z;
  omr = 1 - r;

  return omz * omr * f00 + z * omr * f10 + omz * r * f01 + z * r * f11;
}

/*
 * Extrapolate between two different grid sizes.
 *
 * dest = Destination grid (rgrid2d *; output).
 * src  = Source grid (rgrid2d *; input).
 *
 */

EXPORT void rgrid2d_extrapolate_cyl(rgrid2d *dest, rgrid2d *src) {

  long i, j, nx = dest->nx, ny = dest->ny;
  double step = dest->step, z, r, offset = step * src->ny/2.0;

  for (i = 0; i < nx; i++) {
    z = (i - nx/2.0) * step;
    for (j = 0; j < ny; j++) {
      r = j * step - offset;
      dest->value[i * ny + j] = rgrid2d_value(src, z, r);
    }
  }
}
