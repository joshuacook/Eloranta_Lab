/*
 * Routines for 3D complex cylindrical grids.
 *
 */

#include "grid.h"
#include "private.h"
#include "private3d.h"

/*
 * Access grid point at given index. 3D cylindrical coordinates.
 * Note that the bondaries are fixed: z and phi are periodic
 * whereas r is not.
 *
 * grid = grid to be accessed (cgrid3d *).
 * i    = index along x (long).
 * j    = index along y (long).
 * k    = index along z (long).
 *
 * Returns grid value at index (i, j, k).
 *
 */

EXPORT inline double complex cgrid3d_value_at_index_cyl(const cgrid3d *grid, long i, long j, long k) {

  long nr = grid->nx, nphi = grid->ny, nz = grid->nz;

  if (i < 0 || j < 0 || k < 0 || i >= nr || j >= nphi || k >= nz) {
    if(i < 0) {
      i = abs(i);
      j += nphi/2;
    }
    if(i >= nr) i = nr-1;
    j %= nphi;
    if (j < 0) j = nphi + j;
    k %= nz;
    if (k < 0) k = nz + k;  
  }
  return grid->value[(i*nphi + j)*nz + k];  
}

/*
 * Multiply a given cylindrical grid by a function.
 *
 * grid = destination grid for the operation (cgrid3d *).
 * func = function providing the mapping (double complex (*)(void *, double, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and (r, phi, z) are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_product_func_cyl(cgrid3d *grid, double complex (*func)(void *arg, double r, double phi, double z), void *farg) {

  long i, j, k, ij, ijnz, nr = grid->nx, nphi = grid->ny, nz = grid->nz;
  double z, r, phi, step = grid->step, step_phi = 2.0 * M_PI / (double) nphi;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nz,nr,nphi,step,step_phi,func,value) private(i,j,ij,ijnz,k,z,r,phi) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    j = ij % nphi;
    r = i * step;    
    phi = j * step_phi;
    for(k = 0; k < nz; k++) {
      z = (k - nz/2.0) * step;
      value[ijnz + k] *= func(farg, r, phi, z);
    }
  }
}

/*
 * Map a given function onto 3D cylindrical grid.
 *
 * grid = destination grid for the operation (cgrid3d *).
 * func = function providing the mapping (double complex (*)(void *, double, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and (r, phi, z) are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_map_cyl(cgrid3d *grid, double complex (*func)(void *arg, double r, double phi, double z), void *farg) {

  long i, j, k, ij, ijnz, nr = grid->nx, nphi = grid->ny, nz = grid->nz;
  double z, r, phi, step = grid->step, step_phi = 2.0 * M_PI / (double) nphi;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nz,nr,nphi,step,step_phi,func,value) private(i,j,ij,ijnz,k,z,r,phi) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    j = ij % nphi;
    r = i * step;
    phi = j * step_phi;    
    for(k = 0; k < nz; k++) {
      z = (k - nz/2.0) * step;
      value[ijnz + k] = func(farg, r, phi, z);
    }
  }
}

/*
 * Integrate over a 3D cylindrical grid.
 *
 * grid = grid to be integrated (cgrid3d *).
 *
 * Returns the integral value (double).
 *
 */

EXPORT double complex cgrid3d_integral_cyl(const cgrid3d *grid) {

  long i, ij, k, ijnz, nr = grid->nx, nphi = grid->ny, nz = grid->nz;
  double sum = 0.0, ssum, r, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nr,nz,nphi,value,step) private(i,ij,ijnz,k,ssum,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    r = i * step;
    ssum = 0.0;
    for(k = 0; k < nz; k++)
      ssum += value[ijnz + k];
    sum += r * ssum;
  }
  
  return sum * grid->step * grid->step * (2.0 * M_PI / (double) nphi);
}

/*
 * Integrate over 3D cylindrical grid with limits.
 *
 * grid = grid to be integrated (cgrid3d *).
 * rl   = lower limit for r (double).
 * ru   = upper limit for r (double).
 * phil = lower limit for phi (double).
 * phiu = upper limit for phi (double).
 * zl   = lower limit for z (double).
 * zu   = upper limit for z (double).
 *
 * Returns the integral value (double).
 *
 */

EXPORT double complex cgrid3d_integral_region_cyl(const cgrid3d *grid, double rl, double ru, double phil, double phiu, double zl, double zu) {

  long iu, il, i, ju, jl, j, ku, kl, k, nphi = grid->ny, nz = grid->nz;
  double r, sum;
  double complex *value = grid->value;
  double step = grid->step, step_phi = 2.0 * M_PI / (double) nphi;
   
  il = rl / step;
  iu = ru / step;
  jl = phil / step_phi;
  ju = phiu / step_phi;
  kl = zl / step + nz/2;
  ku = zu / step + nz/2;
  
  sum = 0.0;
#pragma omp parallel for firstprivate(il,iu,jl,ju,kl,ku,nz,nphi,value,step) private(i,j,k,r) reduction(+:sum) default(none) schedule(runtime)
  for (i = il; i < iu; i++) {
    r = step * (double) i;
    for (j = jl; j < ju; j++)
      for (k = kl; k < ku; k++)
	sum += value[(i * nphi + j) * nz  + k] * r;
  }
  return sum * step * step * step_phi; 
}

/* 
 * Integrate over the cylindrical grid squared (int grid^2).
 *
 * grid = grid to be integrated (cgrid3d *).
 *
 * Returns the integral (double).
 *
 */

EXPORT double cgrid3d_integral_of_square_cyl(const cgrid3d *grid) {

  long i, ij, k, ijnz, nr = grid->nx, nphi = grid->ny, nz = grid->nz;
  double sum = 0, ssum = 0, r, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nphi,value,nr,nz,step) private(ij,ijnz,k,ssum,i,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    r = i * step;
    ssum = 0.0;
    for(k = 0; k < nz; k++)
      ssum += sqnorm(value[ijnz + k]);
    sum += ssum * r;
  }
  
  return sum * grid->step * grid->step * (2.0 * M_PI / (double) nphi);
}

/*
 * Calculate overlap between two cylindrical 3D grids (int grida gridb).
 *
 * grida = 1st grid (cgrid3d *).
 * gridb = 2nd grid (cgrid3d *).
 *
 * Returns the value of the overlap integral (double).
 *
 */

EXPORT double complex cgrid3d_integral_of_conjugate_product_cyl(const cgrid3d *grida, const cgrid3d *gridb) {

  long i, ij, k, ijnz, nr = grida->nx, nphi = grida->ny, nz = grida->nz;
  double sum = 0.0, ssum, r, step = grida->step;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nz,nphi,avalue,bvalue,nr,step) private(ij,i,ijnz,k,ssum,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    r = i * step;
    ssum = 0.0;
    for(k = 0; k < nz; k++)
      ssum += conj(avalue[ijnz + k]) * bvalue[ijnz + k];
    sum += r * ssum;
  }

  return sum * grida->step * grida->step * (2.0 * M_PI / (double) nphi);
}

/*
 * Calculate the expectation value of a grid over a 3D cylindrical grid.
 * (int gridb grida gridb = int grida gridb^2).
 *
 * grida = grid giving the probability (gridb^2) (cgrid3d *).
 * gridb = grid to be averaged (cgrid3d *).
 *
 * Returns the average value (double *).
 *
 */

EXPORT double complex cgrid3d_grid_expectation_value_cyl(const cgrid3d *grida, const cgrid3d *gridb) {

  long i, ij, k, ijnz, nr = grida->nx, nphi = grida->ny, nz = grida->nz;
  double sum = 0.0, ssum, r, step = grida->step;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nz,nphi,avalue,bvalue,step,nr) private(ij,ijnz,k,ssum,i,r) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    r = i * step;
    ssum = 0.0;
    for(k = 0; k < nz; k++)
      ssum += sqnorm(avalue[ijnz + k]) * bvalue[ijnz + k];
    sum += r * ssum;
  }
  
  return sum * grida->step * grida->step * (2.0 * M_PI / (double) nphi);
}

/*
 * Calculate the expectation value of a function over a grid.
 * (int grida func grida = int func grida^2).
 *
 * func  = function to be averaged (double complex (*)(void *, double, double, double, double)).
 *         The arguments are: optional arg, grida(r,phi,z), r, phi, z.
 * grida = grid giving the probability (grida^2) (cgrid3d *).
 *
 * Returns the average value (double).
 *
 */
 
EXPORT double complex cgrid3d_grid_expectation_value_func_cyl(void *arg, double complex (*func)(void *arg, double val, double r, double phi, double z), const cgrid3d *grida) {
   
  long ij, i, j, k, ijnz, nr = grida->nx, nphi = grida->ny, nz = grida->nz;
  double sum = 0.0, ssum;
  double complex *avalue = grida->value;
  double z, r, phi, step = grida->step, step_phi = 2.0 * M_PI / (double) nphi;
 
#pragma omp parallel for firstprivate(func,arg,nz,nr,nphi,avalue,step,step_phi) private(z,r,phi,i,j,ij,ijnz,k,ssum) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    j = ij % nphi;
    r = i * step;
    phi = j * step_phi;
    ssum = 0.0;
    for(k = 0; k < nz; k++) {
      z = (k - nz/2.0) * step;
      ssum += sqnorm(avalue[ijnz + k]) * func(arg, avalue[ijnz + k], r, phi, z);
    }
    sum += r * ssum;
  }
  
  return sum * step * step * step_phi;
}

/* 
 * Integrate over the 3D cylindrical grid multiplied by weighting function (int grid w(x)).
 *
 * grid   = grid to be integrated over (cgrid3d *).
 * weight = function defining the weight (double complex (*)(double, double, double)). The arguments are (r,phi,z) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double complex cgrid3d_weighted_integral_cyl(const cgrid3d *grid, double complex (*weight)(void *farg, double r, double phi, double z), void *farg) {

  long i, j, k, ij, ijnz, nz = grid->nz, nr = grid->nx, nphi = grid->ny;
  double z, r, phi, step = grid->step, step_phi = 2.0 * M_PI / (double) nphi;
  double sum = 0.0, ssum;
  double w;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nz,nr,nphi,step,value,weight,farg,step_phi) private(w,i,j,k,ij,ijnz,z,r,phi,ssum) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    j = ij % nphi;
    r = i * step;
    phi = j * step_phi;
    ssum = 0.0;
    for(k = 0; k < nz; k++) {
      z = (k - nz/2.0) * step;
      w = weight(farg, r, phi, z);
      ssum += w * value[ijnz + k];
    }
    sum += r * ssum;
  }
  
  return sum * step * step * step_phi;
}

/* 
 * Integrate over square of the 3D cylindrical grid multiplied by weighting function (int grid^2 w(x)).
 *
 * grid   = grid to be integrated over (cgrid3d *).
 * weight = function defining the weight (double complex (*)(double, double, double)).
 *          The arguments are (r,phi,z) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double cgrid3d_weighted_integral_of_square_cyl(const cgrid3d *grid, double complex (*weight)(void *farg, double r, double phi, double z), void *farg) {

  long i, j, k, ij, ijnz, nr = grid->nx, nphi = grid->ny, nz = grid->nz;
  double z, r, phi, step = grid->step, step_phi = 2.0 * M_PI / (double) nphi;
  double sum = 0, ssum, w;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nz,nr,nphi,step,value,weight,farg,step_phi) private(w,i,j,ij,ijnz,k,z,r,phi,ssum) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    j = ij % nphi;
    r = i * step;
    phi = j * step_phi;      
    ssum = 0.0;
    for(k = 0; k < nz; k++) {
      z = (k - nz/2.0) * step;
      w = weight(farg, r, phi, z);
      ssum += w * sqnorm(value[ijnz + k]);
    }    
    sum += r * ssum;
  }
  
  return sum * step * step * step_phi;
}

/* 
 * Differentiate a grid with respect to r.
 *
 * grid     = grid to be differentiated (cgrid3d *).
 * gradient = differentiated grid output (cgrid3d *).
 * 
 * No return value.
 *
 */

EXPORT void cgrid3d_fd_gradient_cyl_r(const cgrid3d *grid, cgrid3d *gradient) {

  cgrid3d_fd_gradient_x(grid, gradient);
}

/* 
 * Differentiate a grid with respect to phi.
 *
 * grid     = grid to be differentiated (cgrid3d *).
 * gradient = differentiated grid output (cgrid3d *).
 * 
 * No return value.
 *
 */

static double complex tmp_one_over_r(void *NA, double r, double phi, double z) {

  return 1.0 / (r + GRID_EPS);
}

EXPORT void cgrid3d_fd_gradient_cyl_phi(const cgrid3d *grid, cgrid3d *gradient) {

  cgrid3d_fd_gradient_y(grid, gradient);
  cgrid3d_product_func_cyl(gradient, &tmp_one_over_r, NULL);
}

/* 
 * Differentiate a grid with respect to z.
 *
 * grid     = grid to be differentiated (cgrid3d *).
 * gradient = differentiated grid output (cgrid3d *).
 * 
 * No return value.
 *
 */

EXPORT void cgrid3d_fd_gradient_cyl_z(const cgrid3d *grid, cgrid3d *gradient) {

  cgrid3d_fd_gradient_z(grid, gradient);
}

/*
 * Calculate gradient of a cylindrical 3D grid.
 *
 * grid         = grid to be differentiated twice (cgrid3d *).
 * gradient_r   = r output grid for the operation (cgrid3d *).
 * gradient_phi = phi output grid for the operation (cgrid3d *).
 * gradient_z   = z output grid for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_fd_gradient_cyl(const cgrid3d *grid, cgrid3d *gradient_r, cgrid3d *gradient_phi, cgrid3d *gradient_z) {

  cgrid3d_fd_gradient_cyl_r(grid, gradient_r);
  cgrid3d_fd_gradient_cyl_phi(grid, gradient_phi);
  cgrid3d_fd_gradient_cyl_z(grid, gradient_z);
}

/*
 * Calculate laplacian of cylindrical 3D grid.
 *
 * grid    = source grid (cgrid3d *).
 * laplace = output grid for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_fd_laplace_cyl(const cgrid3d *grid, cgrid3d *laplace) {

  long i, j, k, ij, ijnz, nr = grid->nx, nphi = grid->ny, nz = grid->nz;
  double step = grid->step, step_phi = 2.0 * M_PI / (double) nphi, inv_phi_step2 = 1.0 / (step_phi * step_phi), inv_step2 = 1.0 / (step * step), inv_2step = 1.0 / (2.0 * step), r;
  double complex *lvalue = laplace->value;
  
#pragma omp parallel for firstprivate(nz,nr,lvalue,inv_step2,inv_2step,inv_phi_step2,grid,nphi,step) private(ij,ijnz,i,j,k,r) default(none) schedule(runtime)
  for(ij = 0; ij < nr * nphi; ij++) {
    ijnz = ij * nz;
    i = ij / nphi;
    j = ij % nphi;
    r = step * (double) i;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = inv_step2 * (cgrid3d_value_at_index_cyl(grid, i, j, k-1) - 2.0 * cgrid3d_value_at_index_cyl(grid, i, j, k) - cgrid3d_value_at_index_cyl(grid, i, j, k+1)) /* d2/dz^2 */
	+ inv_step2 * (cgrid3d_value_at_index_cyl(grid, i-1, j, k) - 2.0 * cgrid3d_value_at_index_cyl(grid, i, j, k) - cgrid3d_value_at_index_cyl(grid, i+1, j, k)) /* d2/dr^2 */
	+ (inv_2step / (r + GRID_EPS)) * (cgrid3d_value_at_index_cyl(grid, i+1, j, k) - cgrid3d_value_at_index_cyl(grid, i-1, j, k)) /* (1/r) d/dr */
	+ (inv_phi_step2 / (r * r + GRID_EPS)) * (cgrid3d_value_at_index_cyl(grid, i, j-1, k) - 2.0 * cgrid3d_value_at_index_cyl(grid, i, j, k) - cgrid3d_value_at_index_cyl(grid, i, j+1, k)); /* (1/r^2) * d2/dphi^2 */
  }
}

/*
 * Access grid point at given (r,phi,z) point (with linear interpolation).
 *
 * grid = grid to be accessed (cgrid3d *).
 * r    = r value (double).
 * phi  = phi value (double).
 * z    = z value (double).
 *
 * Returns grid value at (r,phi,z).
 *
 * (the most positive elements are missing - rolled over to negative; 
 *  periodic boundaries)
 *
 */

EXPORT inline double complex cgrid3d_value_cyl(const cgrid3d *grid, double r, double phi, double z) {

  double complex f000, f100, f010, f001, f110, f101, f011, f111;
  double omz, omr, omphi, step = grid->step, step_phi = 2.0 * M_PI / (double) grid->ny;
  long i, j, k;
  
  /* i to index and 0 <= r < 1 */
  i = (r /= step);
  r -= i;
  
  /* j to index and 0 <= phi < 1 */
  j = (phi /= step_phi);
  phi -= j;

  /* k to index and 0 <= z < 1 */
  k = (z /= step);
  if (z < 0) k--;
  z -= k;
  k += grid->nz / 2;
  
  /* linear extrapolation 
   *
   * f(r,phi,z) = (1-r) (1-phi) (1-z) f(0,0,0) + r (1-phi) (1-z) f(1,0,0) + (1-r) phi (1-z) f(0,1,0) + (1-r) (1-phi) z f(0,0,1) 
   *            + r       phi   (1-z) f(1,1,0) + r (1-phi)   z   f(1,0,1) + (1-r) phi   z   f(0,1,1) +   r     phi   z f(1,1,1)
   */ 
  f000 = cgrid3d_value_at_index_cyl(grid, i, j, k);
  f100 = cgrid3d_value_at_index_cyl(grid, i+1, j, k);
  f010 = cgrid3d_value_at_index_cyl(grid, i, j+1, k);
  f001 = cgrid3d_value_at_index_cyl(grid, i, j, k+1);
  f110 = cgrid3d_value_at_index_cyl(grid, i+1, j+1, k);
  f101 = cgrid3d_value_at_index_cyl(grid, i+1, j, k+1);
  f011 = cgrid3d_value_at_index_cyl(grid, i, j+1, k+1);
  f111 = cgrid3d_value_at_index_cyl(grid, i+1, j+1, k+1);
  
  omr = 1.0 - r;
  omphi = 1.0 - phi;
  omz = 1.0 - z;

  return omr * omphi * omz * f000 + r * omphi * omz * f100 + omr * phi * omz * f010 + omr * omphi * z * f001
    + r * phi * omz * f110 + r * omphi * z * f101 + omr * phi * z * f011 + r * phi * z * f111;
}

/*
 * Extrapolate between two different grid sizes (3D cylindrical).
 *
 * dest = Destination grid (cgrid3d *; output).
 * src  = Source grid (cgrid3d *; input).
 *
 */

EXPORT void cgrid3d_extrapolate_cyl(cgrid3d *dest, cgrid3d *src) {

  long i, j, k, nr = dest->nx, nphi = dest->ny, nz = dest->nz, nphiz = nphi * nz, tmp;
  double step = dest->step, step_phi = 2.0 * M_PI / (double) nphi, z, r, phi;

  for (i = 0; i < nr; i++) {
    r = i * step;
    for (j = 0; j < nphi; j++) {
      phi = j * step_phi;
      tmp = i * nphiz + j * nz;
      for (k = 0; k < nz; k++) {
	z = (k - nz/2.0) * step;
	dest->value[tmp + k] = cgrid3d_value_cyl(src, r, phi, z);
      }
    }
  }
}

/*
 * Map cylindrical grid onto Cartesian grid (3D).
 *
 * Note that the Cartesian grid must have the following dimensions:
 *
 * Cyl:(r, phi, z) -> Cart:(x, y, z)
 *      nx, ny, nz         nx, ny, nz
 * 
 * Requirement: cyl->nz == cart->nz, 2 * cyl->nx - 1 == cart->nx, 2 * cyl->nx - 1 == cart->ny.
 *
 */

EXPORT void cgrid3d_map_cyl_on_cart(cgrid3d *cart, cgrid3d *cyl) {

  long nx = cart->nx, ny = cart->ny, nz = cart->nz, i, j, k, nyz = ny * nz;
  long nr = cyl->nx;
  double z, r, phi, x, y, step = cart->step;

  if(nz != cyl->nz || (2 * nr) != nx || (2 * nr) != ny) {
    fprintf(stderr, "libgrid: Illegal dimensions in cgrid3d_map_cyl_on_cart().\n");
    exit(1);
  }
  for(i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    for (j = 0; j < ny; j++) {
      y = (j - ny/2.0) * step;
      for (k = 0; k < nz; k++) {
	z = (k - nz/2.0) * step;
	r = sqrt(x*x + y*y);
	// phi: [0,2\pi[ 
	phi = M_PI - atan2(y, -x);
	cart->value[i * nyz + j * nz + k] = cgrid3d_value_cyl(cyl, r, phi, z);
      }
    }
  }
}

/*
 * Map Cartesian grid onto Cylindrical grid (3D).
 *
 */

EXPORT void cgrid3d_map_cart_on_cyl(cgrid3d *cyl, cgrid3d *cart) {

  long i, j, k, nz = cyl->nz, nr = cyl->nx, nphi = cyl->ny, nphiz = nphi * nz;
  long ny = cart->ny, nx = cart->nx;
  double x, y, z, r, phi, step = cyl->step, step_phi = 2.0 * M_PI / (double) nphi;

  if(nz != cart->nz || (2 * nr) != nx || (2 * nr) != ny) {
    fprintf(stderr, "libgrid: Illegal dimensions in cgrid3d_map_cart_on_cyl().\n");
    exit(1);
  }

  for (i = 0; i < nr; i++) {
    r = step * (double) i;
    for (j = 0; j < nphi; j++) {
      phi = step_phi * (double) j;
      for (k = 0; k < nz; k++) {
	z = (k - nz/2.0) * step;
	x = r * cos(phi);
	y = r * sin(phi);
	cyl->value[i * nphiz + j * nz + k] = cgrid3d_value(cart, x, y, z);
      }
    }
  }
}

/*
 * Add cylindrical grid to Cartesian grid (3D).
 *
 * Note that the Cartesian grid must have the following dimensions:
 *
 * Cyl:(r, phi, z) -> Cart:(x, y, z)
 *      nx, ny, nz         nx, ny, nz
 * 
 * Requirement: cyl->nz == cart->nz, 2 * cyl->nx - 1 == cart->nx, 2 * cyl->nx - 1 == cart->ny.
 *
 */

EXPORT void cgrid3d_add_cyl_on_cart(cgrid3d *cart, cgrid3d *cyl) {

  long nx = cart->nx, ny = cart->ny, nz = cart->nz, i, j, k, nyz = ny * nz;
  long nr = cyl->nx;
  double z, r, phi, x, y, step = cart->step;

  if(nz != cyl->nz || (2 * nr) != nx || (2 * nr) != ny) {
    fprintf(stderr, "libgrid: Illegal dimensions in rgrid3d_add_cyl_on_cart().\n");
    exit(1);
  }
  for(i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    for (j = 0; j < ny; j++) {
      y = (j - ny/2.0) * step;
      for (k = 0; k < nz; k++) {
	z = (k - nz/2.0) * step;
	r = sqrt(x*x + y*y);
	// phi: [0,2\pi[ 
	phi = M_PI - atan2(y, -x);
	cart->value[i * nyz + j * nz + k] += cgrid3d_value_cyl(cyl, r, phi, z);
      }
    }
  }
}

/*
 * Add Cartesian grid onto Cylindrical grid (3D).
 *
 */

EXPORT void cgrid3d_add_cart_on_cyl(cgrid3d *cyl, cgrid3d *cart) {

  long i, j, k, nz = cyl->nz, nr = cyl->nx, nphi = cyl->ny, nphiz = nphi * nz;
  long ny = cart->ny, nx = cart->nx;
  double x, y, z, r, phi, step = cyl->step, step_phi = 2.0 * M_PI / (double) nphi;

  if(nz != cart->nz || (2 * nr) != nx || (2 * nr) != ny) {
    fprintf(stderr, "libgrid: Illegal dimensions in rgrid3d_add_cart_on_cyl().\n");
    exit(1);
  }

  for (i = 0; i < nr; i++) {
    r = step * (double) i;
    for (j = 0; j < nphi; j++) {
      phi = step_phi * (double) j;
      for (k = 0; k < nz; k++) {
	z = (k - nz/2.0) * step;
	x = r * cos(phi);
	y = r * sin(phi);
	cyl->value[i * nphiz + j * nz + k] += cgrid3d_value(cart, x, y, z);
      }
    }
  }
}
