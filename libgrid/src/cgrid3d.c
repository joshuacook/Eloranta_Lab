 /*
 * Routines for 3D complex grids.
 *
 * Nx is major index and Nz is minor index (varies most rapidly).
 *
 * To debug memory allocation, define __DEBUG__.
 *
 */

#include "grid.h"
#include "private.h"
#include "private3d.h"

#ifdef __DEBUG__
static int allocated_grids = 0;
#endif /* __DEBUG__ */

/*
 * Allocate 3D grid.
 *
 * nx                 = number of points on the grid along x (long).
 * ny                 = number of points on the grid along y (long).
 * nz                 = number of points on the grid along z (long).
 * step               = spatial step length on the grid (double).
 * value_outside      = condition for accessing boundary points:
 *                      CGRID3D_DIRICHLET_BOUNDARY: Dirichlet boundary
 *                      or CGRID3D_NEUMANN_BOUNDARY: Neumann boundary
 *                      or CGRID3D_PERIODIC_BOUNDARY: Periodic boundary
 *                      or CGRID3D_VORTEX_X_BOUNDARY: Vortex long X
 *                      or CGRID3D_VORTEX_Y_BOUNDARY: Vortex long y
 *                      or CGRID3D_VORTEX_Z_BOUNDARY: Vortex long z
 *                      or user supplied function with pointer to grid and
 *                         grid index as parameters to provide boundary access.
 * outside_params_ptr = pointer for passing parameters for the given boundary
 *                      access function. Use 0 to with the predefined boundary
 *                      functions (void *).
 *
 * Return value: pointer to the allocated grid (cgrid3d *). Returns NULL on
 * error.
 *
 */

EXPORT cgrid3d *cgrid3d_alloc(long nx, long ny, long nz, double step, double complex (*value_outside)(const cgrid3d *grid, long i, long j, long k), void *outside_params_ptr) {

  cgrid3d *grid;
  
  grid = (cgrid3d *) malloc(sizeof(cgrid3d));
  if (!grid) {
    fprintf( stderr, "libgrid: Error in cgrid3d_alloc(). Could not allocate memory for 3d grid.\n");
    return 0;
  }
  
  if (!(grid->value = (double complex *) fftw_malloc(nx * ny * nz * sizeof(double complex)))) {
    fprintf(stderr, "libgrid: Error in cgrid3d_alloc(). Could not allocate memory for cgrid3d->value.\n");
    free(grid);
    return 0;
  }
  
  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->step = step;
  /* Set the origin of coordinates to its default value */
  grid->x0 = 0. ;
  grid->y0 = 0. ;
  grid->z0 = 0. ;
  /* Set the origin of momentum (i.e. frame of reference velocity) to its default value */
  grid->kx0 = 0. ;
  grid->ky0 = 0. ;
  grid->kz0 = 0. ;
  
  if (value_outside)
    grid->value_outside = value_outside;
  else
    grid->value_outside = cgrid3d_value_outside_constantdirichlet;
    
  if (outside_params_ptr)
    grid->outside_params_ptr = outside_params_ptr;
  else {
    grid->default_outside_params = 0.0;
    grid->outside_params_ptr = &grid->default_outside_params;
  }
  
  grid->plan = grid->iplan = grid->implan = grid->iimplan = NULL;
  
#if __DEBUG__
  allocated_grids++;
  fprintf(stderr, "libgrid(debug): %3d 3d complex grids allocated.\n", allocated_grids);
#endif

  cgrid3d_constant(grid, NAN);
  
  return grid;
}

/*
 * Set the origin of coordinates, meaning the coordinates of the grid will be:
 * 	x(i)  = (i - nx/2)* step - x0
 * 	y(j)  = (j - ny/2)* step - y0
 * 	z(k)  = (k - nz/2)* step - z0
 */
EXPORT void cgrid3d_set_origin( cgrid3d *grid , double x0, double y0, double z0){
	grid->x0 = x0 ;
	grid->y0 = y0 ;
	grid->z0 = z0 ;
}

/* Shift the origin */
EXPORT void cgrid3d_shift_origin( cgrid3d *grid , double x0, double y0, double z0){
	grid->x0 += x0 ;
	grid->y0 += y0 ;
	grid->z0 += z0 ;
}

/*
 * Set the origin of momentum space, or the velocity of the frame of reference. 
 *  kx0, ky0 and kz0 can be any real numbers but keep in mind that the grid
 *  will only contain the point k=0 if they are multiples of 
 *  kx0min = 2. * M_PI * HBAR / (NX * STEP * MASS) 
 *  ky0min = 2. * M_PI * HBAR / (NY * STEP * MASS) 
 *  kz0min = 2. * M_PI * HBAR / (NZ * STEP * MASS)
 *
 */
EXPORT void cgrid3d_set_momentum( cgrid3d *grid , double kx0, double ky0, double kz0){
	grid->kx0 = kx0 ;
	grid->ky0 = ky0 ;
	grid->kz0 = kz0 ;
}

/*
 * Free 3D grid.
 *
 * grid = pointer to 3D grid to be freed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_free(cgrid3d *grid) {

  if (grid) {
    if (grid->value) fftw_free(grid->value);
    cgrid3d_fftw_free(grid);
    free(grid);
  }
}

/* 
 * Write 3D grid on disk in binary format.
 *
 * grid = 3D grid to be written (cgrid3d *).
 * out  = file handle for the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_write(cgrid3d *grid, FILE *out) {

  fwrite(&grid->nx, sizeof(long), 1, out);
  fwrite(&grid->ny, sizeof(long), 1, out);
  fwrite(&grid->nz, sizeof(long), 1, out);
  fwrite(&grid->step, sizeof(double), 1, out);
  fwrite(grid->value, sizeof(double complex), grid->nx * grid->ny * grid->nz, out);
}

/*
 * Peek at a grid file to get dimensions.
 *
 * fp = File pointer for operation (FILE *).
 * nx = # of points along x (long *).
 * ny = # of points along y (long *).
 * nz = # of points along z (long *).
 * step = spatial step length (double *).
 *
 * No return value.
 *
 * Notes: - This works for both real and complex grids and hence it is called
 *          grid3d_read_peek().
 *        - This rewinds the fp so that Xgrid3d_read() can be called directly
 *          after this.
 *
 */

EXPORT void grid3d_read_peek(FILE *fp, long *nx, long *ny, long *nz, double *step) {

  fread(nx, sizeof(long), 1, fp);
  fread(ny, sizeof(long), 1, fp);
  fread(nz, sizeof(long), 1, fp);
  fread(step, sizeof(double), 1, fp);
  rewind(fp);
}

/* 
 * Read 3D grid from disk in binary format.
 *
 * grid = 3D grid to be read (cgrid3d *).
 * in   = file handle for reading the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_read(cgrid3d *grid, FILE *in) {

  long nx, ny, nz;
  double step;
  
  fread(&nx, sizeof(long), 1, in);
  fread(&ny, sizeof(long), 1, in);
  fread(&nz, sizeof(long), 1, in);
  fread(&step, sizeof(double), 1, in);
  
  if (nx != grid->nx || ny != grid->ny || nz != grid->nz || step != grid->step) {
    cgrid3d *tmp;

    fprintf(stderr, "libgrid: Grid in file has different size than grid in memory.\n");
    fprintf(stderr, "libgrid: Interpolating between grids.\n");
    if(!(tmp = cgrid3d_alloc(nx, ny, nz, step, grid->value_outside, NULL))) {
      fprintf(stderr, "libgrid: Error allocating grid in cgrid3d_read().\n");
      abort();
    }
    fread(tmp->value, sizeof(double complex), nx * ny * nz, in);
    cgrid3d_extrapolate(grid, tmp);
    cgrid3d_free(tmp);
    return;
  }
  
  fread(grid->value, sizeof(double complex), grid->nx * grid->ny * grid->nz, in);
}

/*
 * Copy 3D grid from one grid to another.
 *
 * copy = destination grid (cgrid3d *).
 * grid = source grid (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_copy(cgrid3d *copy, const cgrid3d *grid) {

  long i, nx = grid->nx, nyz = grid->ny * grid->nz, bytes = grid->ny * grid->nz * sizeof(double complex);
  double complex *gvalue = grid->value;
  double complex *cvalue = copy->value;
  
  copy->nx = grid->nx;
  copy->ny = grid->ny;
  copy->nz = grid->nz;
  copy->step = grid->step;
  
  copy->x0 = grid->x0;
  copy->y0 = grid->y0;
  copy->z0 = grid->z0;
  copy->kx0 = grid->kx0;
  copy->ky0 = grid->ky0;
  copy->kz0 = grid->kz0;
#pragma omp parallel for firstprivate(nx,nyz,bytes,gvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    memmove(&cvalue[i*nyz], &gvalue[i*nyz], bytes);
}

/*
 * Take complex conjugate of 3D grid.
 * 
 * conjugate = destination for complex conjugated grid (cgrid3d *).
 * grid      = source grid for the operation (cgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination may be the same grid.
 * 
 */

EXPORT void cgrid3d_conjugate(cgrid3d *conjugate, const cgrid3d *grid) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double complex *cvalue = conjugate->value;
  double complex *gvalue = grid->value;
  
  conjugate->nx = grid->nx;
  conjugate->ny = grid->ny;
  conjugate->nz = grid->nz;
  conjugate->step = grid->step;
  
#pragma omp parallel for firstprivate(nxy,nz,cvalue,gvalue) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = conj(gvalue[ijnz + k]);
  }
}

/*
 * Shift 3D grid by given amount spatially.
 *
 * shifted = destination grid for the operation (cgrid3d *).
 * grid    = source grid for the operation (cgrid3d *).
 * x       = shift spatially by this amount in x (double).
 * y       = shift spatially by this amount in y (double).
 * z       = shift spatially by this amount in z (double).
 *
 * No return value.
 *
 * NOTE: Source and destination may be the same grid.
 *
 */

EXPORT void cgrid3d_shift(cgrid3d *shifted, const cgrid3d *grid, double x, double y, double z) {

  sShiftParametersc3d params;

  /* shift by (x,y,z) i.e. current grid center to (x,y,z) */
  params.x = x;  params.y = y;  params.z = z;  params.grid = grid;
  cgrid3d_map(shifted, shift_cgrid3d, &params);
}

/* 
 * Zero 3D grid.
 *
 * grid = grid to be zeroed (cgrid3d *).
 *
 * No return value.
 * 
 */

EXPORT void cgrid3d_zero(cgrid3d *grid) { 

  cgrid3d_constant(grid, 0.0); 
}

/* 
 * Set 3D grid to a constant value.
 *
 * grid = grid to be set (cgrid3d *).
 * c    = value (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_constant(cgrid3d *grid, double complex c) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nz,value,c) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++) {
      value[ijnz + k] = c;
    }
  }
}

/*
 * Multiply a given grid by a function.
 *
 * grid = destination grid for the operation (cgrid3d *).
 * func = function providing the mapping (double complex (*)(void *, double, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and x,y,z are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_product_func(cgrid3d *grid, double complex (*func)(void *arg, double x, double y, double z), void *farg) {

  long i, j, k, ij, ijnz, nx = grid->nx , ny = grid->ny, nxy = nx * ny, nz = grid->nz;
  double x,y,z, step = grid->step;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nz,nxy,step,func,value,x0,y0,z0) private(i,j,ij,ijnz,k,x,y,z) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    x = (i - nx/2) * step - x0;
    y = (j - ny/2) * step - y0;    
    for(k = 0; k < nz; k++) {
      z = (k - nz/2) * step - z0 ;
      value[ijnz + k] *= func(farg, x, y, z);
    }
  }
}

/*
 * Map a given function onto 3D grid.
 *
 * grid = destination grid for the operation (cgrid3d *).
 * func = function providing the mapping (double complex (*)(void *, double, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and x,y,z are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_map(cgrid3d *grid, double complex (*func)(void *arg, double x, double y, double z), void *farg) {

  long i, j, k, ij, ijnz, nx = grid->nx , ny = grid->ny, nxy = nx * ny, nz = grid->nz;
  double x,y,z, step = grid->step;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nz,nxy,step,func,value,x0,y0,z0) private(i,j,ij,ijnz,k,x,y,z) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    x = (i - nx/2) * step - x0;
    y = (j - ny/2) * step - y0;
    for(k = 0; k < nz; k++) {
      z = (k - nz/2) * step - z0;
      value[ijnz + k] = func(farg, x, y, z);
    }
  }
}

/*
 * Map a given function onto 3D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid.
 * *
 * grid = destination grid for the operation (cgrid3d *).
 * func = function providing the mapping (double complex (*)(void *, double, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and (x, y, z) is the point (doubles) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 * ns   = number of intermediate points to be used in smoothing (int).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_smooth_map(cgrid3d *grid, double complex (*func)(void *arg, double x, double y, double z), void *farg, int ns) {

  long i, j, k, ij, ijnz, nx = grid->nx , ny = grid->ny, nxy = nx * ny, nz = grid->nz;
  double xc, yc, zc, step = grid->step;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nz,nxy,ns,step,func,value,x0,y0,z0) private(i,j,k,ijnz,xc,yc,zc) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    xc = (i - nx/2) * step - x0;
    yc = (j - ny/2) * step - y0;
    for(k = 0; k < nz; k++) {
      zc = (k - nz/2) * step - z0;
      value[ijnz + k] = linearly_weighted_integralc3d(func, farg, xc, yc, zc, step, ns);
    }
  }
}

/*
 * Map a given function onto 3D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid. Limits for intermediate steps and
 * tolerance can be given.
 *
 * grid   = destination grid for the operation (cgrid3d *).
 * func   = function providing the mapping (double complex (*)(void *, double, double, double)).
 *          The first argument (void *) is for external user specified data
 *          and x,y,z are the coordinates (double) where the function is evaluated.
 * farg   = pointer to user specified data (void *).
 * min_ns = minimum number of intermediate points to be used in smoothing (int).
 * max_ns = maximum number of intermediate points to be used in smoothing (int).
 * tol    = tolerance for weighing (double).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_adaptive_map(cgrid3d *grid, double complex (*func)(void *arg, double x, double y, double z), void *farg, int min_ns, int max_ns, double tol) {

  long i, j, k, ij, ijnz, nx = grid->nx , ny = grid->ny, nxy = nx * ny, nz = grid->nz, ns;
  double xc, yc, zc, step = grid->step;
  double tol2 = tol * tol;
  double complex sum, sump;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  double complex *value = grid->value;
  
  if (min_ns < 1) min_ns = 1;
  if (max_ns < min_ns) max_ns = min_ns;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nz,nxy,min_ns,max_ns,step,func,value,tol2,x0,y0,z0) private(i,j,k,ijnz,ns,xc,yc,zc,sum,sump) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    xc = (i - nx/2) * step - x0;
    yc = (j - ny/2) * step - y0;
    for(k = 0; k < nz; k++) {
      zc = (k - nz/2) * step - z0;
      sum  = func(farg, xc, yc, zc); sump = 0.0;
      for(ns = min_ns; ns <= max_ns; ns *= 2) {
        sum  = linearly_weighted_integralc3d(func, farg, xc, yc, zc, step, ns);
        sump = linearly_weighted_integralc3d(func, farg, xc, yc, zc, step, ns+1);
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
      
      value[ijnz + k] = 0.5 * (sum + sump);
    }
    
    /*fprintf(stderr, "\n");*/
  }
}

/*
 * Add two 3D grids ("gridc = grida + gridb").
 *
 * gridc = destination grid (cgrid3d *).
 * grida = 1st of the grids to be added (cgrid3d *).
 * gridb = 2nd of the grids to be added (cgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid3d_sum(cgrid3d *gridc, const cgrid3d *grida, const cgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,cvalue) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = avalue[ijnz + k] + bvalue[ijnz + k];
  }
}

/* 
 * Subtract two grids ("gridc = grida - gridb").
 *
 * gridc = destination grid (cgrid3d *).
 * grida = 1st source grid (cgrid3d *).
 * gridb = 2nd source grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: both source and destination may be the same.
 *
 */

EXPORT void cgrid3d_difference(cgrid3d *gridc, const cgrid3d *grida, const cgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,cvalue) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = avalue[ijnz + k] - bvalue[ijnz + k];
  }
}

/* 
 * Calculate product of two grids ("gridc = grida * gridb").
 *
 * gridc = destination grid (cgrid3d *).
 * grida = 1st source grid (cgrid3d *).
 * gridb = 2nd source grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid3d_product(cgrid3d *gridc, const cgrid3d *grida, const cgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,cvalue) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = avalue[ijnz + k] * bvalue[ijnz + k];
  }
}

/* 
 * Rise a grid to given power.
 *
 * gridb    = destination grid (cgrid3d *).
 * grida    = 1st source grid (cgrid3d *).
 * exponent = exponent to be used (double).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same.
 *       This routine uses pow() so that the exponent can be
 *       fractional but this is slow! Do not use this for integer
 *       exponents.
 *
 * TODO: Add complex exponentiation later.
 *
 */

EXPORT void cgrid3d_power(cgrid3d *gridb, const cgrid3d *grida, double exponent) {

  long ij, k, ijnz, nxy = gridb->nx * gridb->ny, nz = gridb->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,exponent) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      bvalue[ijnz + k] = pow(avalue[ijnz + k], exponent);
  }
}


/* 
 * Rise absolute value of a grid to given power.
 *
 * gridb    = destination grid (cgrid3d *).
 * grida    = 1st source grid (cgrid3d *).
 * exponent = exponent to be used (double).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same.
 *       This routine uses pow() so that the exponent can be
 *       fractional but this is slow! Do not use this for integer
 *       exponents.
 *
 * TODO: Add complex exponentiation later.
 *
 */

EXPORT void cgrid3d_abs_power(cgrid3d *gridb, const cgrid3d *grida, double exponent) {

  long ij, k, ijnz, nxy = gridb->nx * gridb->ny, nz = gridb->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,exponent) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      bvalue[ijnz + k] = pow(cabs(avalue[ijnz + k]), exponent);
  }
}

/*
 * Divide two grids ("gridc = grida / gridb").
 *
 * gridc = destination grid (cgrid3d *).
 * grida = 1st source grid (cgrid3d *).
 * gridb = 2nd source grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same. EPS added to avoid
 * possible NaNs.
 *
 */

EXPORT void cgrid3d_division(cgrid3d *gridc, const cgrid3d *grida, const cgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,cvalue) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = avalue[ijnz + k] / (bvalue[ijnz + k] + GRID_EPS);
  }
}

/* 
 * Conjugate product of two grids ("gridc = conj(grida) * gridb").
 *
 * gridc = destination grid (cgrid3d *).
 * grida = 1st source grid (cgrid3d *).
 * gridb = 2nd source grid (cgrid3d *).
 *
 * No return value.
 * 
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid3d_conjugate_product(cgrid3d *gridc, const cgrid3d *grida, const cgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,cvalue) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = conj(avalue[ijnz + k]) * bvalue[ijnz + k];
  }
}

/*
 * Add a constant to a 3D grid.
 *
 * grid = grid where the constant is added (cgrid3d *).
 * c    = constant to be added (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_add(cgrid3d *grid, double complex c) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nz,value,c) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      value[ijnz + k] += c;
  }
}

/*
 * Multiply grid by a constant.
 *
 * grid = grid to be multiplied (cgrid3d *).
 * c    = multiplier (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_multiply(cgrid3d *grid, double complex c) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nz,value,c) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      value[ijnz + k] *= c;
  }
}

/* 
 * Add and multiply: grid = (grid + ca) * cm.
 *
 * grid = grid to be operated (cgrid3d *).
 * ca   = constant to be added (double complex).
 * cm   = multiplier (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_add_and_multiply(cgrid3d *grid, double complex ca, double complex cm) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nz,value,ca,cm) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      value[ijnz + k] = (value[ijnz + k] + ca) * cm;
  }
}

/*
 * Multiply and add: grid = cm * grid + ca.
 *
 * grid = grid to be operated (cgrid3d *).
 * ca   = constant to be added (double complex).
 * cm   = multiplier (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_multiply_and_add(cgrid3d *grid, double complex cm, double complex ca) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nz,value,ca,cm) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      value[ijnz + k] = value[ijnz + k] * cm + ca;
  }
}

/* 
 * Add scaled grid (multiply/add): gridc = gridc + d * grida
 *
 * gridc = destination grid for the operation (cgrid3d *).
 * d     = multiplier for the operation (double complex).
 * grida = source grid for the operation (cgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid3d_add_scaled(cgrid3d *gridc, double complex d, const cgrid3d *grida) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double complex *avalue = grida->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(d,nxy,nz,avalue,cvalue) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] += d * avalue[ijnz + k];
  }
}

/*
 * Perform the following operation: gridc = gridc + d * grida * gridb.
 *
 * gridc = destination grid (cgrid3d *).
 * d     = constant multiplier (double complex).
 * grida = 1st source grid (cgrid3d *).
 * gridb = 2nd source grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid3d_add_scaled_product(cgrid3d *gridc, double complex d, const cgrid3d *grida, const cgrid3d *gridb) {
  
  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,cvalue,d) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] += d * avalue[ijnz + k] * bvalue[ijnz + k];
  }
}

/*
 * Operate on a grid by a given operator: gridc = O(grida).
 *
 * gridc    = destination grid (cgrid3d *).
 * grida    = source grid (cgrid3d *).
 * operator = operator (double complex (*)(double complex)).
 *            (i.e. a function mapping a given C-number to another)
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid3d_operate_one(cgrid3d *gridc, const cgrid3d *grida, double complex (*operator)(double complex a)) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double complex *avalue = grida->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,cvalue,operator) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = operator(avalue[ijnz + k]);
  }
}

/* 
 * Operate on two grids and place the result in third: gridc = O(grida, gridb).
 * where O is the operator.
 *
 * gridc    = destination grid (cgrid3d *).
 * grida    = 1s source grid (cgrid3d *).
 * gridb    = 2nd source grid (cgrid3d *).
 * operator = operator mapping grida and gridb (double complex (*)(double
 *            complex, double complex)).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid3d_operate_two(cgrid3d *gridc, const cgrid3d *grida, const cgrid3d *gridb, double complex (*operator)(double complex a, double complex b)) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,cvalue,operator) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = operator(avalue[ijnz + k], bvalue[ijnz + k]);
  }
}

/*
 * Operate on a grid by a given operator.
 *
 * grid     = grid to be operated (cgrid3d *).
 * operator = operator (void (*)(double complex *)).
 * 
 * No return value.
 *
 */

EXPORT void cgrid3d_transform_one(cgrid3d *grid, void (*operator)(double complex *a)) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nz,value,operator) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      operator(&value[ijnz + k]);
  }
}

/*
 * Operate on two separate grids by a given operator.
 *
 * grida    = grid to be operated (cgrid3d *).
 * gridb    = grid to be operated (cgrid3d *).
 * operator = operator (void (*)(double complex *)).
 * 
 * No return value.
 *
 */

EXPORT void cgrid3d_transform_two(cgrid3d *grida, cgrid3d *gridb, void (*operator)(double complex *a, double complex *b)) {

  long ij, k, ijnz, nxy = grida->nx * grida->ny, nz = grida->nz;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,operator) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      operator(&avalue[ijnz + k], &bvalue[ijnz + k]);
  }
}

/*
 * Integrate over a grid.
 *
 * grid = grid to be integrated (cgrid3d *).
 *
 * Returns the integral value (double complex).
 *
 */

EXPORT double complex cgrid3d_integral(const cgrid3d *grid) {

  long i, j, k, nx = grid->nx , ny =  grid->ny, nz = grid->nz;
  double complex sum = 0.0;
  
  // TODO: collapse(2) ? 
#pragma omp parallel for firstprivate(nx,ny,nz,grid) private(i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    for (j = 0; j <= ny; j++)
      for (k = 0; k <= nz; k++)
	sum += cgrid3d_value_at_index(grid, i, j, k);
  }
 
  return sum * grid->step * grid->step * grid->step;
}

/*
 * Integrate over a grid with limits.
 *
 * grid = grid to be integrated (cgrid3d *).
 * xl   = lower limit for x (double).
 * xu   = upper limit for x (double).
 * yl   = lower limit for y (double).
 * yu   = upper limit for y (double).
 * zl   = lower limit for z (double).
 * zu   = upper limit for z (double).
 *
 * Returns the integral value (double complex).
 *
 */

EXPORT double complex cgrid3d_integral_region(const cgrid3d *grid, double xl, double xu, double yl, double yu, double zl, double zu) {

  long iu, il, i, ju, jl, j, ku, kl, k;
  double complex sum;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  double step = grid->step;
   
  il = grid->nx/2 + (xl + x0)/ step ;
  iu = grid->nx/2 + (xu + x0)/ step ;
  jl = grid->ny/2 + (yl + y0)/ step ;
  ju = grid->ny/2 + (yu + y0)/ step ;
  kl = grid->nz/2 + (zl + z0)/ step ;
  ku = grid->nz/2 + (zu + z0)/ step ;
  
  sum = 0.0;
  // TODO: collapse(2) ?
#pragma omp parallel for firstprivate(il,iu,jl,ju,kl,ku,grid) private(i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = il; i <= iu; i++)
    for (j = jl; j <= ju; j++)
      for (k = kl; k <= ku; k++)
	sum += cgrid3d_value_at_index(grid, i, j, k);
 return sum * step * step * step; 
}
 
/* 
 * Integrate over the grid squared (int |grid|^2).
 *
 * grid = grid to be integrated (cgrid3d *).
 *
 * Returns the integral (double complex).
 *
 */

EXPORT double cgrid3d_integral_of_square(const cgrid3d *grid) {

  long i, j, k, nx = grid->nx , ny = grid->ny, nz = grid->nz;
  double sum = 0;
  
  // TODO: collapse(2) ?
#pragma omp parallel for firstprivate(nx,ny,nz,grid) private(i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    for (j = 0; j <= ny; j++)
      for (k = 0; k <= nz; k++)
	sum += sqnorm(cgrid3d_value_at_index(grid, i, j, k));
  }
 
  return sum * grid->step * grid->step * grid->step;
}

/*
 * Calculate overlap between two grids (int grida^*gridb).
 *
 * grida = 1st grid (complex conjugated) (cgrid3d *).
 * gridb = 2nd grid (no complex conjugation) (cgrid3d *).
 *
 * Returns the value of the overlap integral (double complex).
 *
 */

EXPORT double complex cgrid3d_integral_of_conjugate_product(const cgrid3d *grida, const cgrid3d *gridb) {

  long i, j, k, nx = grida->nx , ny = grida->ny, nz = grida->nz;
  double complex sum = 0.0;
  
  // TODO: collapse(2) ?
#pragma omp parallel for firstprivate(nx,ny,nz,grida,gridb) private(i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    for (j = 0; j <= ny; j++)
      for (k = 0; k <= nz; k++)
	sum += conj(cgrid3d_value_at_index(grida, i, j, k)) * cgrid3d_value_at_index(gridb, i, j, k);
  }   
  return sum * grida->step * grida->step * grida->step;
}

/*
 * Calculate the expectation value of a grid over a grid.
 * (int gridb^* grida gridb = int grida |gridb|^2).
 *
 * grida = grid giving the probability (|gridb|^2) (cgrid3d *).
 * gridb = grid to be averaged (cgrid3d *).
 *
 * Returns the average value (double complex).
 *
 */

EXPORT double complex cgrid3d_grid_expectation_value(const cgrid3d *grida, const cgrid3d *gridb) {

  long i, j, k, nx = grida->nx , ny = grida->ny , nz = grida->nz;
  double complex sum = 0.0;
  
  // TODO: collapse(2) ?
#pragma omp parallel for firstprivate(nx,ny,nz,grida,gridb) private(i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    for (j = 0; j <= ny; j++)
      for (k = 0; k <= nz; k++)
	sum += sqnorm(cgrid3d_value_at_index(grida, i, j, k)) * cgrid3d_value_at_index(gridb, i, j, k);
  }
 
  return sum * grida->step * grida->step * grida->step;
}
 
/*
 * Calculate the expectation value of a function over a grid.
 * (int grida^* func grida = int func |grida|^2).
 *
 * func  = function to be averaged (double complex (*)(void *, double complex, double, double, double)).
 *         The arguments are: optional arg, grida(x,y,z), x, y, z.
 * grida = grid giving the probability (|grida|^2) (cgrid3d *).
 *
 * Returns the average value (double complex).
 *
 */
 
EXPORT double complex cgrid3d_grid_expectation_value_func(void *arg, double complex (*func)(void *arg, double complex val, double x, double y, double z), const cgrid3d *grida) {
   
  long i, j, k, nx = grida->nx, ny = grida->ny, nz = grida->nz;
  double complex sum = 0.0, tmp;
  double x0 = grida->x0, y0 = grida->y0, z0 = grida->z0;
  double x, y, z, step = grida->step;
  
  // TODO: collapse(2) ? move x = ... inside the 2nd loop?
#pragma omp parallel for firstprivate(nx,ny,nz,grida,x0,y0,z0,step,func,arg) private(x,y,z,i,j,k,tmp) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    x = (i - nx/2) * step - x0;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step - y0;
      for (k = 0; k <= nz; k++) {
	z = (k - nz/2) * step - z0;
	tmp = cgrid3d_value_at_index(grida, i, j, k);
	sum += sqnorm(tmp) * func(arg, tmp, x, y, z);
      }
    }
  }
 
  return sum * step * step * step;
}


/* 
 * Integrate over the grid multiplied by weighting function (int grid w(x)).
 *
 * grid   = grid to be integrated over (cgrid3d *).
 * weight = function defining the weight (double complex (*)(double, double, double)). The arguments are (x,y,z) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double complex).
 *
 */

EXPORT double complex cgrid3d_weighted_integral(const cgrid3d *grid, double complex (*weight)(void *farg, double x, double y, double z), void *farg) {

  long i,j,k, nx = grid->nx , ny = grid->ny, nz = grid->nz;
  double x,y,z, step = grid->step;
  double complex sum = 0.0;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  
  // TODO: collapse(2) ? move x = ... inside the 2nd loop?
#pragma omp parallel for firstprivate(nx,ny,nz,grid,x0,y0,z0,step,weight,farg) private(x,y,z,i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    x = (i - nx/2) * step - x0;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step - y0;
      for (k = 0; k <= nz; k++) {
	z = (k - nz/2) * step - z0;
	sum += weight(farg, x, y, z) * cgrid3d_value_at_index(grid, i, j, k);
      }
    }
  }
 
  return sum * grid->step * grid->step * grid->step;
}

/* 
 * Integrate over square of the grid multiplied by weighting function (int grid^2 w(x)).
 *
 * grid   = grid to be integrated over (cgrid3d *).
 * weight = function defining the weight (double complex (*)(double, double, double)).
 *          The arguments are (x,y,z) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double complex).
 *
 */

EXPORT double cgrid3d_weighted_integral_of_square(const cgrid3d *grid, double (*weight)(void *farg, double x, double y, double z), void *farg) {

  long i,j,k, nx = grid->nx , ny = grid->ny, nz = grid->nz;
  double x,y,z, step = grid->step;
  double sum = 0;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  
  // TODO: collapse(2) ? move x = ... inside the 2nd loop?
#pragma omp parallel for firstprivate(nx,ny,nz,grid,x0,y0,z0,step,weight,farg) private(x,y,z,i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    x = (i - nx/2) * step - x0;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step - y0;
      for (k = 0; k <= nz; k++) {
	z = (k - nz/2) * step - z0;
	sum += weight(farg, x, y, z) * sqnorm(cgrid3d_value_at_index(grid, i, j, k));
      }
    }
  }
 
  return sum * grid->step * grid->step * grid->step;
}

/* 
 * Differentiate a grid with respect to x (central difference).
 *
 * grid     = grid to be differentiated (cgrid3d *).
 * gradient = differentiated grid output (cgrid3d *).
 * 
 * No return value.
 *
 */

EXPORT void cgrid3d_fd_gradient_x(const cgrid3d *grid, cgrid3d *gradient) {

  long i, j, k, ij, ijnz, ny = grid->ny, nz = grid->nz, nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double complex *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in cgrid3d_fd_gradient_x().\n");
    return;
  }

#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = inv_delta * (cgrid3d_value_at_index(grid, i+1, j, k) - cgrid3d_value_at_index(grid, i-1, j, k));
  }
}

/* 
 * Differentiate a grid with respect to y.
 *
 * grid     = grid to be differentiated (cgrid3d *).
 * gradient = differentiated grid output (cgrid3d *).
 * 
 * No return value.
 *
 */

EXPORT void cgrid3d_fd_gradient_y(const cgrid3d *grid, cgrid3d *gradient) {

  long i, j, k, ij, ijnz, ny = grid->ny, nz = grid->nz, nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double complex *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in cgrid3d_fd_gradient_y().\n");
    return;
  }

#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = inv_delta * (cgrid3d_value_at_index(grid, i, j+1, k) - cgrid3d_value_at_index(grid, i, j-1, k));
  }
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

EXPORT void cgrid3d_fd_gradient_z(const cgrid3d *grid, cgrid3d *gradient) {

  long i, j, k, ij, ijnz, ny = grid->ny, nz = grid->nz, nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double complex *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in cgrid3d_fd_gradient_z().\n");
    return;
  }

#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = inv_delta * (cgrid3d_value_at_index(grid, i, j, k+1) - cgrid3d_value_at_index(grid, i, j, k-1));
  }
}

/*
 * Calculate gradient of a grid.
 *
 * grid       = grid to be differentiated twice (cgrid3d *).
 * gradient_x = x output grid for the operation (cgrid3d *).
 * gradient_y = y output grid for the operation (cgrid3d *).
 * gradient_z = z output grid for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_fd_gradient(const cgrid3d *grid, cgrid3d *gradient_x, cgrid3d *gradient_y, cgrid3d *gradient_z) {

  cgrid3d_fd_gradient_x(grid, gradient_x);
  cgrid3d_fd_gradient_y(grid, gradient_y);
  cgrid3d_fd_gradient_z(grid, gradient_z);
}

/*
 * Calculate laplacian of the grid.
 *
 * grid    = source grid (cgrid3d *).
 * laplace = output grid for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_fd_laplace(const cgrid3d *grid, cgrid3d *laplace) {

  long i, j, k, ij, ijnz;
  long ny = grid->ny, nz = grid->nz;
  long nxy = grid->nx * grid->ny;
  double inv_delta2 = 1.0 / (grid->step * grid->step);
  double complex *lvalue = laplace->value;
  
#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta2,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = inv_delta2 * (-6.0 * cgrid3d_value_at_index(grid, i, j, k) + cgrid3d_value_at_index(grid, i, j, k+1)
				       + cgrid3d_value_at_index(grid, i, j, k-1) + cgrid3d_value_at_index(grid, i, j+1, k) 
				       + cgrid3d_value_at_index(grid, i, j-1, k) + cgrid3d_value_at_index(grid,i+1,j,k) 
				       + cgrid3d_value_at_index(grid,i-1,j,k));
  }
}

/*
 * Calculate dot product of the gradient of the grid.
 *
 * grid          = source grid for the operation (cgrid3d *).
 * grad_dot_grad = destination grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: grid and grad_dot_grad may not be the same grid.
 *
 */

EXPORT void cgrid3d_fd_gradient_dot_gradient(const cgrid3d *grid, cgrid3d *grad_dot_grad) {

  long i, j, k, ij, ijnz;
  long ny = grid->ny, nz = grid->nz;
  long nxy = grid->nx * grid->ny;
  double inv_2delta2 = 1.0 / (2.0*grid->step * 2.0*grid->step);
  double complex *gvalue = grad_dot_grad->value;
  
/*  grad f(x,y,z) dot grad f(x,y,z) = [ |f(+,0,0) - f(-,0,0)|^2 + |f(0,+,0) - f(0,-,0)|^2 + |f(0,0,+) - f(0,0,-)|^2 ] / (2h)^2 */
#pragma omp parallel for firstprivate(ny,nz,nxy,gvalue,inv_2delta2,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      gvalue[ijnz + k] = inv_2delta2 * 
	(sqnorm(cgrid3d_value_at_index(grid, i, j, k+1) - cgrid3d_value_at_index(grid, i, j, k-1)) + sqnorm(cgrid3d_value_at_index(grid, i, j+1, k) - cgrid3d_value_at_index(grid, i, j-1, k))
	 + sqnorm(cgrid3d_value_at_index(grid, i+1, j, k) - cgrid3d_value_at_index(grid, i-1, j, k)));
  }
}

/*
 * Print the grid with both real and imaginary parts into file (ASCII format).
 *
 * grid = grid to be printed out (cgrid3d *).
 * out  = output file pointer (FILE *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_print(const cgrid3d *grid, FILE *out) {

  long i, j, k;

  for(i = 0; i < grid->nx; i++) {
    for(j = 0; j < grid->ny; j++) {
      for(k = 0; k < grid->nz; k++) {
        fprintf(out, "%16.8le %16.8le   ", 
		creal(cgrid3d_value_at_index(grid, i, j, k)),
		cimag(cgrid3d_value_at_index(grid, i, j, k)));
	  }
      fprintf(out, "\n");
    }
    fprintf(out, "\n");
  }
}

/*
 * Perform Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (cgrid3d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       Also no normalization is performed.
 *
 */

EXPORT void cgrid3d_fft(cgrid3d *grid) {

  if (!grid->plan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in cgrid3d_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    cgrid3d_fftw_alloc(grid);
  }
  cgrid3d_fftw(grid);
}

/*
 * Perform inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid3d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       No normalization.
 *
 */

EXPORT void cgrid3d_inverse_fft(cgrid3d *grid) {

  if (!grid->iplan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in cgrid3d_inverse_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    cgrid3d_fftw_alloc(grid);
  }
  cgrid3d_fftw_inv(grid);
}

/*
 * Perform scaled inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid3d *).
 * c    = scaling factor (i.e. the output is multiplied by this constant) (double complex).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void cgrid3d_scaled_inverse_fft(cgrid3d *grid, double complex c) {

  cgrid3d_inverse_fft(grid);
  cgrid3d_multiply(grid, c);  
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by FFT norm.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid3d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void cgrid3d_inverse_fft_norm(cgrid3d *grid) {

  cgrid3d_scaled_inverse_fft(grid, grid->fft_norm);
}

/*
 * Convolue FFT transformed grids. To apply this on two grids (grida and gridb)
 * and place the result in gridc:
 * cgrid3d_fft(grida);
 * cgrid3d_fft(gridb);
 * cgrid3d_convolue(gridc, grida, gridb);
 * cgrid3d_inverse_fft(gridc);
 * gridc now contains the convolution of grida and gridb.
 *
 * grida = 1st grid to be convoluted (cgrid3d *).
 * gridb = 2nd grid to be convoluted (cgrid3d *).
 * gridc = output (cgrid3d *).
 *
 * No return value.
 *
 * Note: the input/output grids may be the same.
 *
 */

EXPORT void cgrid3d_fft_convolute(cgrid3d *gridc, const cgrid3d *grida, const cgrid3d *gridb) {

  long i,j,k, ij, ijnz, nx, ny, nz, nxy;
  double step, norm;
  double complex *cvalue, *bvalue;

  /* int f(r) g(r-r') d^3r' = iF[ F[f] F[g] ] = (step / N)^3 iFFT[ FFT[f] FFT[g] ] */
  nx = gridc->nx;
  ny = gridc->ny;
  nz = gridc->nz;
  nxy = nx * ny;
  step = gridc->step;
  
  if (gridc != grida)
    cgrid3d_copy(gridc, grida);
  
  cvalue = gridc->value;
  bvalue = gridb->value;
  
  norm = step * step * step * grida->fft_norm;         /* David: fft_norm */
  
#pragma omp parallel for firstprivate(nx,ny,nz,nxy,cvalue,bvalue,norm) private(i,j,ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++) {
      /* if odd */
      if ((i + j + k) & 1)
        cvalue[ijnz + k] = -norm * cvalue[ijnz + k] * bvalue[ijnz + k];
      else
        cvalue[ijnz + k] = norm * cvalue[ijnz + k] * bvalue[ijnz + k];
    }
  }
}

/*
 * Differentiate grid in the Fourier space along x.
 *
 * grid       = grid to be differentiated (in Fourier space) (cgrid3d *).
 * gradient_x = output grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: input and output grids may be the same.
 *
 */

EXPORT void cgrid3d_fft_gradient_x(const cgrid3d *grid, cgrid3d *gradient_x) {

  long i, k, ij, ijnz, nx, ny, nz, nxy;
  double kx0 = grid->kx0 ;
  double kx, step, norm;
  double complex *gxvalue = gradient_x->value;
  
  /* f'(x) = iF[ i kx F[f(x)] ] */  
  nx = grid->nx;
  ny = grid->ny;
  nz = grid->nz;
  nxy = nx * ny;
  step = grid->step;
  
  /* David: fft_norm */
  norm = grid->fft_norm;
  
  if (gradient_x != grid)
    cgrid3d_copy(gradient_x, grid);

  if( grid -> value_outside == CGRID3D_NEUMANN_BOUNDARY  ||
      grid -> value_outside == CGRID3D_VORTEX_X_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Y_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Z_BOUNDARY )
    {
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,gxvalue,kx0) private(i,ij,ijnz,k,kx) default(none) schedule(runtime)
      for(ij = 0; ij < nxy; ij++) {
        i = ij / ny;
        ijnz = ij * nz;
	kx = M_PI * i / (nx * step) - kx0;
        for(k = 0; k < nz; k++)	  
          gxvalue[ijnz + k] *= (kx * norm) * I;
      }
    }
    else{
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,gxvalue,kx0) private(i,ij,ijnz,k,kx) default(none) schedule(runtime)
      for(ij = 0; ij < nxy; ij++) {
        i = ij / ny;
        ijnz = ij * nz;
        
        /* 
         * k = 2 pi n / L 
         * if k < n/2, k = k
         * else k = -k
         */
        if (i < nx / 2)
          kx = 2.0 * M_PI * i / (nx * step) - kx0;
        else 
          kx = 2.0 * M_PI * (i - nx) / (nx * step) - kx0;
        
        for(k = 0; k < nz; k++)	  
          gxvalue[ijnz + k] *= (kx * norm) * I;
      }
    }
}

/*
 * Differentiate grid in the Fourier space along y.
 *
 * grid       = grid to be differentiated (in Fourier space) (cgrid3d *).
 * gradient_y = output grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: input and output grids may be the same.
 *
 */

EXPORT void cgrid3d_fft_gradient_y(const cgrid3d *grid, cgrid3d *gradient_y) {

  long j, k, ij, ijnz, nx, ny, nz, nxy;
  double ky, step, norm;
  double ky0 = grid->ky0 ;
  double complex *gyvalue = gradient_y->value;
  
  /* f'(y) = iF[ i ky F[f(y)] ] */  
  nx = grid->nx;
  ny = grid->ny;
  nz = grid->nz;
  nxy = nx * ny;
  step = grid->step;
  
  /* David: fft_norm */
  norm = grid->fft_norm;
  
  if (gradient_y != grid)
    cgrid3d_copy(gradient_y, grid);
  if( grid -> value_outside == CGRID3D_NEUMANN_BOUNDARY  ||
      grid -> value_outside == CGRID3D_VORTEX_X_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Y_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Z_BOUNDARY )
  {
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,gyvalue,ky0) private(j,ij,ijnz,k,ky) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      j = ij % ny;
      ijnz = ij * nz;
        ky = M_PI * j / (ny * step) - ky0 ;
      for(k = 0; k < nz; k++)	  
        gyvalue[ijnz + k] *= ky * norm * I;
    }
  }
  else{
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,gyvalue,ky0) private(j,ij,ijnz,k,ky) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      j = ij % ny;
      ijnz = ij * nz;
      
      /* 
       * k = 2 pi n / L 
       * if k < n/2, k = k
       * else k = -k
       */
      if (j < ny / 2)
        ky = 2.0 * M_PI * j / (ny * step) - ky0 ;
      else 
        ky = 2.0 * M_PI * (j - ny) / (ny * step) - ky0 ;
      
      for(k = 0; k < nz; k++)	  
        gyvalue[ijnz + k] *= ky * norm * I;
    }
  }
}

/*
 * Differentiate grid in the Fourier space along z.
 *
 * grid       = grid to be differentiated (in Fourier space) (cgrid3d *).
 * gradient_z = output grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: input and output grids may be the same.
 *
 */

EXPORT void cgrid3d_fft_gradient_z(const cgrid3d *grid, cgrid3d *gradient_z) {

  long k, ij, ijnz, nx, ny, nz, nxy;
  double kz, lz, step, norm;
  double kz0 = grid->kz0 ;
  double complex *gzvalue = gradient_z->value;
  
  /* f'(z) = iF[ i kz F[f(z)] ] */
  nx = grid->nx;
  ny = grid->ny;
  nz = grid->nz;
  nxy = nx * ny;
  step = grid->step;
  
  /* David: fft_norm */
  norm = grid->fft_norm;
  
  if (gradient_z != grid)
    cgrid3d_copy(gradient_z, grid);
  if( grid -> value_outside == CGRID3D_NEUMANN_BOUNDARY  ||
      grid -> value_outside == CGRID3D_VORTEX_X_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Y_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Z_BOUNDARY )
    {
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,gzvalue,kz0) private(ij,ijnz,k,kz,lz) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      ijnz = ij * nz;
      
      lz = nz * step;
      for(k = 0; k < nz; k++) {
        kz = M_PI * k / lz - kz0;
        gzvalue[ijnz + k] *= kz * norm * I;
      }
    }
  }
  else{
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,gzvalue,kz0) private(ij,ijnz,k,kz,lz) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      ijnz = ij * nz;
      
      /* 
       * k = 2 pi n / L 
       * if k < n/2, k = k
       * else k = -k
       */
      
      lz = nz * step;
      for(k = 0; k < nz; k++) {
        if (k < nz / 2)
          kz = 2.0 * M_PI * k / lz - kz0;
        else 
          kz = 2.0 * M_PI * (k - nz) / lz - kz0;
        
        gzvalue[ijnz + k] *= kz * norm * I;
      }
    }    
  }

}

/* 
 * Calculate gradient of a grid (in Fourier space).
 *
 * grid       = grid to be differentiated (cgrid3d *).
 * gradient_x = x output grid (cgrid3d *).
 * gradient_y = y output grid (cgrid3d *).
 * gradient_z = z output grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: input/output grids may be the same.
 *
 */

EXPORT void cgrid3d_fft_gradient(const cgrid3d *grid, cgrid3d *gradient_x, cgrid3d *gradient_y, cgrid3d *gradient_z) {

  cgrid3d_fft_gradient_x(grid, gradient_x);
  cgrid3d_fft_gradient_y(grid, gradient_y);
  cgrid3d_fft_gradient_z(grid, gradient_z);
}

/* 
 * Calculate second derivative of a grid (in Fourier space).
 *
 * grid    = grid to be differentiated (cgrid3d *).
 * laplace = output grid (cgrid3d *).
 *
 * No return value.
 *
 * Note: input/output grids may be the same.
 *
 */

EXPORT void cgrid3d_fft_laplace(const cgrid3d *grid, cgrid3d *laplace)  {

  long i,j,k, ij, ijnz, nx,ny,nz, nxy;
  double kx0 = grid->kx0 , ky0 = grid->ky0 , kz0 = grid->kz0 ;
  double kx, ky, kz, lz, step, norm;
  double complex *lvalue = laplace->value;
  
  /* f''(x) = iF[ -k^2 F[f(x)] ] */
  nx = grid->nx;
  ny = grid->ny;
  nz = grid->nz;
  nxy = nx * ny;
  step = grid->step;
  
  norm = grid->fft_norm;
  
  if (grid != laplace)
    cgrid3d_copy(laplace, grid);
  
  if( grid -> value_outside == CGRID3D_NEUMANN_BOUNDARY  ||
      grid -> value_outside == CGRID3D_VORTEX_X_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Y_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Z_BOUNDARY )	
  {
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,lvalue,kx0,ky0,kz0) private(i,j,ij,ijnz,k,kx,ky,kz,lz) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      i = ij / ny;
      j = ij % ny;
      ijnz = ij * nz;
      
      kx = M_PI * i / (nx * step) - kx0;
      ky = M_PI * j / (ny * step) - ky0;
      lz = nz * step;
      for(k = 0; k < nz; k++) {
        kz = M_PI * k / lz - kz0;
        lvalue[ijnz + k] *= (-kx*kx -ky*ky -kz*kz) * norm * I;
      }
    }    
  }
  else{
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,lvalue,kx0,ky0,kz0) private(i,j,ij,ijnz,k,kx,ky,kz,lz) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      i = ij / ny;
      j = ij % ny;
      ijnz = ij * nz;
      
      /* 
       * k = 2 pi n / L 
       * if k < n/2, k = k
       * else k = -k
       */
      if (i < nx / 2)
        kx = 2.0 * M_PI * i / (nx * step) - kx0;
      else 
        kx = 2.0 * M_PI * (i - nx) / (nx * step) -kx0;
      
      if (j < ny / 2)
        ky = 2.0 * M_PI * j / (ny * step) - ky0;
      else 
        ky = 2.0 * M_PI * (j - ny) / (ny * step) - ky0;
      
      lz = nz * step;
      for(k = 0; k < nz; k++) {
        if (k < nz / 2)
          kz = 2.0 * M_PI * k / lz - kz0;
        else 
          kz = 2.0 * M_PI * (k - nz) / lz - kz0;
        
        lvalue[ijnz + k] *= (-kx*kx -ky*ky -kz*kz) * norm * I;
      }
    }
  }

}

/*
 * Calculate expectation value of laplace operator in the Fourier space (int grid^* grid'').
 *
 * grid    = source grid for the operation (in Fourier space) (cgrid3d *).
 * laplace = laplacian of the grid (input) (cgrid3d *).
 *
 * Returns the expectation value (double).
 *
 */

EXPORT double cgrid3d_fft_laplace_expectation_value(const cgrid3d *grid, cgrid3d *laplace)  {

  long i,j,k, ij, ijnz, nx,ny,nz, nxy;
  double kx, ky, kz, lz, step, norm, sum = 0, ssum;
  double kx0 = grid->kx0 , ky0 = grid->ky0 , kz0 = grid->kz0 ;
  double complex *lvalue = laplace->value;
  double aux ;
  /* int f*(x) f''(x) dx */
  nx = grid->nx;
  ny = grid->ny;
  nz = grid->nz;
  nxy = nx * ny;
  step = grid->step;
  
  /* int (delta FFT[f(x)] )^2 dk => delta^2 / N delta */
  /* David: Normalization */
  norm = grid->fft_norm;
  norm = step * step * step * grid->fft_norm;
  
  if (grid != laplace)
    cgrid3d_copy(laplace, grid);
  if( grid -> value_outside == CGRID3D_NEUMANN_BOUNDARY  ||
      grid -> value_outside == CGRID3D_VORTEX_X_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Y_BOUNDARY ||
      grid -> value_outside == CGRID3D_VORTEX_Z_BOUNDARY )
    {
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,lvalue,kx0,ky0,kz0) private(i,j,ij,ijnz,k,kx,ky,kz,lz, ssum, aux) reduction(+:sum) default(none) schedule(runtime)
   for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    ijnz = ij * nz;
   
    kx = M_PI * i / (nx * step) -kx0;
    ky = M_PI * j / (ny * step) -ky0;
    
    ssum = 0.0;
    lz = nz * step;
    
    for(k = 0; k < nz; k++) {
      kz = M_PI * k / lz - kz0;
      /* Manual fixing of boundaries: the symmetry points (i=0 or i=nx-1 etc) have 1/2 the weigth in the integral */
      aux = (- kx*kx - ky*ky - kz*kz) * sqnorm(lvalue[ijnz + k]) ;
      if(i==0 || i==nx-1) aux *= 0.5 ;
      if(j==0 || j==ny-1) aux *= 0.5 ;
      if(k==0 || k==nz-1) aux *= 0.5 ;
      ssum += aux ; 
      //ssum += (- kx*kx - ky*ky - kz*kz) * sqnorm(lvalue[ijnz + k]);
    }
    
    sum += ssum;
   }
    }else{  
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,lvalue,kx0,ky0,kz0) private(i,j,ij,ijnz,k,kx,ky,kz,lz, ssum) reduction(+:sum) default(none) schedule(runtime)
   for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    ijnz = ij * nz;
    
    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (i < nx / 2)
      kx = 2.0 * M_PI * i / (nx * step) - kx0;
    else 
      kx = 2.0 * M_PI * (i - nx) / (nx * step) - kx0;
    
    if (j < ny / 2)
      ky = 2.0 * M_PI * j / (ny * step) - ky0;
    else 
      ky = 2.0 * M_PI * (j - ny) / (ny * step) - ky0;
    
    ssum = 0.0;
    lz = nz * step;
    
    for(k = 0; k < nz; k++) {
      if (k < nz / 2)
        kz = 2.0 * M_PI * k / lz - kz0;
      else 
        kz = 2.0 * M_PI * (k - nz) / lz - kz0;
      
      ssum += (- kx*kx - ky*ky - kz*kz) * sqnorm(lvalue[ijnz + k]);
    }
    
    sum += ssum;
   }
  }

  return sum * norm;
}

/* Boundary condition routines */

EXPORT double complex cgrid3d_value_outside_constantdirichlet(const cgrid3d *grid, long i, long j, long k) {

  return *((double complex *) grid->outside_params_ptr);
}

/*
 * The symmetry point are i=0 and i=nx-1 for consistency with the FFT plan FFTW_REDFT00. 
 * If one wants to use REDFT01 the symmetry points are i=-0.5 and i=nx-0.5 
 * TODO: do we want to have nd - i - 1 (FFT compatibility)?
 */
EXPORT double complex cgrid3d_value_outside_neumann(const cgrid3d *grid, long i, long j, long k) {

  long nx = grid->nx, ny = grid->ny, nz = grid->nz, nd;

  nd = nx * 2;
  if (i < 0) i  = -i ;
  if (i >= nd) i %= nd;
  if (i >= nx) i  = nd - i -1;
  
  nd = ny * 2;
  if (j < 0) j  = -j ;
  if (j >= nd) j %= nd;
  if (j >= ny) j  = nd - j -1;
  
  nd = nz * 2;
  if (k < 0) k  = -k ;
  if (k >= nd) k %= nd;
  if (k >= nz) k  = nd - k -1;

  return grid->value[(i*ny + j)*nz + k];
}

EXPORT double complex cgrid3d_value_outside_periodic(const cgrid3d *grid, long i, long j, long k) {

  long nx = grid->nx, ny = grid->ny, nz = grid->nz;
  
  i %= nx;
  if (i < 0) i = nx + i;
  j %= ny;
  if (j < 0) j = ny + j;
  k %= nz;
  if (k < 0) k = nz + k;
  
  return grid->value[(i*ny + j)*nz + k];  
}

/* TODO Fix these */
EXPORT double complex cgrid3d_value_outside_vortex_x(const cgrid3d *grid, long i, long j, long k) {

  long nx = grid->nx, ny = grid->ny, nz = grid->nz, nd;

  /* In vortex direction, map i to {0,nx-1}. This does not introduce any phase */
  nd = nx * 2;
  if (i < 0) i  = -i-1;
  if (i >= nd) i %= nd;
  if (i >= nx) i  = nd -1 -i;
  
  /* First, map j,k to {0,2*nx-1} range
   * (this does not introduce any phase)
   */
  nd = ny * 2;
  if (j < 0) j = j%nd + nd;
  if (j >= nd) j %= nd;

  nd = nz * 2;
  if (k < 0 ) k = k%nd + nd;
  if (k >= nd) k %= nd;

  /* Then, if j has to be mapped to {0,nx-1} return -phi*
   *       if k has to be mapped to {0,ny-1} return +phi*
   *       if both have to be mapped return -phi
   */
  if(j>=ny){
	  j  = 2*ny -1 -j;
	  if(k>=nz){
		  k = nd -1 -k;
		  return -grid->value[(i*ny + j)*nz + k];
	  }
	  else
		  return -conj(grid->value[(i*ny + j)*nz + k]);
  }
  else{
	  if(k>=nz){
		  k = nd -1 -k;
		  return conj(grid->value[(i*ny + j)*nz + k]);
	  }
	  else
		  return grid->value[(i*ny + j)*nz + k];
  }
}

/* TODO Fix these */
EXPORT double complex cgrid3d_value_outside_vortex_y(const cgrid3d *grid, long i, long j, long k) {

  long nx = grid->nx, ny = grid->ny, nz = grid->nz, nd;

  /* In vortex direction, map i to {0,nx-1}. This does not introduce any phase */
  nd = ny * 2;
  if (j < 0) j  = -j-1;
  if (j >= nd) j %= nd;
  if (j >= ny) j  = nd -1 -j;

  /* First, map i,j to {0,2*nx-1} range
   * (this does not introduce any phase)
   */
  nd = nz * 2;
  if (k < 0) k = k%nd + nd;
  if (k >= nd) k %= nd;

  nd = nx * 2;
  if (i < 0 ) i = i%nd + nd;
  if (i >= nd) i %= nd;

  /* Then, if i has to be mapped to {0,nx-1} return -phi*
   *       if j has to be mapped to {0,ny-1} return +phi*
   *       if both have to be mapped return -phi
   */
  if(k>=nz){
	  k  = 2*nz -1 -k;
	  if(i>=nx){
		  i  = nd -1 -i;
		  return -grid->value[(i*ny + j)*nz + k];
	  }
	  else
		  return -conj(grid->value[(i*ny + j)*nz + k]);
  }
  else{
	  if(i>=nx){
		  i  = nd -1 -i;
		  return conj(grid->value[(i*ny + j)*nz + k]);
	  }
	  else
		  return grid->value[(i*ny + j)*nz + k];
  }
}

/* TODO Fix these */
EXPORT double complex cgrid3d_value_outside_vortex_z(const cgrid3d *grid, long i, long j, long k) {

  long nx = grid->nx, ny = grid->ny, nz = grid->nz, nd;

  /* In vortex direction, map k to {0,nx-1}. This does not introduce any phase */
  nd = nz * 2;
  if (k < 0) k  = -k-1;
  if (k >= nd) k %= nd;
  if (k >= nz) k  = nd -1 -k;
  
  /* First, map i,j to {0,2*nx-1} range
   * (this does not introduce any phase)
   */
  nd = nx * 2;
  if (i < 0) i = i%nd + nd;
  if (i >= nd) i %= nd;

  nd = ny * 2;
  if (j < 0 ) j = j%nd + nd;
  if (j >= nd) j %= nd;

  /* Then, if i has to be mapped to {0,nx-1} return -phi*
   *       if j has to be mapped to {0,ny-1} return +phi*
   *       if both have to be mapped return -phi
   */
  if(i>=nx){
	  i  = 2*nx -1 -i;
	  if(j>=ny){
		  j  = nd -1 -j;
		  return -grid->value[(i*ny + j)*nz + k];
	  }
	  else
		  return -conj(grid->value[(i*ny + j)*nz + k]);
  }
  else{
	  if(j>=ny){
		  j  = nd -1 -j;
		  return conj(grid->value[(i*ny + j)*nz + k]);
	  }
	  else
		  return grid->value[(i*ny + j)*nz + k];
  }
}

/* End boundary condition routines */

/*
 * Access grid point at given index.
 *
 * grid = grid to be accessed (cgrid3d *).
 * i    = index along x (long).
 * j    = index along y (long).
 * k    = index along z (long).
 *
 * Returns grid value at index (i, j, k).
 *
 */

EXPORT inline double complex cgrid3d_value_at_index(const cgrid3d *grid, long i, long j, long k) {

  if (i < 0 || j < 0 || k < 0 || i >= grid->nx || j >= grid->ny || k >= grid->nz)
    return grid->value_outside(grid, i, j, k);
  return grid->value[(i*grid->ny + j)*grid->nz + k];
}

/*
 * Access grid point at given (x,y,z) point (with linear interpolation).
 *
 * grid = grid to be accessed (cgrid3d *).
 * x    = x value (double).
 * y    = y value (double).
 * z    = z value (double).
 *
 * Returns grid value at (x,y,z).
 *
 */

EXPORT inline double complex cgrid3d_value(const cgrid3d *grid, double x, double y, double z) {

  double complex f000, f100, f010, f001, f110, f101, f011, f111;
  double omx, omy, omz, step = grid->step;
  long i, j, k;
  
  /* i to index and 0 <= x < 1 */
  x = ( x + grid->x0) / step ;
  i = x ;
  if (x < 0) i--;
  x -= i;
  i += grid->nx/2;
  
  /* j to index and 0 <= y < 1 */
  y = ( y + grid->y0) / step ;
  j = y ;
  if (y < 0) j--;
  y -= j;
  j += grid->ny/2;
  
  /* k to index and 0 <= z < 1 */
  z = ( z + grid->z0) / step ;
  k = z ;
  if (z < 0) k--;
  z -= k;
  k += grid->nz/2;
  
  /* linear extrapolation 
   *
   * f(x,y) = (1-x) (1-y) (1-z) f(0,0,0) + x (1-y) (1-z) f(1,0,0) + (1-x) y (1-z) f(0,1,0) + (1-x) (1-y) z f(0,0,1) 
   *          + x     y   (1-z) f(1,1,0) + x (1-y)   z   f(1,0,1) + (1-x) y   z   f(0,1,1) +   x     y   z f(1,1,1)
   */ 
  f000 = cgrid3d_value_at_index(grid, i, j, k);
  f100 = cgrid3d_value_at_index(grid, i+1, j, k);
  f010 = cgrid3d_value_at_index(grid, i, j+1, k);
  f001 = cgrid3d_value_at_index(grid, i, j, k+1);
  f110 = cgrid3d_value_at_index(grid, i+1, j+1, k);
  f101 = cgrid3d_value_at_index(grid, i+1, j, k+1);
  f011 = cgrid3d_value_at_index(grid, i, j+1, k+1);
  f111 = cgrid3d_value_at_index(grid, i+1, j+1, k+1);
  
  omx = 1.0 - x;
  omy = 1.0 - y;
  omz = 1.0 - z;

  return omx * (omy * (omz * f000 + z * f001) + y * (omz * f010 + z * f011))
    + x * (omy * (omz * f100 + z * f101) + y * (omz * f110 + z * f111));
}

/*
 * Extrapolate between two different grid sizes.
 *
 * dest = Destination grid (cgrid3d *; output).
 * src  = Source grid (cgrid3d *; input).
 *
 */

EXPORT void cgrid3d_extrapolate(cgrid3d *dest, cgrid3d *src) {

  long i, j, k, nx = dest->nx, ny = dest->ny, nz = dest->nz;
  double x0 = dest->x0, y0 = dest->y0, z0 = dest->z0;
  double step = dest->step, x, y, z;

  for (i = 0; i < nx; i++) {
    x = (i - nx/2) * step - x0 ;
    for (j = 0; j < ny; j++) {
      y = (j - ny/2) * step - y0 ;
      for (k = 0; k < nz; k++) {
	z = (k - nz/2) * step - z0 ;
	dest->value[i * ny * nz + j * nz + k] = cgrid3d_value(src, x, y, z);
      }
    }
  }
}

/*
 * Gives the value of the grid rotated by an angle theta.
 * The grid and the angle is specified through the rotation structure,
 * passed as a void.
 * The angle is giving as sin and cos for efficiency.
 *
 * arg is a rotation with the grid (either real or complex) and sin and cos
 * of the rotation angle.
 */

double complex cgrid3d_value_rotate_z(void *arg, double x, double y, double z) {

  /* Unpack the values in arg */ 
  cgrid3d *grid = ((rotation *) arg)->cgrid;
  double sth = ((rotation *) arg)->sinth, cth = ((rotation *) arg)->costh, xp, yp;

  xp = -y * sth + x * cth;
  yp =  y * cth + x * sth;

  return cgrid3d_value(grid, xp, yp, z);
}

/*
 * Rotate a grid by a given angle around the z-axis.
 *  cgrid3d *in : pointer with original grid.
 *  cgrid3d *out : pointer with rotated grid.
 *  double th: angle (radians) of the rotation.
 *
 *  The grid in and out CANNOT be the same.
 */

EXPORT void cgrid3d_rotate_z(cgrid3d *out, cgrid3d *in, double th) {

  rotation *r;

  if (in == out){
	  fprintf(stderr,"libgrid: in and out grids in rgrid3d_rotate_z must be different\n") ;
	  abort() ;
  }
  r = malloc(sizeof(rotation));
  r->cgrid = in;
  r->sinth = sin(-th);  // same direction of rotation as -wLz
  r->costh = cos(th);
  
  cgrid3d_map(out, cgrid3d_value_rotate_z, (void *) r);
  free(r);
}

/*
 * Clear real part of complex grid.
 *
 */

EXPORT void cgrid3d_zero_re(cgrid3d *grid) {

  long i;

  for(i = 0; i < grid->nx * grid->ny * grid->nz; i++)
    grid->value[i] = I * cimag(grid->value[i]);
}

/*
 * Clear imaginary part of complex grid.
 *
 */

EXPORT void cgrid3d_zero_im(cgrid3d *grid) {

  long i;

  for(i = 0; i < grid->nx * grid->ny * grid->nz; i++)
    grid->value[i] = creal(grid->value[i]);
}

/*
 * Extract complex phase factors from a given grid.
 *
 * src = Source grid (cgrid3d *).
 * dst = Dest grid (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_phase(rgrid3d *dst, cgrid3d *src) {

  long i;

  if(dst->nx != src->nx || dst->ny != src->ny || dst->nz != src->nz) {
    fprintf(stderr, "libgrid: incompatible dimensions in cgrid3d_phase().\n");
    exit(1);
  }

  for (i = 0; i < dst->nx * dst->ny * dst->nz; i++) {
    if(cabs(src->value[i]) < GRID_EPS) dst->value[i] = 0.0;
    else dst->value[i] = creal(-I * clog(src->value[i] / cabs(src->value[i])));
    if(dst->value[i] < 0.0) dst->value[i] = 2.0 * M_PI + dst->value[i];
  }
}

/*
 * Add random noise to grid.
 *
 * grid  = Grid where the noise will be added (cgrid3d *).
 * scale = Scaling for random numbers [-1,+1[ (double).
 *
 */

EXPORT void cgrid3d_random(cgrid3d *grid, double scale) {

  static int been_here = 0;
  long i;

  if(!been_here) {
    srand48(time(0));
    been_here = 1;
  }
  for (i = 0; i < grid->nx * grid->ny * grid->nz; i++)
    grid->value[i] += scale * (2.0 * (drand48() - 0.5) + 2.0 * (drand48() - 0.5) * I);
}
