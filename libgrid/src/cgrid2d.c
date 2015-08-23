/*
 * Routines for 2D complex grids.
 *
 * Nx is major index and Ny is minor index (varies most rapidly).
 *
 * To debug memory allocation, define __DEBUG__.
 *
 * TODO: Center of grid stuff (see 3D).
 *
 */

#include "grid.h"
#include "private.h"
#include "private2d.h"

#ifdef __DEBUG__
static int allocated_grids = 0;
#endif /* __DEBUG__ */

/*
 * Allocate 2D grid.
 *
 * nx                 = number of points on the grid along x (long).
 * ny                 = number of points on the grid along y (long).
 * step               = spatial step length on the grid (double).
 * value_outside      = condition for accessing boundary points:
 *                      CGRID2D_DIRICHLET_BOUNDARY: Dirichlet boundary.
 *                      or CGRID2D_NEUMANN_BOUNDARY: Neumann boundary.
 *                      or CGRID2D_PERIODIC_BOUNDARY: Periodic boundary.
 *                      or user supplied function with pointer to grid and
 *                         grid index as parameters to provide boundary access.
 * outside_params_ptr = pointer for passing parameters for the given boundary
 *                      access function. Use 0 to with the predefined boundary
 *                      functions (CGRID2D_* above).
 *
 * Return value: pointer to the allocated grid (cgrid2d *). Returns NULL on
 * error.
 *
 */

EXPORT cgrid2d *cgrid2d_alloc(long nx, long ny, double step, double complex (*value_outside)(const cgrid2d *grid, long i, long j), void *outside_params_ptr) {

  cgrid2d *grid;
  
  grid = (cgrid2d *) malloc(sizeof(cgrid2d));
  if (!grid) {
    fprintf( stderr, "libgrid: Error in cgrid2d_alloc(). Could not allocate memory for 2d grid.\n");
    return 0;
  }
  
    if (!(grid->value = (double complex *) fftw_malloc(nx * ny * sizeof(double complex)))) {
    fprintf(stderr, "libgrid: Error in cgrid2d_alloc(). Could not allocate memory for cgrid2d->value.\n");
    free(grid);
    return 0;
  }
  
  grid->nx = nx;
  grid->ny = ny;
  grid->step = step;
  
  if (value_outside)
    grid->value_outside = value_outside;
  else
    grid->value_outside = cgrid2d_value_outside_constantdirichlet;
  
  if (outside_params_ptr)
    grid->outside_params_ptr = outside_params_ptr;
  else {
    grid->default_outside_params = 0.0;
    grid->outside_params_ptr = &grid->default_outside_params;
  }
  
  grid->plan = grid->iplan = NULL;
  
#if __DEBUG__
  allocated_grids++;
  fprintf(stderr, "libmeas(debug): %3d 2d grids allocated.\n", allocated_grids);
#endif

  cgrid2d_constant(grid, NAN);
  
  return grid;
}

/*
 * Free 2D grid.
 *
 * grid = pointer to 2D grid to be freed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_free(cgrid2d *grid) {

  if (grid) {
    if (grid->value) fftw_free(grid->value);
    cgrid2d_fftw_free(grid);
    free(grid);
  }
}

/* 
 * Write 2D grid on disk in binary format.
 *
 * grid = 2D grid to be written (cgrid2d *).
 * out  = file handle for the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_write(cgrid2d *grid, FILE *out) {

  fwrite(&grid->nx, sizeof(long), 1, out);
  fwrite(&grid->ny, sizeof(long), 1, out);
  fwrite(&grid->step, sizeof(double), 1, out);
  fwrite(grid->value, sizeof(double complex), grid->nx * grid->ny, out);
}

/* 
 * Read 2D grid from disk in binary format.
 *
 * grid = 2D grid to be read (cgrid2d *).
 * in   = file handle for reading the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 * Note: This works for both Cartesian and cylindrical grids.
 *
 */

EXPORT void cgrid2d_read(cgrid2d *grid, FILE *in) {

  long nx, ny;
  double step;
  
  fread(&nx, sizeof(long), 1, in);
  fread(&ny, sizeof(long), 1, in);
  fread(&step, sizeof(double), 1, in);
  
  if (nx != grid->nx || ny != grid->ny || step != grid->step) {
    cgrid2d *tmp;

    fprintf(stderr, "libgrid: Grid in file has different size than grid in memory.\n");
    fprintf(stderr, "libgrid: Interpolating between grids.\n");
    if(!(tmp = cgrid2d_alloc(nx, ny, step, grid->value_outside, NULL))) {
      fprintf(stderr, "libgrid: Error allocating grid in cgrid2d_read().\n");
      abort();
    }
    fread(tmp->value, sizeof(double complex), nx * ny, in);
    cgrid2d_extrapolate(grid, tmp);
    cgrid2d_free(tmp);
    return;
  }
  
  fread(grid->value, sizeof(double complex), grid->nx * grid->ny, in);
}

/*
 * Copy 2D grid from one grid to another.
 *
 * copy = destination grid (cgrid2d *).
 * grid = source grid (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_copy(cgrid2d *copy, const cgrid2d *grid) {

  long i, nx = grid->nx, ny = grid->ny, bytes = grid->ny * sizeof(double complex);
  double complex *gvalue = grid->value;
  double complex *cvalue = copy->value;
  
  copy->nx = grid->nx;
  copy->ny = grid->ny;
  copy->step = grid->step;
  
#pragma omp parallel for firstprivate(nx,ny,bytes,gvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    memmove(&cvalue[i * ny], &gvalue[i * ny], bytes);
}

/*
 * Take complex conjugate of 2D grid.
 * 
 * conjugate = destination for complex conjugated grid (cgrid2d *).
 * grid      = source grid for the operation (cgrid2d *).
 *
 * No return value.
 *
 * Note: source and destination may be the same grid.
 * 
 */

EXPORT void cgrid2d_conjugate(cgrid2d *conjugate, const cgrid2d *grid) {

  long ij, nxy = grid->nx * grid->ny;
  double complex *cvalue = conjugate->value;
  double complex *gvalue = grid->value;
  
  conjugate->nx = grid->nx;
  conjugate->ny = grid->ny;
  conjugate->step = grid->step;
  
#pragma omp parallel for firstprivate(nxy,cvalue,gvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = conj(gvalue[ij]);
}

/*
 * Shift 2D grid by given amount spatially.
 *
 * shifted = destination grid for the operation (cgrid2d *).
 * grid    = source grid for the operation (cgrid2d *).
 * x       = shift spatially by this amount in x (double).
 * y       = shift spatially by this amount in y (double).
 *
 * No return value.
 *
 * NOTE: Source and destination may be the same grid.
 *
 */

EXPORT void cgrid2d_shift(cgrid2d *shifted, const cgrid2d *grid, double x, double y) {

  sShiftParametersc2d params;

  /* shift by (x,y) i.e. current grid center to (x,y) */
  params.x = x;  params.y = y; params.grid = grid;
  cgrid2d_map(shifted, shift_cgrid2d, &params);
}

/* 
 * Zero 2D grid.
 *
 * grid = grid to be zeroed (cgrid2d *).
 *
 * No return value.
 * 
 */

EXPORT void cgrid2d_zero(cgrid2d *grid) { 

  cgrid2d_constant(grid, 0.0); 
}

/* 
 * Set 2D grid to a constant value.
 *
 * grid = grid to be set (cgrid2d *).
 * c    = value (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_constant(cgrid2d *grid, double complex c) {

  long ij, nxy = grid->nx * grid->ny;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,c) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] = c;
}

/*
 * Multiply a given grid by a function.
 *
 * grid = destination grid for the operation (cgrid2d *).
 * func = function providing the mapping (double complex (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and x,y are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_product_func(cgrid2d *grid, double complex (*func)(void *arg, double x, double y), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double x, y, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,step,func,value) private(i,j,ij,x,y) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    x = (i - nx/2.0) * step;
    y = (j - ny/2.0) * step;
    value[ij] *= func(farg, x, y);
  }
}

/*
 * Map a given function onto 2D grid.
 *
 * grid = destination grid for the operation (cgrid2d *).
 * func = function providing the mapping (double complex (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and x,y are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_map(cgrid2d *grid, double complex (*func)(void *arg, double x, double y), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double x, y, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,step,func,value) private(i,j,ij,x,y) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    x = (i - nx/2.0) * step;
    y = (j - ny/2.0) * step;
    value[ij] = func(farg, x, y);
  }
}

/*
 * Map a given function onto 2D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid.
 * *
 * grid = destination grid for the operation (cgrid2d *).
 * func = function providing the mapping (double complex (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and (x,y) is the point (doubles) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 * ns   = number of intermediate points to be used in smoothing (int).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_smooth_map(cgrid2d *grid, double complex (*func)(void *arg, double x, double y), void *farg, int ns) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double xc, yc, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,ns,step,func,value) private(i,j,xc,yc) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;    
    xc = (i - nx/2.0) * step;
    yc = (j - ny/2.0) * step;
    value[ij] = linearly_weighted_integralc2d(func, farg, xc, yc, step, ns);
  }
}

/*
 * Map a given function onto 2D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid. Limits for intermediate steps and
 * tolerance can be given.
 *
 * grid   = destination grid for the operation (cgrid2d *).
 * func   = function providing the mapping (double complex (*)(void *, double, double, double)).
 *          The first argument (void *) is for external user specified data
 *          and (x,y) is the point (doubles) where the function is evaluated.
 * farg   = pointer to user specified data (void *).
 * min_ns = minimum number of intermediate points to be used in smoothing (int).
 * max_ns = maximum number of intermediate points to be used in smoothing (int).
 * tol    = tolerance for weighing (double).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_adaptive_map(cgrid2d *grid, double complex (*func)(void *arg, double x, double y), void *farg, int min_ns, int max_ns, double tol) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny, ns;
  double xc, yc, step = grid->step;
  double tol2 = tol * tol;
  double complex sum, sump;
  double complex *value = grid->value;
  
  if (min_ns < 1) min_ns = 1;
  if (max_ns < min_ns) max_ns = min_ns;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,min_ns,max_ns,step,func,value,tol2) private(i,j,ns,xc,yc,sum,sump) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;    
    xc = (i - nx/2.0) * step;
    yc = (j - ny/2.0) * step;          
    sum  = func(farg, xc, yc); sump = 0.0;
    for(ns = min_ns; ns <= max_ns; ns *= 2) {
      sum  = linearly_weighted_integralc2d(func, farg, xc, yc, step, ns);
      sump = linearly_weighted_integralc2d(func, farg, xc, yc, step, ns+1);
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
 * Add two 2D grids ("gridc = grida + gridb").
 *
 * gridc = destination grid (cgrid2d *).
 * grida = 1st of the grids to be added (cgrid2d *).
 * gridb = 2nd of the grids to be added (cgrid2d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid2d_sum(cgrid2d *gridc, const cgrid2d *grida, const cgrid2d *gridb) {

  long ij, nxy = gridc->nx * gridc->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = avalue[ij] + bvalue[ij];
}

/* 
 * Subtract two grids ("gridc = grida - gridb").
 *
 * gridc = destination grid (cgrid2d *).
 * grida = 1st source grid (cgrid2d *).
 * gridb = 2nd source grid (cgrid2d *).
 *
 * No return value.
 *
 * Note: both source and destination may be the same.
 *
 */

EXPORT void cgrid2d_difference(cgrid2d *gridc, const cgrid2d *grida, const cgrid2d *gridb) {

  long ij, nxy = gridc->nx * gridc->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = avalue[ij] - bvalue[ij];
}

/* 
 * Calculate product of two grids ("gridc = grida * gridb").
 *
 * gridc = destination grid (cgrid2d *).
 * grida = 1st source grid (cgrid2d *).
 * gridb = 2nd source grid (cgrid2d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid2d_product(cgrid2d *gridc, const cgrid2d *grida, const cgrid2d *gridb) {

  long ij, nxy = gridc->nx * gridc->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = avalue[ij] * bvalue[ij];
}

/* 
 * Rise a grid to given power.
 *
 * gridb    = destination grid (cgrid2d *).
 * grida    = 1st source grid (cgrid2d *).
 * exponent = exponent to be used (double).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *       This routine uses pow() so that the exponent can be
 *       fractional but this is slow! Do not use this for integer
 *       exponents.
 *
 * TODO: Add complex exponentiation later.
 *
 */

EXPORT void cgrid2d_power(cgrid2d *gridb, const cgrid2d *grida, double exponent) {

  long ij, nxy = gridb->nx * gridb->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,exponent) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    bvalue[ij] = pow(avalue[ij], exponent);
}

/* 
 * Rise absolute value of a grid to given power.
 *
 * gridb    = destination grid (cgrid2d *).
 * grida    = 1st source grid (cgrid2d *).
 * exponent = exponent to be used (double).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *       This routine uses pow() so that the exponent can be
 *       fractional but this is slow! Do not use this for integer
 *       exponents.
 *
 * TODO: Add complex exponentiation later.
 *
 */

EXPORT void cgrid2d_abs_power(cgrid2d *gridb, const cgrid2d *grida, double exponent) {

  long ij, nxy = gridb->nx * gridb->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,exponent) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    bvalue[ij] = pow(cabs(avalue[ij]), exponent);
}

/*
 * Divide two grids ("gridc = grida / gridb").
 *
 * gridc = destination grid (cgrid2d *).
 * grida = 1st source grid (cgrid2d *).
 * gridb = 2nd source grid (cgrid2d *).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same.
 *
 */

EXPORT void cgrid2d_division(cgrid2d *gridc, const cgrid2d *grida, const cgrid2d *gridb) {
  
  long ij, nxy = gridc->nx * gridc->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = avalue[ij] / bvalue[ij];
}

/* 
 * Conjugate product of two grids ("gridc = conj(grida) * gridb").
 *
 * gridc = destination grid (cgrid2d *).
 * grida = 1st source grid (cgrid2d *).
 * gridb = 2nd source grid (cgrid2d *).
 *
 * No return value.
 * 
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid2d_conjugate_product(cgrid2d *gridc, const cgrid2d *grida, const cgrid2d *gridb) {

  long ij, nxy = gridc->nx * gridc->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = conj(avalue[ij]) * bvalue[ij];
}

/*
 * Add a constant to a 2D grid.
 *
 * grid = grid where the constant is added (cgrid2d *).
 * c    = constant to be added (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_add(cgrid2d *grid, double complex c) {

  long ij, nxy = grid->nx * grid->ny;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,c) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] += c;
}

/*
 * Multiply grid by a constant.
 *
 * grid = grid to be multiplied (cgrid2d *).
 * c    = multiplier (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_multiply(cgrid2d *grid, double complex c) {

  long ij, nxy = grid->nx * grid->ny;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,c) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] *= c;
}

/* 
 * Add and multiply: grid = (grid + ca) * cm.
 *
 * grid = grid to be operated (cgrid2d *).
 * ca   = constant to be added (double complex).
 * cm   = multiplier (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_add_and_multiply(cgrid2d *grid, double complex ca, double complex cm) {

  long ij, nxy = grid->nx * grid->ny;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,ca,cm) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] = (value[ij] + ca) * cm;
}

/*
 * Multiply and add: grid = cm * grid + ca.
 *
 * grid = grid to be operated (cgrid2d *).
 * ca   = constant to be added (double complex).
 * cm   = multiplier (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_multiply_and_add(cgrid2d *grid, double complex cm, double complex ca) {

  long ij, nxy = grid->nx * grid->ny;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,ca,cm) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] = value[ij] * cm + ca;
}

/* 
 * Add scaled grid: gridc = gridc + d * grida
 *
 * gridc = destination grid for the operation (cgrid2d *).
 * d     = multiplier for the operation (double complex).
 * grida = source grid for the operation (cgrid2d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid2d_add_scaled(cgrid2d *gridc, double complex d, const cgrid2d *grida) {

  long ij, nxy = gridc->nx * gridc->ny;
  double complex *avalue = grida->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(d,nxy,avalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] += d * avalue[ij];
}

/*
 * Add multiply two grids and a constant: gridc = gridc + c * grida * gridb.
 *
 * gridc = destination grid (cgrid2d *).
 * d     = constant multiplier (double complex).
 * grida = 1st source grid (cgrid2d *).
 * gridb = 2nd source grid (cgrid2d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid2d_add_scaled_product(cgrid2d *gridc, double complex d, const cgrid2d *grida, const cgrid2d *gridb) {
  
  long ij, nxy = gridc->nx * gridc->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue,d) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] += d * avalue[ij] * bvalue[ij];
}

/*
 * Operate on a grid by a given operator: gridc = O(grida).
 *
 * gridc    = destination grid (cgrid2d *).
 * grida    = source grid (cgrid2d *).
 * operator = operator (double complex (*)(double complex)).
 *            (i.e. a function mapping a given C-number to another)
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid2d_operate_one(cgrid2d *gridc, const cgrid2d *grida, double complex (*operator)(double complex a)) {

  long ij, nxy = gridc->nx * gridc->ny;
  double complex *avalue = grida->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,cvalue,operator) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = operator(avalue[ij]);
}

/* 
 * Operate on two grids and place the result in third: gridc = O(grida, gridb).
 * where O is the operator.
 *
 * gridc    = destination grid (cgrid2d *).
 * grida    = 1s source grid (cgrid2d *).
 * gridb    = 2nd source grid (cgrid2d *).
 * operator = operator mapping grida and gridb (double complex (*)(double
 *            complex, double complex)).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid2d_operate_two(cgrid2d *gridc, const cgrid2d *grida, const cgrid2d *gridb, double complex (*operator)(double complex a, double complex b)) {

  long ij, nxy = gridc->nx * gridc->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue,operator) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = operator(avalue[ij], bvalue[ij]);
}

/*
 * Operate on a grid by a given operator.
 *
 * grid     = grid to be operated (cgrid2d *).
 * operator = operator (void (*)(double complex *)).
 * 
 * No return value.
 *
 */

EXPORT void cgrid2d_transform_one(cgrid2d *grid, void (*operator)(double complex *a)) {

  long ij, nxy = grid->nx * grid->ny;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,operator) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    operator(&value[ij]);
}

/*
 * Operate on two separate grids by a given operator.
 *
 * grida    = grid to be operated (cgrid2d *).
 * gridb    = grid to be operated (cgrid2d *).
 * operator = operator (void (*)(double complex *)).
 * 
 * No return value.
 *
 */

EXPORT void cgrid2d_transform_two(cgrid2d *grida, cgrid2d *gridb, void (*operator)(double complex *a, double complex *b)) {

  long ij, nxy = grida->nx * grida->ny;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,operator) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    operator(&avalue[ij], &bvalue[ij]);
}

/*
 * Integrate over a grid.
 *
 * grid = grid to be integrated (cgrid2d *).
 *
 * Returns the integral value (double complex).
 *
 */

EXPORT double complex cgrid2d_integral(const cgrid2d *grid) {

  long i, j, nx = grid->nx, ny = grid->ny;
  double complex sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,grid) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i <= nx; i++)
    for(j = 0; j <= ny; j++)
      sum += cgrid2d_value_at_index(grid, i, j);  
  return sum * grid->step * grid->step;
}

/*
 * Integrate over a grid with limits.
 *
 * grid = grid to be integrated (cgrid2d *).
 * xl   = lower limit for x (double).
 * xu   = upper limit for x (double).
 * yl   = lower limit for y (double).
 * yu   = upper limit for y (double).
 *
 * Returns the integral value (double complex).
 *
 */

EXPORT double complex cgrid2d_integral_region(const cgrid2d *grid, double xl, double xu, double yl, double yu) {

  long iu, il, i, ju, jl, j, nx = grid->nx, ny = grid->ny;
  double complex sum;
  double step = grid->step;
   
  il = xl / step + nx/2;
  iu = xu / step + nx/2;
  jl = yl / step + ny/2;
  ju = yu / step + ny/2;
  
  sum = 0.0;
#pragma omp parallel for firstprivate(il,iu,jl,ju,nx,ny,grid) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for (i = il; i <= iu; i++)
    for (j = jl; j <= ju; j++)
      sum += cgrid2d_value_at_index(grid, i, j);
  return sum * step * step;
}

/* 
 * Integrate over the grid squared (int |grid|^2).
 *
 * grid = grid to be integrated (cgrid2d *).
 *
 * Returns the integral (double complex).
 *
 */

EXPORT double cgrid2d_integral_of_square(const cgrid2d *grid) {

  long i, j, nx = grid->nx, ny = grid->ny;
  double sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,grid) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i <= nx; i++)
    for(j = 0; j <= ny; j++)
      sum += sqnorm(cgrid2d_value_at_index(grid, i, j));  
  return sum * grid->step * grid->step;
}

/*
 * Calculate overlap between two grids (int grida^*gridb).
 *
 * grida = 1st grid (complex conjugated) (cgrid2d *).
 * gridb = 2nd grid (no complex conjugation) (cgrid2d *).
 *
 * Returns the value of the overlap integral (double complex).
 *
 */

EXPORT double complex cgrid2d_integral_of_conjugate_product(const cgrid2d *grida, const cgrid2d *gridb) {

  long i, j, nx = grida->nx, ny = grida->ny;
  double complex sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,grida,gridb) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i <= nx; i++)
    for(j = 0; j <= ny; j++)
      sum += conj(cgrid2d_value_at_index(grida, i, j)) * cgrid2d_value_at_index(gridb, i, j);
  return sum * grida->step * grida->step;
}

/*
 * Calculate expectation value of a grid over the probability density given by the other.
 * (int grida |gridb|^2).
 *
 * grida = grid giving the probability (|gridb|^2) (cgrid2d *).
 * gridb = grid to be averaged (cgrid2d *).
 *
 * Returns the average value (double complex).
 *
 */

EXPORT double complex cgrid2d_grid_expectation_value(const cgrid2d *grida, const cgrid2d *gridb) {

  long i, j, nx = grida->nx, ny = grida->ny;
  double complex sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,grida,gridb) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      sum += sqnorm(cgrid2d_value_at_index(grida,i,j)) * cgrid2d_value_at_index(gridb,i,j);  
  return sum * grida->step * grida->step;
}

/*
 * Calculate the expectation value of a function over a grid.
 * (int grida^* func grida = int func |grida|^2).
 *
 * func  = function to be averaged (double complex (*)(void *, double complex, double, double)).
 *         The arguments are: optional arg, grida(x,y), x, y.
 * grida = grid giving the probability (|grida|^2) (cgrid2d *).
 *
 * Returns the average value (double complex).
 *
 */
 
EXPORT double complex cgrid2d_grid_expectation_value_func(void *arg, double complex (*func)(void *arg, double complex val, double x, double y), const cgrid2d *grida) {
   
  long i, j, nx = grida->nx, ny = grida->ny;
  double complex sum = 0.0, tmp;
  double x, y, step = grida->step;
  
#pragma omp parallel for firstprivate(func,arg,nx,ny,grida,step) private(x,y,i,j,tmp) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i <= nx; i++) {
    x = (i - nx/2) * step;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step;
      tmp = cgrid2d_value_at_index(grida, i, j);
      sum += sqnorm(tmp) * func(arg, tmp, x, y);
    }
  }
  return sum * step * step;
}

/* 
 * Integrate over the grid multiplied by weighting function (int grid w(x)).
 *
 * grid   = grid to be integrated over (cgrid2d *).
 * weight = function defining the weight (double complex (*)(double, double)). The arguments are (x,y) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double complex).
 *
 */

EXPORT double complex cgrid2d_weighted_integral(const cgrid2d *grid, double complex (*weight)(void *farg, double x, double y), void *farg) {

  long i, j, ny = grid->ny, nx = grid->nx;
  double x, y, step = grid->step;
  double complex sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,grid,step,weight,farg) private(i,j,x,y) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i <= nx; i++) {
    x = (i - nx/2) * step;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step;
      sum += weight(farg, x, y) * cgrid2d_value_at_index(grid, i, j);
    }
  }
  return sum * step * step;
}

/* 
 * Integrate over square of the grid multiplied by weighting function (int grid^2 w(x)).
 *
 * grid   = grid to be integrated over (cgrid2d *).
 * weight = function defining the weight (double complex (*)(double, double)).
 *          The arguments are (x,y) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double complex).
 *
 */

EXPORT double cgrid2d_weighted_integral_of_square(const cgrid2d *grid, double (*weight)(void *farg, double x, double y), void *farg) {

  long i, j, nx = grid->nx, ny = grid->ny;
  double x, y, step = grid->step;
  double complex sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,step,grid,weight,farg) private(i,j,x,y) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {    
    x = (i - nx/2) * step;
    for(j = 0; j < ny; j++) {
      y = (j - ny/2) * step;
      sum += weight(farg, x, y) * sqnorm(cgrid2d_value_at_index(grid, i, j));
    }
  }  
  return sum * step * step;
}

/* 
 * Differentiate a grid with respect to x.
 *
 * grid     = grid to be differentiated (cgrid2d *).
 * gradient = differentiated grid output (cgrid2d *).
 * 
 * No return value.
 *
 */

EXPORT void cgrid2d_fd_gradient_x(const cgrid2d *grid, cgrid2d *gradient) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double complex *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in cgrid2d_fd_gradient_x().\n");
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
 * Differentiate a grid with respect to y.
 *
 * grid     = grid to be differentiated (cgrid2d *).
 * gradient = differentiated grid output (cgrid2d *).
 * 
 * No return value.
 *
 */

EXPORT void cgrid2d_fd_gradient_y(const cgrid2d *grid, cgrid2d *gradient) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double complex *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in cgrid2d_fd_gradient_y().\n");
    return;
  }
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta,grid) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = inv_delta * (cgrid2d_value_at_index(grid, i, j+1) - cgrid2d_value_at_index(grid, i, j-1));
  }
}

/*
 * Calculate gradient of a grid.
 *
 * grid       = grid to be differentiated twice (cgrid2d *).
 * gradient_x = x output grid for the operation (cgrid2d *).
 * gradient_y = y output grid for the operation (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_fd_gradient(const cgrid2d *grid, cgrid2d *gradient_x, cgrid2d *gradient_y) {

  cgrid2d_fd_gradient_x(grid, gradient_x);
  cgrid2d_fd_gradient_y(grid, gradient_y);
}

/*
 * Calculate laplacian of the grid.
 *
 * grid    = source grid (cgrid2d *).
 * laplace = output grid for the operation (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_fd_laplace(const cgrid2d *grid, cgrid2d *laplace) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_delta2 = 1.0 / (grid->step * grid->step);
  double complex *lvalue = laplace->value;
  
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta2,grid) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = inv_delta2 * (-4.0 * cgrid2d_value_at_index(grid, i, j) + cgrid2d_value_at_index(grid, i, j+1) + cgrid2d_value_at_index(grid, i, j-1) + cgrid2d_value_at_index(grid,i+1,j) + cgrid2d_value_at_index(grid,i-1,j));
  }
}

/*
 * Calculate dot product of the gradient of the grid.
 *
 * grid          = source grid for the operation (cgrid2d *).
 * grad_dot_grad = destination grid (cgrid2d *).
 *
 * No return value.
 *
 * Note: grid and grad_dot_grad may not be the same grid.
 *
 */

EXPORT void cgrid2d_fd_gradient_dot_gradient(const cgrid2d *grid, cgrid2d *grad_dot_grad) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_2delta2 = 1.0 / (2.0 * grid->step * 2.0 * grid->step);
  double complex *gvalue = grad_dot_grad->value;
  
/*  grad f(x,y,z) dot grad f(x,y,z) = [ |f(+,0) - f(-,0)|^2 + |f(0,+) - f(0,-)|^2] / (2h)^2 */
#pragma omp parallel for firstprivate(ny,nxy,gvalue,inv_2delta2,grid) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    gvalue[ij] = inv_2delta2 * (sqnorm(cgrid2d_value_at_index(grid, i, j+1) - cgrid2d_value_at_index(grid, i, j-1)) + sqnorm(cgrid2d_value_at_index(grid, i+1, j) - cgrid2d_value_at_index(grid, i-1, j)));
  }
}

/*
 * Print the grid with both real and imaginary parts into file (ASCII format).
 *
 * grid = grid to be printed out (cgrid2d *).
 * out  = output file pointer (FILE *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_print(const cgrid2d *grid, FILE *out) {

  long i, j;

  for( i = 0; i < grid->nx; i++ ) {
    for( j = 0; j < grid->ny; j++ ) {
      fprintf(out, "%16.8le %16.8le   ", 
	      creal(cgrid2d_value_at_index(grid, i, j)),
	      cimag(cgrid2d_value_at_index(grid, i, j)));
    }
    fprintf(out, "\n");
  }
}

/*
 * Perform Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (cgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       Also no normalization is performed.
 *
 */

EXPORT void cgrid2d_fft(cgrid2d *grid) {

  if (!grid->plan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in cgrid2d_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    cgrid2d_fftw_alloc(grid);
  }  
  cgrid2d_fftw(grid);
}

/*
 * Perform normalized Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (cgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void cgrid2d_fourier_transform(cgrid2d *grid) {

  cgrid2d_fft(grid);
  cgrid2d_multiply(grid, (double complex) (grid->step * grid->step));
}

/*
 * Perform inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       No normalization.
 *
 */

EXPORT void cgrid2d_inverse_fft(cgrid2d *grid) {

  if (!grid->iplan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in cgrid2d_inverse_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    cgrid2d_fftw_alloc(grid);
  }
  cgrid2d_fftw_inv(grid);
}


/*
 * Perform scaled inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid2d *).
 * c    = scaling factor (i.e. the output is multiplied by this constant) (double complex).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void cgrid2d_scaled_inverse_fft(cgrid2d *grid, double complex c) {

  cgrid2d_inverse_fft(grid);
  cgrid2d_multiply(grid, c);  
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by 1 / (Nx * Ny).
 * (Ni = number of grid points along axis i)
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void cgrid2d_inverse_fft_norm(cgrid2d *grid) {

  cgrid2d_scaled_inverse_fft(grid, (double complex) ((1.0 / grid->nx) * (1.0 / grid->ny)));
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by 1 / (Nx*step * Ny*step).
 * (Ni = number of grid points along axis i, step = grid step length)
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       This is the inverse routine to be used with cgrid2d_fourier_transform().
 *
 */

EXPORT void cgrid2d_inverse_fourier_transform(cgrid2d *grid) {

  cgrid2d_scaled_inverse_fft(grid, (double complex) (1.0 / ((grid->nx * grid->step) *
							   (grid->ny * grid->step))));
}

/*
 * Convolue FFT transformed grids. To apply this on two grids (grida and gridb)
 * and place the result in gridc:
 * cgrid2d_fft(grida);
 * cgrid2d_fft(gridb);
 * cgrid2d_convolue(gridc, grida, gridb);
 * cgrid2d_inverse_fft(gridc);
 * gridc now contains the convolution of grida and gridb.
 *
 * grida = 1st grid to be convoluted (cgrid2d *).
 * gridb = 2nd grid to be convoluted (cgrid2d *).
 * gridc = output (cgrid2d *).
 *
 * No return value.
 *
 * Note: the input/output grids may be the same.
 *
 */

EXPORT void cgrid2d_fft_convolute(cgrid2d *gridc, const cgrid2d *grida, const cgrid2d *gridb) {

  long i, j, ij, nx, ny, nxy;
  double step, norm;
  double complex *cvalue, *bvalue;
  
  /* int f(r) g(r-r') d^2r' = iF[ F[f] F[g] ] = (step / N)^2 iFFT[ FFT[f] FFT[g] ] */
  nx = gridc->nx;
  ny = gridc->ny;
  nxy = nx * ny;
  step = gridc->step;
  
  if (gridc != grida)
    cgrid2d_copy(gridc, grida);
  
  cvalue = gridc->value;
  bvalue = gridb->value;
  
  norm = (step / nx) * (step / ny);
  
#pragma omp parallel for firstprivate(nx,ny,nxy,cvalue,bvalue,norm) private(i,j,ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    /* if odd */
    if ((i + j) & 1)
      cvalue[ij] = -norm * cvalue[ij] * bvalue[ij];
    else
      cvalue[ij] = norm * cvalue[ij] * bvalue[ij];
  }
}

/*
 * Differentiate grid in the Fourier space along x.
 *
 * grid       = grid to be differentiated (in Fourier space) (cgrid2d *).
 * gradient_x = output grid (cgrid2d *).
 *
 * No return value.
 *
 * Note: input and output grids may be the same.
 *
 */

EXPORT void cgrid2d_fft_gradient_x(const cgrid2d *grid, cgrid2d *gradient_x) {

  long i, ij, nx, ny, nxy;
  double kx, step, norm;
  double complex *gxvalue = gradient_x->value;
  
  /* f'(x) = iF[ i kx F[f(x)] ] */  
  nx = grid->nx;
  ny = grid->ny;
  nxy = nx * ny;
  step = grid->step;
  
  norm = (1.0/ nx) * (1.0/ ny);
  
  if (gradient_x != grid)
    cgrid2d_copy(gradient_x, grid);
  
#pragma omp parallel for firstprivate(norm,nx,ny,nxy,step,gxvalue) private(i,ij,kx) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (i < nx / 2)
      kx = 2.0 * M_PI * i / (nx * step);
    else 
      kx = 2.0* M_PI * (i - nx) / (nx * step);
    
    gxvalue[ij] *= kx * norm * I;
  }
}

/*
 * Differentiate grid in the Fourier space along y.
 *
 * grid       = grid to be differentiated (in Fourier space) (cgrid2d *).
 * gradient_y = output grid (cgrid2d *).
 *
 * No return value.
 *
 * Note: input and output grids may be the same.
 *
 */

EXPORT void cgrid2d_fft_gradient_y(const cgrid2d *grid, cgrid2d *gradient_y) {

  long j, ij, nx, ny, nxy;
  double ky, step, norm;
  double complex *gyvalue = gradient_y->value;
  
  /* f'(y) = iF[ i ky F[f(y)] ] */  
  nx = grid->nx;
  ny = grid->ny;
  nxy = nx * ny;
  step = grid->step;
  
  norm = (1.0 / nx) * (1.0 / ny);
  
  if (gradient_y != grid)
    cgrid2d_copy(gradient_y, grid);
  
#pragma omp parallel for firstprivate(norm,nx,ny,nxy,step,gyvalue) private(j,ij,ky) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    j = ij % ny;    
    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (j < ny / 2)
      ky = 2.0 * M_PI * j / (ny * step);
    else 
      ky = 2.0 * M_PI * (j - ny) / (ny * step);    
    gyvalue[ij] *= ky * norm * I;
  }
}

/* 
 * Calculate gradient of a grid (in Fourier space).
 *
 * grid       = grid to be differentiated (cgrid2d *).
 * gradient_x = x output grid (cgrid2d *).
 * gradient_y = y output grid (cgrid2d *).
 *
 * No return value.
 *
 * Note: input/output grids may be the same.
 *
 */

EXPORT void cgrid2d_fft_gradient(const cgrid2d *grid, cgrid2d *gradient_x, cgrid2d *gradient_y) {

  cgrid2d_fft_gradient_x(grid, gradient_x);
  cgrid2d_fft_gradient_y(grid, gradient_y);
}

/* 
 * Calculate second derivative of a grid (in Fourier space).
 *
 * grid    = grid to be differentiated (cgrid2d *).
 * laplace = output grid (cgrid2d *).
 *
 * No return value.
 *
 * Note: input/output grids may be the same.
 *
 */

EXPORT void cgrid2d_fft_laplace(const cgrid2d *grid, cgrid2d *laplace)  {

  long i, j, ij, nx, ny, nxy;
  double kx, ky, step, norm;
  double complex *lvalue = laplace->value;
  
  /* f''(x) = iF[ -k^2 F[f(x)] ] */
  nx = grid->nx;
  ny = grid->ny;
  nxy = nx * ny;
  step = grid->step;
  
  norm = (1.0 / nx) * (1.0 / ny);
  
  if (grid != laplace)
    cgrid2d_copy(laplace, grid);
  
#pragma omp parallel for firstprivate(norm,nx,ny,nxy,step,lvalue) private(i,j,ij,kx,ky) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;    
    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (i < nx / 2)
      kx = 2.0 * M_PI * i / (nx * step);
    else 
      kx = 2.0 * M_PI * (i - nx) / (nx * step);
    
    if (j < ny / 2)
      ky = 2.0 * M_PI * j / (ny * step);
    else 
      ky = 2.0 * M_PI * (j - ny) / (ny * step);
    
    lvalue[ij] *= (-kx*kx -ky*ky) * norm;
  }
}

/*
 * Calculate expectation value of laplace operator in the Fourier space (int grid^* grid'').
 *
 * grid    = source grid for the operation (in Fourier space) (cgrid2d *).
 * laplace = laplacian of the grid (input) (cgrid2d *).
 *
 * Returns the expectation value (double).
 *
 */

EXPORT double cgrid2d_fft_laplace_expectation_value(const cgrid2d *grid, cgrid2d *laplace)  {

  long i, j, ij, nx, ny, nxy;
  double kx, ky, step, norm, sum = 0.0;
  double complex *lvalue = laplace->value;
  
  /* int f*(x) f''(x) dx */
  nx = grid->nx;
  ny = grid->ny;
  nxy = nx * ny;
  step = grid->step;
  
  /* int ( delta FFT[f(x)] )^2 dk => delta^2 / N delta */
  norm = (step / nx) * (step / ny);
  
  if (grid != laplace)
    cgrid2d_copy(laplace, grid);
  
#pragma omp parallel for firstprivate(norm,nx,ny,nxy,step,lvalue) private(i,j,ij,kx,ky) reduction(+:sum) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (i < nx / 2)
      kx = 2.0 * M_PI * i / (nx * step);
    else 
      kx = 2.0 * M_PI * (i - nx) / (nx * step);
    
    if (j < ny / 2)
      ky = 2.0 * M_PI * j / (ny * step);
    else 
      ky = 2.0 * M_PI * (j - ny) / (ny * step);
    
    sum += (- kx*kx - ky*ky) * sqnorm(lvalue[ij]);
  }
  
  return sum * norm;
}

/* Boundary condition routines */

EXPORT double complex cgrid2d_value_outside_constantdirichlet(const cgrid2d *grid, long i, long j) {

  return *((double complex *) grid->outside_params_ptr);
}

EXPORT double complex cgrid2d_value_outside_neumann(const cgrid2d *grid, long i, long j) {

  long nx = grid->nx, ny = grid->ny;

  /* The new code corresponds to continuum boundaries */
  if(i < 0) i = 1;
  if(i > nx-1) i = nx - 2;
  if(j < 0) j = 1;
  if(j > ny-1) j = ny - 2;
  
  return grid->value[i * ny + j];  
  /* This is the old mirror Neumann code */
#if 0  

  nd = nx * 2;
  if (i < 0) i  = -i;
  if (i >= nd) i %= nd;
  if (i >= nx) i  = nd - i;
  
  nd = ny * 2;
  if (j < 0) j  = - j;
  if (j >= nd) j %= nd;
  if (j >= ny) j  = nd - j;
  
  return grid->value[i * ny + j];
#endif
}

EXPORT double complex cgrid2d_value_outside_periodic(const cgrid2d *grid, long i, long j) {

  long nx = grid->nx, ny = grid->ny;
  
  i %= nx;
  if (i < 0) i = nx + i;
  j %= ny;
  if (j < 0) j = ny + j;
  
  return grid->value[i * ny + j];  
}

/* End boundary condition routines */

/*
 * Access grid point at given index.
 *
 * grid = grid to be accessed (cgrid2d *).
 * i    = index along x (long).
 * j    = index along y (long).
 *
 * Returns grid value at index (i, j).
 *
 */

EXPORT inline double complex cgrid2d_value_at_index(const cgrid2d *grid, long i, long j) {

  if (i < 0 || j < 0 || i >= grid->nx || j >= grid->ny)
    return grid->value_outside(grid, i, j);
  return grid->value[i * grid->ny + j];
}

/*
 * Access grid point at given (x,y) point (with linear interpolation).
 *
 * grid = grid to be accessed (cgrid2d *).
 * x    = x value (double).
 * y    = y value (double).
 *
 * Returns grid value at (x,y).
 *
 */

EXPORT inline double complex cgrid2d_value(const cgrid2d *grid, double x, double y) {

  double complex f00, f10, f01, f11;
  double omx, omy, step = grid->step;
  long i, j;
  
  /* i to index and 0 <= x < 1 */
  i = (x /= step);
  if (x < 0) i--;
  x -= i;
  i += grid->nx / 2;
  
  /* j to index and 0 <= y < 1 */
  j = (y /= step);
  if (y < 0) j--;
  y -= j;
  j += grid->ny / 2;
  
  /* linear extrapolation 
   *
   * f(x,y) = (1-x) (1-y) f(0,0) + x (1-y)f(1,0) + (1-x) y f(0,1)
   *          + x     y   f(1,1)
   */ 
  f00 = cgrid2d_value_at_index(grid, i, j);
  f10 = cgrid2d_value_at_index(grid, i+1, j);
  f01 = cgrid2d_value_at_index(grid, i, j+1);
  f11 = cgrid2d_value_at_index(grid, i+1, j+1);
  
  omx = 1 - x;
  omy = 1 - y;

  return omx * omy * f00 + x * omy * f10 + omx * y * f01 + x * y * f11;
}

/*
 * Extrapolate between two different grid sizes.
 *
 * dest = Destination grid (cgrid2d *; output).
 * src  = Source grid (cgrid2d *; input).
 *
 */

EXPORT void cgrid2d_extrapolate(cgrid2d *dest, cgrid2d *src) {

  long i, j, nx = dest->nx, ny = dest->ny;
  double step = dest->step, x, y;

  for (i = 0; i < nx; i++) {
    x = (i - nx/2) * step;
    for (j = 0; j < ny; j++) {
      y = (j - ny/2.0) * step;
      dest->value[i * ny + j] = cgrid2d_value(src, x, y);
    }
  }
}

/*
 * Clear real part of complex grid.
 *
 */

EXPORT void cgrid2d_zero_re(cgrid2d *grid) {

  long i;

  for(i = 0; i < grid->nx * grid->ny; i++)
    grid->value[i] = I * cimag(grid->value[i]);
}

/*
 * Clear imaginary part of complex grid.
 *
 */

EXPORT void cgrid2d_zero_im(cgrid2d *grid) {

  long i;

  for(i = 0; i < grid->nx * grid->ny; i++)
    grid->value[i] = creal(grid->value[i]);
}
