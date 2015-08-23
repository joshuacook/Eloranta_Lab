/*
 * Routines for 1D real grids.
 *
 * To debug memory allocation, define __DEBUG__.
 *
 */

#include "grid.h"
#include "private.h"
#include "private1d.h"

#ifdef __DEBUG__
static int allocated_grids = 0;
#endif /* __DEBUG__ */

/*
 * Allocate 1D grid.
 *
 * nx                 = number of points on the grid (long).
 * step               = spatial step length on the grid (double).
 * value_outside      = condition for accessing boundary points:
 *                      RGRID1D_DIRICHLET_BOUNDARY: Dirichlet boundary.
 *                      or RGRID1D_NEUMANN_BOUNDARY: Neumann boundary.
 *                      or RGRID1D_PERIODIC_BOUNDARY: Periodic boundary.
 *                      or user supplied function with pointer to grid and
 *                         grid index as parameters to provide boundary access.
 * outside_params_ptr = pointer for passing parameters for the given boundary
 *                      access function. Use 0 to with the predefined boundary
 *                      functions (RGRID1D_* above).
 *
 * Return value: pointer to the allocated grid (rgrid1d *). Returns NULL on
 * error.
 *
 */

EXPORT rgrid1d *rgrid1d_alloc(long nx, double step, double (*value_outside)(const rgrid1d *grid, long i), void *outside_params_ptr) {

  rgrid1d *grid;
  
  grid = (rgrid1d *) malloc(sizeof(rgrid1d));
  if (!grid) {
    fprintf(stderr, "libgrid: Error in rgrid1d_alloc(). Could not allocate memory for 1d grid.\n");
    return 0;
  }
  
  if(!(grid->value = (double *) fftw_malloc(2 * (nx/2 + 1) * sizeof(double)))) {
    fprintf(stderr, "libgrid: Error in rgrid1d_alloc(). Could not allocate memory for rgrid1d->value.\n");
    free(grid);
    return 0;
  }
  if(!(grid->cint = (cgrid1d *) malloc(sizeof(cgrid1d)))) {
    fprintf(stderr, "libtrid: Error in rgrid1d_alloc(). Could not allocate memory for rgrid1d->cint->value.\n");
    free(grid);
  }
  
  grid->nx = nx;
  grid->step = step;
  /* complex structure interface (after FFT) */
  grid->cint->value = (double complex *) grid->value;
  grid->cint->step = grid->step;
  grid->cint->nx = grid->nx/2 + 1;
  /* value outside not set */
  
  if (value_outside)
    grid->value_outside = value_outside;
  else
    grid->value_outside = rgrid1d_value_outside_constantdirichlet;
  
  if (outside_params_ptr)
    grid->outside_params_ptr = outside_params_ptr;
  else {
    grid->default_outside_params = 0.0;
    grid->outside_params_ptr = &grid->default_outside_params;
  }
  
  grid->plan = grid->iplan = NULL;
  
#if __DEBUG__
  allocated_grids++;
  fprintf(stderr, "libgrid(debug): %3d 1d grids allocated.\n", allocated_grids);
#endif

  rgrid1d_constant(grid, NAN);
  
  return grid;
}

/*
 * Free 1D grid.
 *
 * grid = pointer to 1D grid to be freed (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_free(rgrid1d *grid) {

  if (grid) {
    if (grid->value) fftw_free(grid->value);
    if (grid->cint) free(grid->cint);
    rgrid1d_fftw_free(grid);
    free(grid);
  }
}

/* 
 * Write 1D grid on disk in binary format.
 *
 * grid = 1D grid to be written (rgrid1d *).
 * out  = file handle for the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_write(rgrid1d *grid, FILE *out) {

  fwrite(&grid->nx, sizeof(long), 1, out);
  fwrite(&grid->step, sizeof(double), 1, out);
  fwrite(grid->value, sizeof(double), grid->nx, out);
}

/* 
 * Read 1D grid from disk in binary format.
 *
 * grid = 1D grid to be read (rgrid1d *).
 * in   = file handle for reading the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_read(rgrid1d *grid, FILE *in) {

  long nx;
  double step;
  
  fread(&nx, sizeof(long), 1, in);
  fread(&step, sizeof(double), 1, in);
  
  if (nx != grid->nx || step != grid->step) {
    rgrid1d *tmp;

    fprintf(stderr, "libgrid: Grid in file has different size than grid in memory.\n");
    fprintf(stderr, "libgrid: Interpolating between grids.\n");
    if(!(tmp = rgrid1d_alloc(nx, step, grid->value_outside, NULL))) {
      fprintf(stderr, "libgrid: Error allocating grid in rgrid1d_read().\n");
      abort();
    }
    fread(tmp->value, sizeof(double complex), nx, in);
    rgrid1d_extrapolate(grid, tmp);
    rgrid1d_free(tmp);
    return;
  }
  
  fread(grid->value, sizeof(double), grid->nx, in);
}

/*
 * Copy 1D grid from one grid to another.
 *
 * copy = destination grid (rgrid1d *).
 * grid = source grid (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_copy(rgrid1d *copy, const rgrid1d *grid) {

  long i, nx = grid->nx;
  double *gvalue = grid->value;
  double *cvalue = copy->value;
  
  copy->nx = grid->nx;
  copy->step = grid->step;
  
#pragma omp parallel for firstprivate(nx,gvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = gvalue[i];
}

/*
 * Shift 1D grid by given amount spatially.
 *
 * shifted = destination grid for the operation (rgrid1d *).
 * grid    = source grid for the operation (rgrid1d *).
 * x       = shift spatially by this amount (double).
 *
 * No return value.
 *
 * NOTE: Source and destination may be the same grid.
 *
 */

EXPORT void rgrid1d_shift(rgrid1d *shifted, const rgrid1d *grid, double x) {

  sShiftParametersr1d params;

  /* shift by (x) i.e. current grid center to (x) */
  params.x = x; 
  params.grid = grid;
  rgrid1d_map(shifted, shift_rgrid1d, &params);
}

/* 
 * Zero 1D grid.
 *
 * grid = grid to be zeroed (rgrid1d *).
 *
 * No return value.
 * 
 */

EXPORT void rgrid1d_zero(rgrid1d *grid) { 

  rgrid1d_constant(grid, 0.0); 
}

/* 
 * Set 1D grid to a constant value.
 *
 * grid = grid to be set (rgrid1d *).
 * c    = value (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_constant(rgrid1d *grid, double c) {

  long i, nx = grid->nx;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(value,c,nx) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = c;
}

/*
 * Multiply a grid by a given function.
 *
 * grid = destination grid for the operation (rgrid1d *).
 * func = function providing the mapping (double (*)(void *, double)).
 *        The first argument (void *) is for external user specified data
 *        and x is the coordinate (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_product_func(rgrid1d *grid, double (*func)(void *arg, double x), void *farg) {

  long i, nx = grid->nx;
  double x, step = grid->step;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,step,func,value) private(i,x) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    value[i] *= func(farg, x);
  }
}

/*
 * Map a given function onto 1D grid.
 *
 * grid = destination grid for the operation (rgrid1d *).
 * func = function providing the mapping (double (*)(void *, double)).
 *        The first argument (void *) is for external user specified data
 *        and x is the coordinate (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_map(rgrid1d *grid, double (*func)(void *arg, double x), void *farg) {

  long i, nx = grid->nx;
  double x, step = grid->step;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,step,func,value) private(i,x) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    value[i] = func(farg, x);
  }
}

/*
 * Map a given function onto 1D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid.
 * *
 * grid = destination grid for the operation (rgrid1d *).
 * func = function providing the mapping (double (*)(void *, double)).
 *        The first argument (void *) is for external user specified data
 *        and x is the coordinate (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 * ns   = number of intermediate points to be used in smoothing (int).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_smooth_map(rgrid1d *grid, double (*func)(void *arg, double x), void *farg, int ns) {

  long i, nx = grid->nx;
  double xc, step = grid->step;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ns,step,func,value) private(i,xc) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {
    xc = (i - nx/2.0) * step;
    value[i] = linearly_weighted_integralr1d(func, farg, xc, step, ns);
  }
}

/*
 * Map a given function onto 1D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid. Limits for intermediate steps and
 * tolerance can be given.
 *
 * grid   = destination grid for the operation (rgrid1d *).
 * func   = function providing the mapping (double (*)(void *, double)).
 *          The first argument (void *) is for external user specified data
 *          and x is the coordinate (double) where the function is evaluated.
 * farg   = pointer to user specified data (void *).
 * min_ns = minimum number of intermediate points to be used in smoothing (int).
 * max_ns = maximum number of intermediate points to be used in smoothing (int).
 * tol    = tolerance for weighing (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_adaptive_map(rgrid1d *grid, double (*func)(void *arg, double x), void *farg, int min_ns, int max_ns, double tol) {

  long i, nx = grid->nx, ns;
  double xc, step = grid->step;
  double tol2 = tol * tol;
  double sum, sump;
  double *value = grid->value;
  
  if (min_ns < 1) min_ns = 1;
  if (max_ns < min_ns) max_ns = min_ns;
  
#pragma omp parallel for firstprivate(farg,nx,min_ns,max_ns,step,func,value,tol2) private(i,ns,xc,sum,sump) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {      
    xc = (i - nx/2.0) * step;
    sum  = func(farg, xc);
    for(ns = min_ns; ns <= max_ns; ns *= 2) {
      sum  = linearly_weighted_integralr1d(func, farg, xc, step, ns);
      sump = linearly_weighted_integralr1d(func, farg, xc, step, ns+1);
      if(sqnorm(sum - sump) < tol2) break;
      value[i] = 0.5 * (sum + sump);
    }
  }
}

/*
 * Add two 1D grids ("gridc = grida + gridb").
 *
 * gridc = destination grid (rgrid1d *).
 * grida = 1st of the grids to be added (rgrid1d *).
 * gridb = 2nd of the grids to be added (rgrid1d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid1d_sum(rgrid1d *gridc, const rgrid1d *grida, const rgrid1d *gridb) {

  long i, nx = gridc->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    cvalue[i] = avalue[i] + bvalue[i];    
}

/* 
 * Subtract two grids ("gridc = grida - gridb").
 *
 * gridc = destination grid (rgrid1d *).
 * grida = 1st source grid (rgrid1d *).
 * gridb = 2nd source grid (rgrid1d *).
 *
 * No return value.
 *
 * Note: both source and destination may be the same.
 *
 */

EXPORT void rgrid1d_difference(rgrid1d *gridc, const rgrid1d *grida, const rgrid1d *gridb) {

  long i, nx = gridc->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = avalue[i] - bvalue[i];
}

/* 
 * Calculate product of two grids ("gridc = grida * gridb").
 *
 * gridc = destination grid (rgrid1d *).
 * grida = 1st source grid (rgrid1d *).
 * gridb = 2nd source grid (rgrid1d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid1d_product(rgrid1d *gridc, const rgrid1d *grida, const rgrid1d *gridb) {

  long i, nx = gridc->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = avalue[i] * bvalue[i];
}

/* 
 * Rise a grid to given power.
 *
 * gridb    = destination grid (rgrid1d *).
 * grida    = 1st source grid (rgrid1d *).
 * exponent = exponent to be used (double).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *       This routine uses pow() so that the exponent can be
 *       fractional but this is slow! Do not use this for integer
 *       exponents.
 *
 */

EXPORT void rgrid1d_power(rgrid1d *gridb, const rgrid1d *grida, double exponent) {

  long i, nx = gridb->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,exponent) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    bvalue[i] = pow(avalue[i], exponent);
}

/* 
 * Rise absolute value of a grid to given power.
 *
 * gridb    = destination grid (rgrid1d *).
 * grida    = 1st source grid (rgrid1d *).
 * exponent = exponent to be used (double).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *       This routine uses pow() so that the exponent can be
 *       fractional but this is slow! Do not use this for integer
 *       exponents.
 *
 */

EXPORT void rgrid1d_abs_power(rgrid1d *gridb, const rgrid1d *grida, double exponent) {

  long i, nx = gridb->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,exponent) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    bvalue[i] = pow(fabs(avalue[i]), exponent);
}

/*
 * Divide two grids ("gridc = grida / gridb").
 *
 * gridc = destination grid (rgrid1d *).
 * grida = 1st source grid (rgrid1d *).
 * gridb = 2nd source grid (rgrid1d *).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same.
 *
 */

EXPORT void rgrid1d_division(rgrid1d *gridc, const rgrid1d *grida, const rgrid1d *gridb) {

  long i, nx = gridc->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = avalue[i] / bvalue[i];
}

/* 
 * Conjugate product of two grids ("gridc = conj(grida) * gridb").
 *
 * gridc = destination grid (rgrid1d *).
 * grida = 1st source grid (rgrid1d *).
 * gridb = 2nd source grid (rgrid1d *).
 *
 * No return value.
 * 
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid1d_conjugate_product(rgrid1d *gridc, const rgrid1d *grida, const rgrid1d *gridb) {

  long i, nx = gridc->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = conj(avalue[i]) * bvalue[i];
}

/*
 * Add a constant to a 1D grid.
 *
 * grid = grid where the constant is added (rgrid1d *).
 * c    = constant to be added (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_add(rgrid1d *grid, double c) {

  long i, nx = grid->nx;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,c) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = value[i] + c;
}

/*
 * Multiply grid by a constant.
 *
 * grid = grid to be multiplied (rgrid1d *).
 * c    = multiplier (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_multiply(rgrid1d *grid, double c) {

  long i, nx = grid->nx;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,c) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = value[i] * c;
}

/* 
 * Add and multiply: grid = (grid + ca) * cm.
 *
 * grid = grid to be operated (rgrid1d *).
 * ca   = constant to be added (double).
 * cm   = multiplier (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_add_and_multiply(rgrid1d *grid, double ca, double cm) {

  long i, nx = grid->nx;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,ca,cm) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = (value[i] + ca) * cm;
}

/*
 * Multiply and add: grid = cm * grid + ca.
 *
 * grid = grid to be operated (rgrid1d *).
 * ca   = constant to be added (double).
 * cm   = multiplier (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_multiply_and_add(rgrid1d *grid, double cm, double ca) {

  long i, nx = grid->nx;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,ca,cm) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = value[i] * cm + ca;
}
 
/* 
 * Add scaled grid: gridc = gridc + d * grida
 *
 * gridc = destination grid for the operation (rgrid1d *).
 * d     = multiplier for the operation (double).
 * grida = source grid for the operation (rgrid1d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid1d_add_scaled(rgrid1d *gridc, double d, const rgrid1d *grida) {

  long i, nx = gridc->nx;
  double *avalue = grida->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(d,nx,avalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] += d * avalue[i];
}

/*
 * Add multiply two grids and a constant: gridc = gridc + c * grida * gridb.
 *
 * gridc = destination grid (rgrid1d *).
 * d     = constant multiplier (double).
 * grida = 1st source grid (rgrid1d *).
 * gridb = 2nd source grid (rgrid1d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid1d_add_scaled_product(rgrid1d *gridc, double d, const rgrid1d *grida, const rgrid1d *gridb) {

  long i, nx = gridc->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue,d) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] += d * avalue[i] * bvalue[i];
}

/*
 * Operate on a grid by a given operator: gridc = O(grida).
 *
 * gridc    = destination grid (rgrid1d *).
 * grida    = source grid (rgrid1d *).
 * operator = operator (double (*)(double)).
 *            (i.e., a function mapping a given C-number to another)
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid1d_operate_one(rgrid1d *gridc, const rgrid1d *grida, double (*operator)(double a)) {

  long i, nx = gridc->nx;
  double *avalue = grida->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,cvalue,operator) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = operator(avalue[i]);
}
 
/* 
 * Operate on two grids and place the result in third: gridc = O(grida, gridb).
 * where O is the operator.
 *
 * gridc    = destination grid (rgrid1d *).
 * grida    = 1s source grid (rgrid1d *).
 * gridb    = 2nd source grid (rgrid1d *).
 * operator = operator mapping grida and gridb (double (*)(double, double)).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid1d_operate_two(rgrid1d *gridc, const rgrid1d *grida, const rgrid1d *gridb, double (*operator)(double a, double b)) {

  long i, nx = gridc->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue,operator) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = operator(avalue[i], bvalue[i]);
}

/*
 * Operate on a grid by a given operator.
 *
 * grid     = grid to be operated (rgrid1d *).
 * operator = operator (void (*)(double *)).
 * 
 * No return value.
 *
 */

EXPORT void rgrid1d_transform_one(rgrid1d *grid, void (*operator)(double *a)) {

  long i, nx = grid->nx;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,operator) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    operator(&value[i]);
}

/*
 * Operate on two separate grids by a given operator.
 *
 * grida    = grid to be operated (rgrid1d *).
 * gridb    = grid to be operated (rgrid1d *).
 * operator = operator (void (*)(double *)).
 * 
 * No return value.
 *
 */

EXPORT void rgrid1d_transform_two(rgrid1d *grida, rgrid1d *gridb, void (*operator)(double *a, double *b)) {

  long i, nx = grida->nx;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,operator) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    operator(&avalue[i], &bvalue[i]);
}

/*
 * Integrate over a grid.
 *
 * grid = grid to be integrated (rgrid1d *).
 *
 * Returns the integral value (double).
 *
 */

EXPORT double rgrid1d_integral(const rgrid1d *grid) {

  long i, nx = grid->nx;
  double sum = 0.0;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value) private(i) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    sum += value[i];
  
  return sum * grid->step;
}

/*
 * Integrate over a grid with limits.
 *
 * grid = grid to be integrated (grid2d *).
 * xl   = lower limit for x (double).
 * xu   = upper limit for x (double).
 *
 * Returns the integral value (double).
 *
 */

EXPORT double rgrid1d_integral_region(const rgrid1d *grid, double xl, double xu) {

  long iu, il, i, nx = grid->nx;
  double *value = grid->value, sum;
  double step = grid->step;
   
  il = xl / step + nx/2;
  iu = xu / step + nx/2;
  
  sum = 0.0;
#pragma omp parallel for firstprivate(il,iu,nx,value) private(i) reduction(+:sum) default(none) schedule(runtime)
  for (i = il; i < iu; i++)
    sum += value[i];
  return sum * step;
}

/* 
 * Integrate over the grid squared (int grid^2).
 *
 * grid = grid to be integrated (rgrid1d *).
 *
 * Returns the integral (double ).
 *
 */

EXPORT double rgrid1d_integral_of_square(const rgrid1d *grid) {

  long i, nx = grid->nx;
  double sum = 0;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value) private(i) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
      sum += sqnorm(value[i]);
  
  return sum * grid->step;
}

/*
 * Calculate overlap between two grids (int grida^*gridb).
 *
 * grida = 1st grid (rgrid1d *).
 * gridb = 2nd grid (rgrid1d *).
 *
 * Returns the value of the overlap integral (double).
 *
 */

EXPORT double rgrid1d_integral_of_product(const rgrid1d *grida, const rgrid1d *gridb) {

  long i, nx = grida->nx;
  double sum = 0.0;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue) private(i) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    sum += conj(avalue[i]) * bvalue[i];
  
  return sum * grida->step;
}

/*
 * Calculate average value of a grid over the probability density given by the other.
 * (int grida gridb^2).
 *
 * grida = grid giving the probability (gridb^2) (rgrid1d *).
 * gridb = grid to be averaged (rgrid1d *).
 *
 * Returns the average value (double).
 *
 */

EXPORT double rgrid1d_grid_expectation_value(const rgrid1d *grida, const rgrid1d *gridb) {

  long i, nx = grida->nx;
  double sum = 0.0;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue) private(i) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    sum += sqnorm(avalue[i]) * bvalue[i];

  return sum * grida->step;
}

/*
 * Calculate the expectation value of a function over a grid.
 * (int grida^* func grida = int func grida^2).
 *
 * func  = function to be averaged (double (*)(void *, double, double)).
 *         The arguments are: optional arg, grida(x), x.
 * grida = grid giving the probability (grida^2) (rgrid1d *).
 *
 * Returns the average value (double).
 *
 */
 
EXPORT double rgrid1d_grid_expectation_value_func(void *arg, double (*func)(void *arg, double val, double x), const rgrid1d *grida) {
   
  long i, nx = grida->nx;
  double sum = 0.0;
  double *avalue = grida->value;
  double x, step = grida->step;
  
#pragma omp parallel for firstprivate(func,arg,nx,avalue,step) private(x,i) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    sum += sqnorm(avalue[i]) * func(arg, avalue[i], x);
  }
  
  return sum * step;
}

/* 
 * Integrate over the grid multiplied by weighting function (int grid w(x)).
 *
 * grid   = grid to be integrated over (rgrid1d *).
 * weight = function defining the weight (double (*)(double)).
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double rgrid1d_weighted_integral(const rgrid1d *grid, double (*weight)(void *farg, double x), void *farg) {

  long i, nx = grid->nx;
  double x, step = grid->step;
  double sum = 0.0;
  double w;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,step,value,weight,farg) private(w,i,x) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    w = weight(farg, x);
    sum += w * value[i];
  }

  return sum * grid->step;
}

/* 
 * Integrate over square of the grid multiplied by weighting function (int grid^2 w(x)).
 *
 * grid   = grid to be integrated over (rgrid1d *).
 * weight = function defining the weight (double (*)(double)).
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double rgrid1d_weighted_integral_of_square(const rgrid1d *grid, double (*weight)(void *farg, double x), void *farg) {

  long i, nx = grid->nx;
  double x, step = grid->step;
  double sum = 0, w;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,step,value,weight,farg) private(w,i,x) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;    
    w = weight(farg, x);
    sum += w * sqnorm(value[i]);
  }
  
  return sum * grid->step;
}

/* 
 * Differentiate a grid with respect to x.
 *
 * grid     = grid to be differentiated (rgrid1d *).
 * gradient = differentiated grid output (rgrid1d *).
 * 
 * No return value.
 *
 */

EXPORT void rgrid1d_fd_gradient(const rgrid1d *grid, rgrid1d *gradient) {

  rgrid1d_fd_gradient_x(grid, gradient);
}

EXPORT void rgrid1d_fd_gradient_x(const rgrid1d *grid, rgrid1d *gradient) {

  long i, nx = grid->nx;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double *lvalue = gradient->value;

  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in rgrid1d_fd_gradient_x().\n");
    return;
  }

#pragma omp parallel for firstprivate(nx,lvalue,inv_delta,grid) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    lvalue[i] = inv_delta * (rgrid1d_value_at_index(grid, i+1) - rgrid1d_value_at_index(grid, i-1));
}

/*
 * Calculate second derivative of a grid with respect to x.
 *
 * grid    = grid to be differentiated twice (rgrid1d *).
 * laplace = output grid for the operation (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_fd_laplace(const rgrid1d *grid, rgrid1d *laplace) {

  long i;
  long nx = grid->nx;
  double inv_delta2 = 1./ (grid->step * grid->step);
  double *lvalue = laplace->value;
  
#pragma omp parallel for firstprivate(nx,lvalue,inv_delta2,grid) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    lvalue[i] = inv_delta2 * (-2.0 * rgrid1d_value_at_index(grid, i) + rgrid1d_value_at_index(grid, i-1) + rgrid1d_value_at_index(grid, i+1));
}

/*
 * Calculate dot product of the gradient of the grid.
 *
 * grid          = source grid for the operation (rgrid1d *).
 * grad_dot_grad = destination grid (rgrid1d *).
 *
 * No return value.
 *
 * Note: grid and grad_dot_grad may not be the same grid.
 *
 */

EXPORT void rgrid1d_fd_gradient_dot_gradient(const rgrid1d *grid, rgrid1d *grad_dot_grad){

  long i;
  long nx = grid->nx;
  double inv_2delta2 = 1./ (2.*grid->step * 2.*grid->step);
  double *gvalue = grad_dot_grad->value;
  
  /*  grad f(x,y,z) dot grad f(x,y,z) = [ |f(+,0,0) - f(-,0,0)|^2 + |f(0,+,0) - f(0,-,0)|^2 + |f(0,0,+) - f(0,0,-)|^2 ] / (2h)^2 */
#pragma omp parallel for firstprivate(nx,gvalue,inv_2delta2,grid) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    gvalue[i] = inv_2delta2 * sqnorm(rgrid1d_value_at_index(grid,i+1) - rgrid1d_value_at_index(grid,i-1));
}

/*
 * Print the grid with both real and imaginary parts into file (ASCII format).
 *
 * grid = grid to be printed out (rgrid1d *).
 * out  = output file pointer (FILE *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_print(const rgrid1d *grid, FILE *out) {

  long i;

  for(i = 0; i < grid->nx; i++)
    fprintf(out, "%16.8le\n", grid->value[i]);
}

/*
 * Perform Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (rgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       Also no normalization is performed.
 *
 */

EXPORT void rgrid1d_fft(rgrid1d *grid) {

  if (!grid->plan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in rgrid1d_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    rgrid1d_fftw_alloc(grid);
  }
  rgrid1d_fftw(grid);
}

/*
 * Perform normalized Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (rgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void rgrid1d_fourier_transform(rgrid1d *grid) {

  rgrid1d_fft(grid);
  cgrid1d_multiply(grid->cint, grid->step);
}

/*
 * Perform inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       No normalization.
 *
 */

EXPORT void rgrid1d_inverse_fft(rgrid1d *grid) {

  if (!grid->iplan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in rgrid1d_inverse_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    rgrid1d_fftw_alloc(grid);
  }
  rgrid1d_fftw_inv(grid);
}

/*
 * Perform scaled inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid1d *).
 * c    = scaling factor (i.e. the output is multiplied by this constant) (double).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void rgrid1d_scaled_inverse_fft(rgrid1d *grid, double c) {

  rgrid1d_inverse_fft(grid);
  rgrid1d_multiply(grid, c);  
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by 1 / N.
 * (N = number of grid points)
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void rgrid1d_inverse_fft_norm(rgrid1d *grid) {

  rgrid1d_scaled_inverse_fft(grid, 1.0 / grid->nx);
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by 1 / (N step).
 * (N = number of grid points, step = grid step length)
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       This is the inverse routine to be used with rgrid1d_fourier_transform().
 *
 */

EXPORT void rgrid1d_inverse_fourier_transform(rgrid1d *grid) {

  rgrid1d_scaled_inverse_fft(grid, 1.0 / (grid->nx * grid->step));
}

/*
 * Convolute FFT transformed grids. To apply this on two grids (grida and gridb)
 * and place the result in gridc:
 * rgrid1d_fft(grida);
 * rgrid1d_fft(gridb);
 * rgrid1d_convolue(gridc, grida, gridb);
 * rgrid1d_inverse_fft(gridc);
 * gridc now contains the convolution of grida and gridb.
 *
 * grida = 1st grid to be convoluted (rgrid1d *).
 * gridb = 2nd grid to be convoluted (rgrid1d *).
 * gridc = output (rgrid1d *).
 *
 * No return value.
 *
 * Note: the input/output grids may be the same.
 *
 */

EXPORT void rgrid1d_fft_convolute(rgrid1d *gridc, const rgrid1d *grida, const rgrid1d *gridb) {

  long i, nx;
  double step, norm;
  double complex *avalue, *bvalue, *cvalue;
  
  /* int f(r) g(r-r') dr' = iF[ F[f] F[g] ] = (step / N) iFFT[ FFT[f] FFT[g] ] */
  nx = gridc->cint->nx;
  step = gridc->cint->step;
  
  avalue = grida->cint->value;
  bvalue = gridb->cint->value;
  cvalue = gridc->cint->value;
  
  norm = step / (double) gridc->nx;
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue,norm) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    if (i & 1)
      cvalue[i] = -norm * avalue[i] * bvalue[i];
    else
      cvalue[i] = norm * avalue[i] * bvalue[i];
}

/* Boundary condition routines */

EXPORT double rgrid1d_value_outside_constantdirichlet(const rgrid1d *grid, long i) {

#if __DEBUG__
  if (i >= 0 && i < grid->nx) {
    fprintf(stderr, "libgrid: Error in rgrid1d_value_outside_constantdirichlet(). "
	    "Index (%ld) is not outside the grid.", i);
    abort();
  }
#endif
  
  return  *((double *) grid->outside_params_ptr);
}

EXPORT double rgrid1d_value_outside_neumann(const rgrid1d *grid, long i) {

  long nx = grid->nx;

#if __DEBUG__
  if (i >= 0 && i < grid->nx) {
    fprintf(stderr, "libgrid: Error in rgrid1d_value_outside_neumann(). "
	    "Index (%ld) is not outside the grid.", i);
    abort();
  }
#endif

  if(i < 0) i = 1;
  if(i > nx-1) i = nx - 2;

  return grid->value[i];
  
  /* This is the old mirror Neumann code */
#if 0
  nd = nx * 2;
  if (i < 0)   i  = -i;
  if (i >= nd) i %= nd;
  if (i >= nx) i  = nd - i;
  
  return grid->value[i];
#endif
}

EXPORT double rgrid1d_value_outside_periodic(const rgrid1d *grid, long i) {

  long nx = grid->nx;
  
#if __DEBUG__
  if (i >= 0 && i < grid->nx) {
    fprintf(stderr, "libgrid: Error in rgrid1d_value_outside_periodic(). "
	    "Index (%ld) is not outside the grid.", i);
    abort();
  }
#endif
  
  i %= nx;
  if (i < 0) i = nx + i;
  
  return grid->value[i];  
}

/* End boundary condition routines */

/*
 * Access grid point at given index.
 *
 * grid = grid to be accessed (rgrid1d *).
 * i    = index (long).
 *
 * Returns grid value at index i.
 *
 */

EXPORT inline double rgrid1d_value_at_index(const rgrid1d *grid, long i) {

  if (i < 0 || i >= grid->nx)
    return grid->value_outside(grid, i);
  return grid->value[i];
}

/*
 * Access grid point at given x value (with linear interpolation).
 *
 * grid = grid to be accessed (rgrid1d *).
 * x    = x value (double).
 *
 * Returns grid value at x.
 *
 */

EXPORT inline double rgrid1d_value(const rgrid1d *grid, double x) {

  double f0, f1;
  double omx, step = grid->step;
  long i;
  
  /* i to index and 0 <= x < 1 */
  i = (x /= step);
  if (x < 0) i--;
  x -= i;
  i += grid->nx / 2;

  /* linear extrapolation 
   *
   * f(x) = (1-x) f(0) + x f(1) */ 
  f0 = rgrid1d_value_at_index(grid, i);
  f1 = rgrid1d_value_at_index(grid, i+1);
  
  omx = 1 - x;

  return omx * f0 + x * f1;
}

/*
 * Copy a real grid to a complex grid (to real part).
 *
 * dest   = destination grid (cgrid1d *).
 * source = source grid (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_real_to_complex_re(cgrid1d *dest, rgrid1d *source) {
  
  long i, nx = source->nx;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nx,dst,src) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    dst[i] = (double complex) src[i];
}

/*
 * Copy a real grid to a complex grid (to imaginary part).
 *
 * dest   = destination grid (cgrid1d *).
 * source = source grid (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_real_to_complex_im(cgrid1d *dest, rgrid1d *source) {
  
  long i, nx = source->nx;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nx,dst,src) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    dst[i] = I * (double complex) src[i];
}

/*
 * Product of a real grid with a complex grid.
 *
 * dest   = destination grid (cgrid1d *).
 * source = source grid (rgrid1d *).
 *
 * "dest(complex) = dest(complex) * source(real)"
 *
 * No return value.
 *
 */

EXPORT void grid1d_product_complex_with_real(cgrid1d *dest, rgrid1d *source) {
  
  long i, nx = source->nx;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nx,dst,src) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    dst[i] *= (double complex) src[i];
}

/*
 * Add a real grid to a complex grid (to real part).
 *
 * dest   = destination grid (cgrid1d *).
 * source = source grid (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_add_real_to_complex_re(cgrid1d *dest, rgrid1d *source) {
  
  long i, nx = source->nx;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nx,dst,src) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    dst[i] += (double complex) src[i];
}

/*
 * Add a real grid to a complex grid (to imaginary part).
 *
 * dest   = destination grid (cgrid1d *).
 * source = source grid (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_add_real_to_complex_im(cgrid1d *dest, rgrid1d *source) {
  
  long i, nx = source->nx;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nx,dst,src) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    dst[i] += I * (double complex) src[i];
}

/*
 * Copy imaginary part of a complex grid to a real grid.
 *
 * dest   = destination grid (rgrid1d *).
 * source = source grid (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_complex_im_to_real(rgrid1d *dest, cgrid1d *source) {
  
  long i, nx = source->nx;
  double complex *src = source->value;
  double *dst = dest->value;
  
  dest->nx = source->nx;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nx,dst,src) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    dst[i] = cimag(src[i]);
}

/*
 * Copy real part of a complex grid to a real grid.
 *
 * dest   = destination grid (rgrid1d *).
 * source = source grid (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_complex_re_to_real(rgrid1d *dest, cgrid1d *source) {
  
  long i, nx = source->nx;
  double complex *src = source->value;
  double *dst = dest->value;
  
  dest->nx = source->nx;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nx,dst,src) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    dst[i] = creal(src[i]);
}

/*
 * Extrapolate between two different grid sizes.
 *
 * dest = Destination grid (rgrid1d *; output).
 * src  = Source grid (rgrid1d *; input).
 *
 */

EXPORT void rgrid1d_extrapolate(rgrid1d *dest, rgrid1d *src) {

  long i, nx = dest->nx;
  double step = dest->step, x;

  for (i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    dest->value[i] = rgrid1d_value(src, x);
  }
}


/*
 * Get the largest value contained in a grid
 */
EXPORT double rgrid1d_max(rgrid1d *grid){
	long nx = grid->nx ;
	double *val = grid->value , max_val = val[0] ;
	long i ;

#pragma omp parallel for firstprivate(nx, val) private(i) reduction(max: max_val) default(none) schedule(runtime)
	for(i=0 ; i< nx ; i++){
		if(val[i] > max_val)
			max_val = val[i] ;
	}
	return max_val ;
}

/*
 * Get the lowest value contained in a grid
 */
EXPORT double rgrid1d_min(rgrid1d *grid){
	long nx = grid->nx ;
	double *val = grid->value , min_val = val[0] ;
	long i ;

#pragma omp parallel for firstprivate(nx, val) private(i) reduction(min: min_val) default(none) schedule(runtime)
	for(i=0 ; i< nx ; i++){
		if(val[i] < min_val)
			min_val = val[i] ;
	}
	return min_val ;
}
