/*
 * Routines for 1D complex grids.
 *
 * To debug memory allocation, define __DEBUG__.
 *
 * Periodic grid: the most positive value is not stored (rolled over to -).
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
 *                      CGRID1D_DIRICHLET_BOUNDARY: Dirichlet boundary.
 *                      or CGRID1D_NEUMANN_BOUNDARY: Neumann boundary.
 *                      or CGRID1D_PERIODIC_BOUNDARY: Periodic boundary.
 *                      or user supplied function with pointer to grid and
 *                         grid index as parameters to provide boundary access.
 * outside_params_ptr = pointer for passing parameters for the given boundary
 *                      access function. Use 0 to with the predefined boundary
 *                      functions (CGRID1D_* above).
 *
 * Return value: pointer to the allocated grid (cgrid1d *). Returns NULL on
 * error.
 *
 */

EXPORT cgrid1d *cgrid1d_alloc(long nx, double step, double complex (*value_outside)(const cgrid1d *grid, long i), void *outside_params_ptr) {

  cgrid1d *grid;
  
  grid = (cgrid1d *) malloc(sizeof(cgrid1d));
  if (!grid) {
    fprintf(stderr, "libgrid: Error in cgrid1d_alloc(). Could not allocate memory for 1d grid.\n");
    return 0;
  }
  
  if(!(grid->value = (double complex *) fftw_malloc(nx * sizeof(double complex)))) {
    fprintf(stderr, "libgrid: Error in cgrid1d_alloc(). Could not allocate memory for cgrid1d->value.\n");
    free(grid);
    return 0;
  }
  
  grid->nx = nx;
  grid->step = step;
  
  if (value_outside)
    grid->value_outside = value_outside;
  else
    grid->value_outside = cgrid1d_value_outside_constantdirichlet;
    
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

  cgrid1d_constant(grid, NAN);
  
  return grid;
}

/*
 * Free 1D grid.
 *
 * grid = pointer to 1D grid to be freed (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_free(cgrid1d *grid) {

  if (grid) {
    if(grid->value) fftw_free(grid->value);
    cgrid1d_fftw_free(grid);
    free(grid);
  }
}

/* 
 * Write 1D grid on disk in binary format.
 *
 * grid = 1D grid to be written (cgrid1d *).
 * out  = file handle for the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_write(cgrid1d *grid, FILE *out) {

  fwrite(&grid->nx, sizeof(long), 1, out);
  fwrite(&grid->step, sizeof(double), 1, out);
  fwrite(grid->value, sizeof(double complex), grid->nx, out);
}

/* 
 * Read 1D grid from disk in binary format.
 *
 * grid = 1D grid to be read (cgrid1d *).
 * in   = file handle for reading the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_read(cgrid1d *grid, FILE *in) {

  long nx;
  double step;
  
  fread(&nx, sizeof(long), 1, in);
  fread(&step, sizeof(double), 1, in);
  
  if (nx != grid->nx || step != grid->step) {
    cgrid1d *tmp;

    fprintf(stderr, "libgrid: Grid in file has different size than grid in memory.\n");
    fprintf(stderr, "libgrid: Interpolating between grids.\n");
    if(!(tmp = cgrid1d_alloc(nx, step, grid->value_outside, NULL))) {
      fprintf(stderr, "libgrid: Error allocating grid in cgrid1d_read().\n");
      abort();
    }
    fread(tmp->value, sizeof(double complex), nx, in);
    cgrid1d_extrapolate(grid, tmp);
    cgrid1d_free(tmp);
    return;
  }
  
  fread(grid->value, sizeof(double complex), grid->nx, in);
}

/*
 * Copy 1D grid from one grid to another.
 *
 * copy = destination grid (cgrid1d *).
 * grid = source grid (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_copy(cgrid1d *copy, const cgrid1d *grid) {

  long i, nx = grid->nx;
  double complex *gvalue = grid->value;
  double complex *cvalue = copy->value;
  
  copy->nx = grid->nx;
  copy->step = grid->step;
  
#pragma omp parallel for firstprivate(nx,gvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = gvalue[i];
}

/*
 * Take complex conjugate of 1D grid.
 * 
 * conjugate = destination for complex conjugated grid (cgrid1d *).
 * grid      = source grid for the operation (cgrid1d *).
 *
 * No return value.
 *
 * Note: source and destination may be the same grid.
 * 
 */

EXPORT void cgrid1d_conjugate(cgrid1d *conjugate, const cgrid1d *grid) {

  long i, nx = grid->nx;
  double complex *cvalue = conjugate->value;
  double complex *gvalue = grid->value;
  
  conjugate->nx = nx;
  conjugate->step = grid->step;
  
#pragma omp parallel for firstprivate(cvalue,gvalue,nx) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    cvalue[i] = conj(gvalue[i]);
}

/*
 * Shift 1D grid by given amount spatially.
 *
 * shifted = destination grid for the operation (cgrid1d *).
 * grid    = source grid for the operation (cgrid1d *).
 * x       = shift spatially by this amount (double).
 *
 * No return value.
 *
 * NOTE: Source and destination may be the same grid.
 *
 */

EXPORT void cgrid1d_shift(cgrid1d *shifted, const cgrid1d *grid, double x) {

  sShiftParametersc1d params;

  /* shift by (x) i.e. current grid center to (x) */
  params.x = x; 
  params.grid = grid;
  cgrid1d_map(shifted, shift_cgrid1d, &params);
}

/* 
 * Zero 1D grid.
 *
 * grid = grid to be zeroed (cgrid1d *).
 *
 * No return value.
 * 
 */

EXPORT void cgrid1d_zero(cgrid1d *grid) { 

  cgrid1d_constant(grid, 0.0); 
}

/* 
 * Set 1D grid to a constant value.
 *
 * grid = grid to be set (cgrid1d *).
 * c    = value (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_constant(cgrid1d *grid, double complex c) {

  long i, nx = grid->nx;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(value,c,nx) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = c;
}

/*
 * Multiply a grid by a given function.
 *
 * grid = destination grid for the operation (cgrid1d *).
 * func = function providing the mapping (double complex (*)(void *, double)).
 *        The first argument (void *) is for external user specified data
 *        and x is the coordinate (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

 EXPORT void cgrid1d_product_func(cgrid1d *grid, double complex (*func)(void *arg, double x), void *farg) {

  long i, nx = grid->nx;
  double x, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,step,func,value) private(i,x) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    value[i] *= func(farg, x);
  }
}

/*
 * Map a given function onto 1D grid.
 *
 * grid = destination grid for the operation (cgrid1d *).
 * func = function providing the mapping (double complex (*)(void *, double)).
 *        The first argument (void *) is for external user specified data
 *        and x is the coordinate (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_map(cgrid1d *grid, double complex (*func)(void *arg, double x), void *farg) {

  long i, nx = grid->nx;
  double x, step = grid->step;
  double complex *value = grid->value;
  
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
 * grid = destination grid for the operation (cgrid1d *).
 * func = function providing the mapping (double complex (*)(void *, double)).
 *        The first argument (void *) is for external user specified data
 *        and x is the coordinate (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 * ns   = number of intermediate points to be used in smoothing (int).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_smooth_map(cgrid1d *grid, double complex (*func)(void *arg, double x), void *farg, int ns) {

  long i, nx = grid->nx;
  double xc, step = grid->step;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ns,step,func,value) private(i,xc) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {
    xc = (i - nx/2.0) * step;
    value[i] = linearly_weighted_integralc1d(func, farg, xc, step, ns);
  }
}

/*
 * Map a given function onto 1D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid. Limits for intermediate steps and
 * tolerance can be given.
 *
 * grid   = destination grid for the operation (cgrid1d *).
 * func   = function providing the mapping (double complex (*)(void *, double)).
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

EXPORT void cgrid1d_adaptive_map(cgrid1d *grid, double complex (*func)(void *arg, double x), void *farg, int min_ns, int max_ns, double tol) {

  long i, nx = grid->nx, ns;
  double xc, step = grid->step;
  double tol2 = tol * tol;
  double complex sum, sump;
  double complex *value = grid->value;
  
  if (min_ns < 1) min_ns = 1;
  if (max_ns < min_ns) max_ns = min_ns;
  
#pragma omp parallel for firstprivate(farg,nx,min_ns,max_ns,step,func,value,tol2) private(i,ns,xc,sum,sump) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {      
    xc = (i - nx/2.0) * step;
    sum  = func(farg, xc);
    for(ns = min_ns; ns <= max_ns; ns *= 2) {
      sum  = linearly_weighted_integralc1d(func, farg, xc, step, ns);
      sump = linearly_weighted_integralc1d(func, farg, xc, step, ns+1);
      if(sqnorm(sum - sump) < tol2) break;
      value[i] = 0.5 * (sum + sump);
    }
  }
}

/*
 * Add two 1D grids ("gridc = grida + gridb").
 *
 * gridc = destination grid (cgrid1d *).
 * grida = 1st of the grids to be added (cgrid1d *).
 * gridb = 2nd of the grids to be added (cgrid1d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid1d_sum(cgrid1d *gridc, const cgrid1d *grida, const cgrid1d *gridb) {

  long i, nx = gridc->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    cvalue[i] = avalue[i] + bvalue[i];    
}

/* 
 * Subtract two grids ("gridc = grida - gridb").
 *
 * gridc = destination grid (cgrid1d *).
 * grida = 1st source grid (cgrid1d *).
 * gridb = 2nd source grid (cgrid1d *).
 *
 * No return value.
 *
 * Note: both source and destination may be the same.
 *
 */

EXPORT void cgrid1d_difference(cgrid1d *gridc, const cgrid1d *grida, const cgrid1d *gridb) {

  long i, nx = gridc->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = avalue[i] - bvalue[i];
}

/* 
 * Calculate product of two grids ("gridc = grida * gridb").
 *
 * gridc = destination grid (cgrid1d *).
 * grida = 1st source grid (cgrid1d *).
 * gridb = 2nd source grid (cgrid1d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid1d_product(cgrid1d *gridc, const cgrid1d *grida, const cgrid1d *gridb) {

  long i, nx = gridc->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = avalue[i] * bvalue[i];
}

/* 
 * Rise a grid to given power.
 *
 * gridb    = destination grid (cgrid1d *).
 * grida    = 1st source grid (cgrid1d *).
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

EXPORT void cgrid1d_power(cgrid1d *gridb, const cgrid1d *grida, double exponent) {

  long i, nx = gridb->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,exponent) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    bvalue[i] = pow(avalue[i], exponent);
}

/* 
 * Rise absolute value of a grid to given power.
 *
 * gridb    = destination grid (cgrid1d *).
 * grida    = 1st source grid (cgrid1d *).
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

EXPORT void cgrid1d_abs_power(cgrid1d *gridb, const cgrid1d *grida, double exponent) {

  long i, nx = gridb->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,exponent) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    bvalue[i] = pow(cabs(avalue[i]), exponent);
}

/*
 * Divide two grids ("gridc = grida / gridb").
 *
 * gridc = destination grid (cgrid1d *).
 * grida = 1st source grid (cgrid1d *).
 * gridb = 2nd source grid (cgrid1d *).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same.
 *
 */

EXPORT void cgrid1d_division(cgrid1d *gridc, const cgrid1d *grida, const cgrid1d *gridb) {

  long i, nx = gridc->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = avalue[i] / bvalue[i];
}

/* 
 * Conjugate product of two grids ("gridc = conj(grida) * gridb").
 *
 * gridc = destination grid (cgrid1d *).
 * grida = 1st source grid (cgrid1d *).
 * gridb = 2nd source grid (cgrid1d *).
 *
 * No return value.
 * 
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid1d_conjugate_product(cgrid1d *gridc, const cgrid1d *grida, const cgrid1d *gridb) {

  long i, nx = gridc->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = conj(avalue[i]) * bvalue[i];
}

/*
 * Add a constant to a 1D grid.
 *
 * grid = grid where the constant is added (cgrid1d *).
 * c    = constant to be added (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_add(cgrid1d *grid, double complex c) {

  long i, nx = grid->nx;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,c) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = value[i] + c;
}

/*
 * Multiply grid by a constant.
 *
 * grid = grid to be multiplied (cgrid1d *).
 * c    = multiplier (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_multiply(cgrid1d *grid, double complex c) {

  long i, nx = grid->nx;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,c) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = value[i] * c;
}

/* 
 * Add and multiply: grid = (grid + ca) * cm.
 *
 * grid = grid to be operated (cgrid1d *).
 * ca   = constant to be added (double complex).
 * cm   = multiplier (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_add_and_multiply(cgrid1d *grid, double complex ca, double complex cm) {

  long i, nx = grid->nx;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,ca,cm) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = (value[i] + ca) * cm;
}

/*
 * Multiply and add: grid = cm * grid + ca.
 *
 * grid = grid to be operated (cgrid1d *).
 * ca   = constant to be added (double complex).
 * cm   = multiplier (double complex).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_multiply_and_add(cgrid1d *grid, double complex cm, double complex ca) {

  long i, nx = grid->nx;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,ca,cm) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    value[i] = value[i] * cm + ca;
}
 
/* 
 * Add scaled grid: gridc = gridc + d * grida
 *
 * gridc = destination grid for the operation (cgrid1d *).
 * d     = multiplier for the operation (double complex).
 * grida = source grid for the operation (cgrid1d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid1d_add_scaled(cgrid1d *gridc, double complex d, const cgrid1d *grida) {

  long i, nx = gridc->nx;
  double complex *avalue = grida->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(d,nx,avalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] += d * avalue[i];
}

/*
 * Add multiply two grids and a constant: gridc = gridc + c * grida * gridb.
 *
 * gridc = destination grid (cgrid1d *).
 * d     = constant multiplier (double complex).
 * grida = 1st source grid (cgrid1d *).
 * gridb = 2nd source grid (cgrid1d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid1d_add_scaled_product(cgrid1d *gridc, double complex d, const cgrid1d *grida, const cgrid1d *gridb) {

  long i, nx = gridc->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue,d) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] += d * avalue[i] * bvalue[i];
}

/*
 * Operate on a grid by a given operator: gridc = O(grida).
 *
 * gridc    = destination grid (cgrid1d *).
 * grida    = source grid (cgrid1d *).
 * operator = operator (double complex (*)(double complex)).
 *            (i.e. a function mapping a given C-number to another)
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid1d_operate_one(cgrid1d *gridc, const cgrid1d *grida, double complex (*operator)(double complex a)) {

  long i, nx = gridc->nx;
  double complex *avalue = grida->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,cvalue,operator) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = operator(avalue[i]);
}
 
/* 
 * Operate on two grids and place the result in third: gridc = O(grida, gridb).
 * where O is the operator.
 *
 * gridc    = destination grid (cgrid1d *).
 * grida    = 1s source grid (cgrid1d *).
 * gridb    = 2nd source grid (cgrid1d *).
 * operator = operator mapping grida and gridb (double complex (*)(double
 *            complex, double complex)).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void cgrid1d_operate_two(cgrid1d *gridc, const cgrid1d *grida, const cgrid1d *gridb, double complex (*operator)(double complex a, double complex b)) {

  long i, nx = gridc->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  double complex *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,cvalue,operator) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = operator(avalue[i], bvalue[i]);
}

/*
 * Operate on a grid by a given operator.
 *
 * grid     = grid to be operated (cgrid1d *).
 * operator = operator (void (*)(double complex *)).
 * 
 * No return value.
 *
 */

EXPORT void cgrid1d_transform_one(cgrid1d *grid, void (*operator)(double complex *a)) {

  long i, nx = grid->nx;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value,operator) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    operator(&value[i]);
}

/*
 * Operate on two separate grids by a given operator.
 *
 * grida    = grid to be operated (cgrid1d *).
 * gridb    = grid to be operated (cgrid1d *).
 * operator = operator (void (*)(double complex *)).
 * 
 * No return value.
 *
 */

EXPORT void cgrid1d_transform_two(cgrid1d *grida, cgrid1d *gridb, void (*operator)(double complex *a, double complex *b)) {

  long i, nx = grida->nx;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue,operator) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    operator(&avalue[i], &bvalue[i]);
}

/*
 * Integrate over a grid.
 *
 * grid = grid to be integrated (cgrid1d *).
 *
 * Returns the integral value (double complex).
 *
 */

EXPORT double complex cgrid1d_integral(const cgrid1d *grid) {

  long i, nx = grid->nx;
  double complex sum = 0.0;
  double complex *value = grid->value;
  
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
 * Returns the integral value (double complex).
 *
 */

EXPORT double complex cgrid1d_integral_region(const cgrid1d *grid, double xl, double xu) {

  long iu, il, i, nx = grid->nx;
  double complex *value = grid->value, sum;
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
 * Integrate over the grid squared (int |grid|^2).
 *
 * grid = grid to be integrated (cgrid1d *).
 *
 * Returns the integral (double complex).
 *
 */

EXPORT double cgrid1d_integral_of_square(const cgrid1d *grid) {

  long i, nx = grid->nx;
  double sum = 0;
  double complex *value = grid->value;
  
#pragma omp parallel for firstprivate(nx,value) private(i) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
      sum += sqnorm(value[i]);
  
  return sum * grid->step;
}

/*
 * Calculate overlap between two grids (int grida^*gridb).
 *
 * grida = 1st grid (complex conjugated) (cgrid1d *).
 * gridb = 2nd grid (no complex conjugation) (cgrid1d *).
 *
 * Returns the value of the overlap integral (double complex).
 *
 */

EXPORT double complex cgrid1d_integral_of_conjugate_product(const cgrid1d *grida, const cgrid1d *gridb) {

  long i, nx = grida->nx;
  double complex sum = 0.0;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue) private(i) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    sum += conj(avalue[i]) * bvalue[i];
  
  return sum * grida->step;
}

/*
 * Calculate average value of a grid over the probability density given by the other.
 * (int grida |gridb|^2).
 *
 * grida = grid giving the probability (|gridb|^2) (cgrid1d *).
 * gridb = grid to be averaged (cgrid1d *).
 *
 * Returns the average value (double complex).
 *
 */

EXPORT double complex cgrid1d_grid_expectation_value(const cgrid1d *grida, const cgrid1d *gridb) {

  long i, nx = grida->nx;
  double complex sum = 0.0;
  double complex *avalue = grida->value;
  double complex *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nx,avalue,bvalue) private(i) reduction(+:sum) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    sum += sqnorm(avalue[i]) * bvalue[i];

  return sum * grida->step;
}
/*
 * Calculate the expectation value of a function over a grid.
 * (int grida^* func grida = int func |grida|^2).
 *
 * func  = function to be averaged (double complex (*)(void *, double complex, double)).
 *         The arguments are: optional arg, grida(x), x.
 * grida = grid giving the probability (|grida|^2) (cgrid1d *).
 *
 * Returns the average value (double complex).
 *
 */
 
EXPORT double complex cgrid1d_grid_expectation_value_func(void *arg, double complex (*func)(void *arg, double complex val, double x), const cgrid1d *grida) {
   
  long i, nx = grida->nx;
  double complex sum = 0.0;
  double complex *avalue = grida->value;
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
 * grid   = grid to be integrated over (cgrid1d *).
 * weight = function defining the weight (double complex (*)(double)).
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double complex).
 *
 */

EXPORT double complex cgrid1d_weighted_integral(const cgrid1d *grid, double complex (*weight)(void *farg, double x), void *farg) {

  long i, nx = grid->nx;
  double x, step = grid->step;
  double complex sum = 0.0;
  double complex w;
  double complex *value = grid->value;
  
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
 * grid   = grid to be integrated over (cgrid1d *).
 * weight = function defining the weight (double complex (*)(double)).
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double complex).
 *
 */

EXPORT double cgrid1d_weighted_integral_of_square(const cgrid1d *grid, double (*weight)(void *farg, double x), void *farg) {

  long i, nx = grid->nx;
  double x, step = grid->step;
  double sum = 0, w;
  double complex *value = grid->value;
  
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
 * grid     = grid to be differentiated (cgrid1d *).
 * gradient = differentiated grid output (cgrid1d *).
 * 
 * No return value.
 *
 */

EXPORT void cgrid1d_fd_gradient(const cgrid1d *grid, cgrid1d *gradient) {

  cgrid1d_fd_gradient_x(grid, gradient);
}

EXPORT void cgrid1d_fd_gradient_x(const cgrid1d *grid, cgrid1d *gradient) {

  long i;
  long nx = grid->nx;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double complex *lvalue = gradient->value;

  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in cgrid1d_fd_gradient_x().\n");
    return;
  }
#pragma omp parallel for firstprivate(nx,lvalue,inv_delta,grid) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    lvalue[i] = inv_delta * (cgrid1d_value_at_index(grid, i+1) - cgrid1d_value_at_index(grid, i-1));
}

/*
 * Calculate second derivative of a grid with respect to x.
 *
 * grid    = grid to be differentiated twice (cgrid1d *).
 * laplace = output grid for the operation (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_fd_laplace(const cgrid1d *grid, cgrid1d *laplace) {

  long i;
  long nx = grid->nx;
  double inv_delta2 = 1./ (grid->step * grid->step);
  double complex *lvalue = laplace->value;
  
#pragma omp parallel for firstprivate(nx,lvalue,inv_delta2,grid) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    lvalue[i] = inv_delta2 * (-2.0 * cgrid1d_value_at_index(grid, i) + cgrid1d_value_at_index(grid, i-1) + cgrid1d_value_at_index(grid, i+1));
}

/*
 * Calculate dot product of the gradient of the grid.
 *
 * grid          = source grid for the operation (cgrid1d *).
 * grad_dot_grad = destination grid (cgrid1d *).
 *
 * No return value.
 *
 * Note: grid and grad_dot_grad may not be the same grid.
 *
 */

EXPORT void cgrid1d_fd_gradient_dot_gradient(const cgrid1d *grid, cgrid1d *grad_dot_grad){

  long i;
  long nx = grid->nx;
  double inv_2delta2 = 1./ (2.*grid->step * 2.*grid->step);
  double complex *gvalue = grad_dot_grad->value;
  
  /*  grad f(x,y,z) dot grad f(x,y,z) = [ |f(+,0,0) - f(-,0,0)|^2 + |f(0,+,0) - f(0,-,0)|^2 + |f(0,0,+) - f(0,0,-)|^2 ] / (2h)^2 */
#pragma omp parallel for firstprivate(nx,gvalue,inv_2delta2,grid) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    gvalue[i] = inv_2delta2 * sqnorm(cgrid1d_value_at_index(grid,i+1) - cgrid1d_value_at_index(grid,i-1));
}

/*
 * Print the grid with both real and imaginary parts into file (ASCII format).
 *
 * grid = grid to be printed out (cgrid1d *).
 * out  = output file pointer (FILE *).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_print(const cgrid1d *grid, FILE *out) {

  long i;

  for(i = 0; i < grid->nx; i++)
    fprintf(out, "%16.8le %16.8le\n", creal(grid->value[i]), cimag(grid->value[i]));
}

/*
 * Perform Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (cgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       Also no normalization is performed.
 *
 */

EXPORT void cgrid1d_fft(cgrid1d *grid) {

  if (!grid->plan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in cgrid1d_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    cgrid1d_fftw_alloc(grid);
  }
  cgrid1d_fftw(grid);
}

/*
 * Perform normalized Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (cgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void cgrid1d_fourier_transform(cgrid1d *grid) {

  cgrid1d_fft(grid);
  cgrid1d_multiply(grid, grid->step);
}

/*
 * Perform inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       No normalization.
 *
 */

EXPORT void cgrid1d_inverse_fft(cgrid1d *grid) {

  if (!grid->iplan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in cgrid1d_inverse_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    cgrid1d_fftw_alloc(grid);
  }
  cgrid1d_fftw_inv(grid);
}

/*
 * Perform scaled inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid1d *).
 * c    = scaling factor (i.e. the output is multiplied by this constant) (double complex).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void cgrid1d_scaled_inverse_fft(cgrid1d *grid, double complex c) {

  cgrid1d_inverse_fft(grid);
  cgrid1d_multiply(grid, c);  
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by 1 / N.
 * (N = number of grid points)
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void cgrid1d_inverse_fft_norm(cgrid1d *grid) {

  cgrid1d_scaled_inverse_fft(grid, 1.0 / grid->nx);
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by 1 / (N step).
 * (N = number of grid points, step = grid step length)
 *
 * grid = grid to be inverse Fourier transformed (input/output) (cgrid1d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       This is the inverse routine to be used with cgrid1d_fourier_transform().
 *
 */

EXPORT void cgrid1d_inverse_fourier_transform(cgrid1d *grid) {

  cgrid1d_scaled_inverse_fft(grid, 1.0 / (grid->nx * grid->step));
}

/*
 * Convolue FFT transformed grids. To apply this on two grids (grida and gridb)
 * and place the result in gridc:
 * cgrid1d_fft(grida);
 * cgrid1d_fft(gridb);
 * cgrid1d_convolue(gridc, grida, gridb);
 * cgrid1d_inverse_fft(gridc);
 * gridc now contains the convolution of grida and gridb.
 *
 * grida = 1st grid to be convoluted (cgrid1d *).
 * gridb = 2nd grid to be convoluted (cgrid1d *).
 * gridc = output (cgrid1d *).
 *
 * No return value.
 *
 * Note: the input/output grids may be the same.
 *
 */

EXPORT void cgrid1d_fft_convolute(cgrid1d *gridc, const cgrid1d *grida, const cgrid1d *gridb) {

  long i, nx;
  double step, norm;
  double complex *cvalue, *bvalue;
  
  /* int f(r) g(r-r') dr' = iF[ F[f] F[g] ] = (step / N) iFFT[ FFT[f] FFT[g] ] */
  nx = gridc->nx;
  step = gridc->step;
  
  if (gridc != grida)
    cgrid1d_copy(gridc, grida);
  
  cvalue = gridc->value;
  bvalue = gridb->value;
  
  norm = step / (double) nx;
  
#pragma omp parallel for firstprivate(nx,cvalue,bvalue,norm) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    if (i & 1)
      cvalue[i] *= -norm * bvalue[i];
    else
      cvalue[i] *= norm * bvalue[i];
}

/*
 * Differentiate grid in the Fourier space.
 *
 * grid       = grid to be differentiated (in Fourier space) (cgrid1d *).
 * gradient_x = output grid (cgrid1d *).
 *
 * No return value.
 *
 * Note: input and output grids may be the same.
 *
 */

EXPORT void cgrid1d_fft_gradient(const cgrid1d *grid, cgrid1d *gradient_x) {

  cgrid1d_fft_gradient_x(grid, gradient_x);
}

EXPORT void cgrid1d_fft_gradient_x(const cgrid1d *grid, cgrid1d *gradient_x) {

  long i, nx;
  double kx, step, norm;
  double complex *gxvalue = gradient_x->value;
  
  /* f'(x) = iF[ i kx F[f(x)] ] */
  nx = grid->nx;
  step = grid->step;
  
  norm = 1.0 / (double) nx;
  
  if (gradient_x != grid)
    cgrid1d_copy(gradient_x, grid);
  
#pragma omp parallel for firstprivate(norm,nx,step,gxvalue) private(i,kx) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {
    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (i < nx / 2)
      kx = 2.0 * M_PI * i / (nx * step);
    else 
      kx = 2.0 * M_PI * (i - nx) / (nx * step);
    
    gxvalue[i] *= kx * norm * I;
  }
}

/* 
 * Calculate second derivative of a grid (in Fourier space).
 *
 * grid    = grid to be differentiated (cgrid1d *).
 * laplace = output grid (cgrid1d *).
 *
 * No return value.
 *
 * Note: input/output grids may be the same.
 *
 */

EXPORT void cgrid1d_fft_laplace(const cgrid1d *grid, cgrid1d *laplace) {

  long i, nx;
  double kx, step, norm;
  double complex *lvalue = laplace->value;
  
  nx = grid->nx;
  step = grid->step;
  
  norm = 1.0 / nx;
  
  if (grid != laplace)
    cgrid1d_copy(laplace, grid);
  
#pragma omp parallel for firstprivate(norm,nx,step,lvalue) private(i,kx) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {
    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (i < nx / 2)
      kx = 2.0 * M_PI * i / (nx * step);
    else 
      kx = 2.0 * M_PI * (i - nx) / (nx * step);
    lvalue[i] *= -kx * kx * norm;
  }
}

/*
 * Calculate expectation value of laplace operator in the Fourier space (int grid^* grid'').
 *
 * grid    = source grid for the operation (in Fourier space) (cgrid1d *).
 * laplace = laplacian of the grid (input) (cgrid1d *).
 *
 * Returns the expectation value (double).
 *
 */

EXPORT double cgrid1d_fft_laplace_expectation_value(const cgrid1d *grid, cgrid1d *laplace) {

  long i, nx;
  double kx, step, norm, sum = 0;
  double complex *lvalue = laplace->value;
  
  nx = grid->nx;
  step = grid->step;
  
  /* int (delta FFT[f(x)])^2 dk => delta^2 / N delta */
  norm = step / nx;
  
  if (grid != laplace)
    cgrid1d_copy(laplace, grid);
  
#pragma omp parallel for firstprivate(norm,nx,step,lvalue) private(i,kx) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {
    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (i < nx / 2)
      kx = 2.0 * M_PI * i / (nx * step);
    else 
      kx = 2.0 * M_PI * (i - nx) / (nx * step);
    sum += -kx * kx * sqnorm(lvalue[i]);
  }

  return sum * norm;
}

/* Boundary condition routines */

EXPORT double complex cgrid1d_value_outside_constantdirichlet(const cgrid1d *grid, long i) {

#if __DEBUG__
  if (i >= 0 && i < grid->nx) {
    fprintf(stderr, "libgrid: Error in cgrid1d_value_outside_constantdirichlet(). "
	    "Index (%ld) is not outside the grid.", i);
    abort();
  }
#endif
  
  return  *((double complex *) grid->outside_params_ptr);
}

EXPORT double complex cgrid1d_value_outside_neumann(const cgrid1d *grid, long i) {

  long nx = grid->nx;

#if __DEBUG__
  if (i >= 0 && i < grid->nx) {
    fprintf(stderr, "libgrid: Error in cgrid1d_value_outside_neumann(). "
	    "Index (%ld) is not outside the grid.", i);
    abort();
  }
#endif

  if(i < 0) i = 1;
  if(i > nx-1) i = nx - 2;

  return grid->value[i];
  
  /* This was the old mirror Neumann code */
#if 0
  nd = nx * 2;
  if (i < 0)   i  = -i;
  if (i >= nd) i %= nd;
  if (i >= nx) i  = nd - i;
  
  return grid->value[i];
#endif
}

EXPORT double complex cgrid1d_value_outside_periodic(const cgrid1d *grid, long i) {

  long nx = grid->nx;
  
#if __DEBUG__
  if (i >= 0 && i < grid->nx) {
    fprintf(stderr, "libgrid: Error in cgrid1d_value_outside_periodic(). "
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
 * grid = grid to be accessed (cgrid1d *).
 * i    = index (long).
 *
 * Returns grid value at index i.
 *
 */

EXPORT inline double complex cgrid1d_value_at_index(const cgrid1d *grid, long i) {

  if (i < 0 || i >= grid->nx)
    return grid->value_outside(grid, i);
  return grid->value[i];
}

/*
 * Access grid point at given x value (with linear interpolation).
 *
 * grid = grid to be accessed (cgrid1d *).
 * x    = x value (double).
 *
 * Returns grid value at x.
 *
 */

EXPORT inline double complex cgrid1d_value(const cgrid1d *grid, double x) {

  double complex f0, f1;
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
  f0 = cgrid1d_value_at_index(grid, i);
  f1 = cgrid1d_value_at_index(grid, i+1);
  
  omx = 1 - x;

  return omx * f0 + x * f1;
}

/*
 * Extrapolate between two different grid sizes.
 *
 * dest = Destination grid (cgrid1d *; output).
 * src  = Source grid (cgrid1d *; input).
 *
 */

EXPORT void cgrid1d_extrapolate(cgrid1d *dest, cgrid1d *src) {

  long i, nx = dest->nx;
  double step = dest->step, x;

  for (i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    dest->value[i] = cgrid1d_value(src, x);
  }
}

/*
 * Clear real part of complex grid.
 *
 */

EXPORT void cgrid1d_zero_re(cgrid1d *grid) {

  long i;

  for(i = 0; i < grid->nx; i++)
    grid->value[i] = I * cimag(grid->value[i]);
}

/*
 * Clear imaginary part of complex grid.
 *
 */

EXPORT void cgrid1d_zero_im(cgrid1d *grid) {

  long i;

  for(i = 0; i < grid->nx; i++)
    grid->value[i] = creal(grid->value[i]);
}
