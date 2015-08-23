/*
 * Routines for 2D real grids.
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
 *                      RGRID2D_DIRICHLET_BOUNDARY: Dirichlet boundary.
 *                      or RGRID2D_NEUMANN_BOUNDARY: Neumann boundary.
 *                      or RGRID2D_PERIODIC_BOUNDARY: Periodic boundary.
 *                      or user supplied function with pointer to grid and
 *                         grid index as parameters to provide boundary access.
 * outside_params_ptr = pointer for passing parameters for the given boundary
 *                      access function. Use 0 to with the predefined boundary
 *                      functions (RGRID2D_* above).
 *
 * Return value: pointer to the allocated grid (rgrid2d *). Returns NULL on
 * error.
 *
 */

EXPORT rgrid2d *rgrid2d_alloc(long nx, long ny, double step, double (*value_outside)(const rgrid2d *grid, long i, long j), void *outside_params_ptr) {

  rgrid2d *grid;
  
  grid = (rgrid2d *) malloc(sizeof(rgrid2d));
  if (!grid) {
    fprintf( stderr, "libgrid: Error in rgrid2d_alloc(). Could not allocate memory for 2d grid.\n");
    return 0;
  }
  
  if (!(grid->value = (double *) fftw_malloc(2 * nx * (ny/2 + 1) * sizeof(double)))) {
    fprintf(stderr, "libgrid: Error in rgrid2d_alloc(). Could not allocate memory for rgrid2d->value.\n");
    free(grid);
    return 0;
  }
  if(!(grid->cint = (cgrid2d *) malloc(sizeof(cgrid2d)))) {
    fprintf(stderr, "libtrid: Error in rgrid2d_alloc(). Could not allocate memory for rgrid2d->cint->value.\n");
    free(grid);
  }

  grid->nx = nx;
  grid->ny = ny;
  grid->step = step;
  /* complex structure interface (after FFT) */
  grid->cint->value = (double complex *) grid->value;
  grid->cint->step = grid->step;
  grid->cint->nx = grid->nx;
  grid->cint->ny = grid->ny/2 + 1;
  /* value outside not set */
  
  if (value_outside)
    grid->value_outside = value_outside;
  else
    grid->value_outside = rgrid2d_value_outside_constantdirichlet;
  
  if (outside_params_ptr)
    grid->outside_params_ptr = outside_params_ptr;
  else {
    grid->default_outside_params = 0.0;
    grid->outside_params_ptr = &grid->default_outside_params;
  }
  
  grid->plan = grid->iplan = NULL;
  grid->plan_cyl = grid->iplan_cyl = NULL;
  grid->gh = NULL;
  
#if __DEBUG__
  allocated_grids++;
  fprintf(stderr, "libmeas(debug): %3d 2d grids allocated.\n", allocated_grids);
#endif

  rgrid2d_constant(grid, NAN);
  
  return grid;
}

/*
 * Free 2D grid.
 *
 * grid = pointer to 2D grid to be freed (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_free(rgrid2d *grid) {

  if (grid) {
    if (grid->value) fftw_free(grid->value);
    if (grid->cint) free(grid->cint);
    rgrid2d_fftw_free(grid);
    rgrid2d_fft_cylindrical_free(grid);
    free(grid);
  }
}

/* 
 * Write 2D grid on disk in binary format.
 *
 * grid = 2D grid to be written (rgrid2d *).
 * out  = file handle for the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_write(rgrid2d *grid, FILE *out) {

  fwrite(&grid->nx, sizeof(long), 1, out);
  fwrite(&grid->ny, sizeof(long), 1, out);
  fwrite(&grid->step, sizeof(double), 1, out);
  fwrite(grid->value, sizeof(double), grid->nx * grid->ny, out);
}

/* 
 * Read 2D grid from disk in binary format.
 *
 * grid = 2D grid to be read (rgrid2d *).
 * in   = file handle for reading the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 * Note: works for both Cartesian and Cylindrical 2D grids.
 *
 */

EXPORT void rgrid2d_read(rgrid2d *grid, FILE *in) {

  long nx, ny;
  double step;
  
  fread(&nx, sizeof(long), 1, in);
  fread(&ny, sizeof(long), 1, in);
  fread(&step, sizeof(double), 1, in);
  
  if (nx != grid->nx || ny != grid->ny || step != grid->step) {
    rgrid2d *tmp;

    fprintf(stderr, "libgrid: Grid in file has different size than grid in memory.\n");
    fprintf(stderr, "libgrid: Interpolating between grids.\n");
    if(!(tmp = rgrid2d_alloc(nx, ny, step, grid->value_outside, NULL))) {
      fprintf(stderr, "libgrid: Error allocating grid in rgrid2d_read().\n");
      abort();
    }
    fread(tmp->value, sizeof(double), nx * ny, in);
    rgrid2d_extrapolate(grid, tmp);
    rgrid2d_free(tmp);
    return;
  }
  
  fread(grid->value, sizeof(double), grid->nx * grid->ny, in);
}

/*
 * Copy 2D grid from one grid to another.
 *
 * copy = destination grid (rgrid2d *).
 * grid = source grid (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_copy(rgrid2d *copy, const rgrid2d *grid) {

  long i, nx = grid->nx, ny = grid->ny, bytes = grid->ny * sizeof(double);
  double *gvalue = grid->value;
  double *cvalue = copy->value;
  
  copy->nx = grid->nx;
  copy->ny = grid->ny;
  copy->step = grid->step;
  
#pragma omp parallel for firstprivate(nx,ny,bytes,gvalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    memmove(&cvalue[i * ny], &gvalue[i * ny], bytes);
}

/*
 * Shift 2D grid by given amount spatially.
 *
 * shifted = destination grid for the operation (rgrid2d *).
 * grid    = source grid for the operation (rgrid2d *).
 * x       = shift spatially by this amount in x (double).
 * y       = shift spatially by this amount in y (double).
 *
 * No return value.
 *
 * NOTE: Source and destination may be the same grid.
 *
 */

EXPORT void rgrid2d_shift(rgrid2d *shifted, const rgrid2d *grid, double x, double y) {

  sShiftParametersr2d params;

  /* shift by (x,y) i.e. current grid center to (x,y) */
  params.x = x;  params.y = y; params.grid = grid;
  rgrid2d_map(shifted, shift_rgrid2d, &params);
}

/* 
 * Zero 2D grid.
 *
 * grid = grid to be zeroed (rgrid2d *).
 *
 * No return value.
 * 
 */

EXPORT void rgrid2d_zero(rgrid2d *grid) { 

  rgrid2d_constant(grid, 0.0); 
}

/* 
 * Set 2D grid to a constant value.
 *
 * grid = grid to be set (rgrid2d *).
 * c    = value (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_constant(rgrid2d *grid, double c) {

  long ij, nxy = grid->nx * grid->ny;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,c) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] = c;
}

/*
 * Multiply a given grid by a function.
 *
 * grid = destination grid for the operation (rgrid2d *).
 * func = function providing the mapping (double (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and x,y are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_product_func(rgrid2d *grid, double (*func)(void *arg, double x, double y), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double x, y, step = grid->step;
  double *value = grid->value;
  
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
 * grid = destination grid for the operation (rgrid2d *).
 * func = function providing the mapping (double (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and x,y are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_map(rgrid2d *grid, double (*func)(void *arg, double x, double y), void *farg) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double x, y, step = grid->step;
  double *value = grid->value;
  
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
 * grid = destination grid for the operation (rgrid2d *).
 * func = function providing the mapping (double (*)(void *, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and (x,y) is the point (doubles) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 * ns   = number of intermediate points to be used in smoothing (int).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_smooth_map(rgrid2d *grid, double (*func)(void *arg, double x, double y), void *farg, int ns) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny;
  double xc, yc, step = grid->step;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nxy,ns,step,func,value) private(i,j,xc,yc) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;    
    xc = (i - nx/2.0) * step;
    yc = (j - ny/2.0) * step;
    value[ij] = linearly_weighted_integralr2d(func, farg, xc, yc, step, ns);
  }
}

/*
 * Map a given function onto 2D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid. Limits for intermediate steps and
 * tolerance can be given.
 *
 * grid   = destination grid for the operation (rgrid2d *).
 * func   = function providing the mapping (double (*)(void *, double, double, double)).
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

EXPORT void rgrid2d_adaptive_map(rgrid2d *grid, double (*func)(void *arg, double x, double y), void *farg, int min_ns, int max_ns, double tol) {

  long i, j, ij, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny, ns;
  double xc, yc, step = grid->step;
  double tol2 = tol * tol;
  double sum, sump;
  double *value = grid->value;
  
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
      sum  = linearly_weighted_integralr2d(func, farg, xc, yc, step, ns);
      sump = linearly_weighted_integralr2d(func, farg, xc, yc, step, ns+1);
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
 * gridc = destination grid (rgrid2d *).
 * grida = 1st of the grids to be added (rgrid2d *).
 * gridb = 2nd of the grids to be added (rgrid2d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid2d_sum(rgrid2d *gridc, const rgrid2d *grida, const rgrid2d *gridb) {

  long ij, nxy = gridc->nx * gridc->ny;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = avalue[ij] + bvalue[ij];
}

/* 
 * Subtract two grids ("gridc = grida - gridb").
 *
 * gridc = destination grid (rgrid2d *).
 * grida = 1st source grid (rgrid2d *).
 * gridb = 2nd source grid (rgrid2d *).
 *
 * No return value.
 *
 * Note: both source and destination may be the same.
 *
 */

EXPORT void rgrid2d_difference(rgrid2d *gridc, const rgrid2d *grida, const rgrid2d *gridb) {

  long ij, nxy = gridc->nx * gridc->ny;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = avalue[ij] - bvalue[ij];
}

/* 
 * Calculate product of two grids ("gridc = grida * gridb").
 *
 * gridc = destination grid (rgrid2d *).
 * grida = 1st source grid (rgrid2d *).
 * gridb = 2nd source grid (rgrid2d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid2d_product(rgrid2d *gridc, const rgrid2d *grida, const rgrid2d *gridb) {

  long ij, nxy = gridc->nx * gridc->ny;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = avalue[ij] * bvalue[ij];
}

/* 
 * Rise a grid to given power.
 *
 * gridb    = destination grid (rgrid2d *).
 * grida    = 1st source grid (rgrid2d *).
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

EXPORT void rgrid2d_power(rgrid2d *gridb, const rgrid2d *grida, double exponent) {

  long ij, nxy = gridb->nx * gridb->ny;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,exponent) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    bvalue[ij] = pow(avalue[ij], exponent);
}

/* 
 * Rise absolute value of a grid to given power.
 *
 * gridb    = destination grid (rgrid2d *).
 * grida    = 1st source grid (rgrid2d *).
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

EXPORT void rgrid2d_abs_power(rgrid2d *gridb, const rgrid2d *grida, double exponent) {

  long ij, nxy = gridb->nx * gridb->ny;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,exponent) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    bvalue[ij] = pow(fabs(avalue[ij]), exponent);
}

/*
 * Divide two grids ("gridc = grida / gridb").
 *
 * gridc = destination grid (rgrid2d *).
 * grida = 1st source grid (rgrid2d *).
 * gridb = 2nd source grid (rgrid2d *).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same.
 *
 */

EXPORT void rgrid2d_division(rgrid2d *gridc, const rgrid2d *grida, const rgrid2d *gridb) {
  
  long ij, nxy = gridc->nx * gridc->ny;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = avalue[ij] / bvalue[ij];
}

/*
 * Add a constant to a 2D grid.
 *
 * grid = grid where the constant is added (rgrid2d *).
 * c    = constant to be added (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_add(rgrid2d *grid, double c) {

  long ij, nxy = grid->nx * grid->ny;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,c) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] += c;
}

/*
 * Multiply grid by a constant.
 *
 * grid = grid to be multiplied (rgrid2d *).
 * c    = multiplier (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_multiply(rgrid2d *grid, double c) {

  long ij, nxy = grid->nx * grid->ny;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,c) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] *= c;
}

/* 
 * Add and multiply: grid = (grid + ca) * cm.
 *
 * grid = grid to be operated (rgrid2d *).
 * ca   = constant to be added (double).
 * cm   = multiplier (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_add_and_multiply(rgrid2d *grid, double ca, double cm) {

  long ij, nxy = grid->nx * grid->ny;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,ca,cm) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] = (value[ij] + ca) * cm;
}

/*
 * Multiply and add: grid = cm * grid + ca.
 *
 * grid = grid to be operated (rgrid2d *).
 * ca   = constant to be added (double).
 * cm   = multiplier (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_multiply_and_add(rgrid2d *grid, double cm, double ca) {

  long ij, nxy = grid->nx * grid->ny;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,ca,cm) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    value[ij] = value[ij] * cm + ca;
}

/* 
 * Add scaled grid: gridc = gridc + d * grida
 *
 * gridc = destination grid for the operation (rgrid2d *).
 * d     = multiplier for the operation (double).
 * grida = source grid for the operation (rgrid2d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid2d_add_scaled(rgrid2d *gridc, double d, const rgrid2d *grida) {

  long ij, nxy = gridc->nx * gridc->ny;
  double *avalue = grida->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(d,nxy,avalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] += d * avalue[ij];
}

/*
 * Add multiply two grids and a constant: gridc = gridc + c * grida * gridb.
 *
 * gridc = destination grid (rgrid2d *).
 * d     = constant multiplier (double).
 * grida = 1st source grid (rgrid2d *).
 * gridb = 2nd source grid (rgrid2d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid2d_add_scaled_product(rgrid2d *gridc, double d, const rgrid2d *grida, const rgrid2d *gridb) {
  
  long ij, nxy = gridc->nx * gridc->ny;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue,d) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] += d * avalue[ij] * bvalue[ij];
}

/*
 * Operate on a grid by a given operator: gridc = O(grida).
 *
 * gridc    = destination grid (rgrid2d *).
 * grida    = source grid (rgrid2d *).
 * operator = operator (double (*)(double)).
 *            (i.e. a function mapping a given R-number to another)
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid2d_operate_one(rgrid2d *gridc, const rgrid2d *grida, double (*operator)(double a)) {

  long ij, nxy = gridc->nx * gridc->ny;
  double *avalue = grida->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,cvalue,operator) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = operator(avalue[ij]);
}

/* 
 * Operate on two grids and place the result in third: gridc = O(grida, gridb).
 * where O is the operator.
 *
 * gridc    = destination grid (rgrid2d *).
 * grida    = 1s source grid (rgrid2d *).
 * gridb    = 2nd source grid (rgrid2d *).
 * operator = operator mapping grida and gridb (double (*)(double, double)).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid2d_operate_two(rgrid2d *gridc, const rgrid2d *grida, const rgrid2d *gridb, double (*operator)(double a, double b)) {

  long ij, nxy = gridc->nx * gridc->ny;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,cvalue,operator) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = operator(avalue[ij], bvalue[ij]);
}

/*
 * Operate on a grid by a given operator.
 *
 * grid     = grid to be operated (rgrid2d *).
 * operator = operator (void (*)(double *)).
 * 
 * No return value.
 *
 */

EXPORT void rgrid2d_transform_one(rgrid2d *grid, void (*operator)(double *a)) {

  long ij, nxy = grid->nx * grid->ny;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,operator) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    operator(&value[ij]);
}

/*
 * Operate on two separate grids by a given operator.
 *
 * grida    = grid to be operated (rgrid2d *).
 * gridb    = grid to be operated (rgrid2d *).
 * operator = operator (void (*)(double *)).
 * 
 * No return value.
 *
 */

EXPORT void rgrid2d_transform_two(rgrid2d *grida, rgrid2d *gridb, void (*operator)(double *a, double *b)) {

  long ij, nxy = grida->nx * grida->ny;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,bvalue,operator) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    operator(&avalue[ij], &bvalue[ij]);
}

/*
 * Integrate over a grid.
 *
 * grid = grid to be integrated (rgrid2d *).
 *
 * Returns the integral value (double).
 *
 */

EXPORT double rgrid2d_integral(const rgrid2d *grid) {

  long i, j, nx = grid->nx, ny = grid->ny;
  double sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,grid) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
	sum += rgrid2d_value_at_index(grid, i, j);
  return sum * grid->step * grid->step;
}

/*
 * Integrate over a grid with limits.
 *
 * grid = grid to be integrated (rgrid2d *).
 * xl   = lower limit for x (double).
 * xu   = upper limit for x (double).
 * yl   = lower limit for y (double).
 * yu   = upper limit for y (double).
 *
 * Returns the integral value (double).
 *
 */

EXPORT double rgrid2d_integral_region(const rgrid2d *grid, double xl, double xu, double yl, double yu) {

  long iu, il, i, ju, jl, j, nx = grid->nx, ny = grid->ny;
  double sum;
  double step = grid->step;
   
  il = xl / step + nx/2;
  iu = xu / step + nx/2;
  jl = yl / step + ny/2;
  ju = yu / step + ny/2;
  
  sum = 0.0;
#pragma omp parallel for firstprivate(il,iu,jl,ju,nx,ny,grid) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for (i = il; i <= iu; i++)
    for (j = jl; j <= ju; j++)
      sum += rgrid2d_value_at_index(grid, i, j);
  return sum * step * step;
}

/* 
 * Integrate over the grid squared (int grid^2).
 *
 * grid = grid to be integrated (rgrid2d *).
 *
 * Returns the integral (double).
 *
 */

EXPORT double rgrid2d_integral_of_square(const rgrid2d *grid) {

  long i, j, nx = grid->nx, ny = grid->ny;
  double sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,grid) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      sum += sqnorm(rgrid2d_value_at_index(grid, i, j));
  return sum * grid->step * grid->step;
}

/*
 * Calculate overlap between two grids (int grida gridb).
 *
 * grida = 1st grid (rgrid2d *).
 * gridb = 2nd grid (rgrid2d *).
 *
 * Returns the value of the overlap integral (double).
 *
 */

EXPORT double rgrid2d_integral_of_product(const rgrid2d *grida, const rgrid2d *gridb) {

  long i, j, nx = grida->nx, ny = grida->ny;
  double sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,grida,gridb) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      sum += rgrid2d_value_at_index(grida, i, j) * rgrid2d_value_at_index(gridb, i, j);
  return sum * grida->step * grida->step;
}

/*
 * Calculate expectation value of a grid over the probability density given by the other.
 * (int gridb grida^2).
 *
 * grida = grid giving the probability (grida^2) (rgrid2d *).
 * gridb = grid to be averaged (rgrid2d *).
 *
 * Returns the average value (double).
 *
 */

EXPORT double rgrid2d_grid_expectation_value(const rgrid2d *grida, const rgrid2d *gridb) {

  long i, j, nx = grida->nx, ny = grida->ny;
  double sum = 0.0;
  
#pragma omp parallel for firstprivate(nx,ny,grida,gridb) private(i,j) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      sum += sqnorm(rgrid2d_value_at_index(grida, i, j)) * rgrid2d_value_at_index(gridb, i, j);
  return sum * grida->step * grida->step;
}

/*
 * Calculate the expectation value of a function over a grid.
 * (int grida func grida = int func |grida|^2).
 *
 * func  = function to be averaged (double (*)(void *, double, double, double)).
 *         The arguments are: optional arg, grida(x,y), x, y.
 * grida = grid giving the probability (grida^2) (rgrid2d *).
 *
 * Returns the average value (double).
 *
 */
 
EXPORT double rgrid2d_grid_expectation_value_func(void *arg, double (*func)(void *arg, double val, double x, double y), const rgrid2d *grida) {
   
  long i, j, nx = grida->nx, ny = grida->ny;
  double sum = 0.0, tmp, step = grida->step, x, y;
  
#pragma omp parallel for firstprivate(nx,ny,grida,step,func,arg) private(x,y,i,j,tmp) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    x = (i - nx/2) * step;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step;
      tmp = rgrid2d_value_at_index(grida, i, j);
      sum += sqnorm(tmp) * func(arg, tmp, x, y);
    }
  }
  return sum * step * step;
}

/* 
 * Integrate over the grid multiplied by weighting function (int grid w(x)).
 *
 * grid   = grid to be integrated over (rgrid2d *).
 * weight = function defining the weight (double (*)(double, double)). The arguments are (x,y) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double rgrid2d_weighted_integral(const rgrid2d *grid, double (*weight)(void *farg, double x, double y), void *farg) {

  long i, j, nx = grid->nx, ny = grid->ny;
  double sum = 0.0, step = grid->step, x, y;
  
#pragma omp parallel for firstprivate(nx,ny,grid,step,weight,farg) private(x,y,i,j) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    x = (i - nx/2) * step;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step;
      sum += weight(farg, x, y) * rgrid2d_value_at_index(grid, i, j);
    }
  }
  return sum * step * step;
}

/* 
 * Integrate over square of the grid multiplied by weighting function (int grid^2 w(x,y)).
 *
 * grid   = grid to be integrated over (rgrid2d *).
 * weight = function defining the weight (double (*)(double, double)).
 *          The arguments are (x,y) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double rgrid2d_weighted_integral_of_square(const rgrid2d *grid, double (*weight)(void *farg, double x, double y), void *farg) {

  long i, j, nx = grid->nx, ny = grid->ny;
  double sum = 0.0, step = grid->step, x, y;
  
#pragma omp parallel for firstprivate(nx,ny,grid,step,weight,farg) private(x,y,i,j) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    x = (i - nx/2) * step;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step;
      sum += weight(farg, x, y) * sqnorm(rgrid2d_value_at_index(grid, i, j));
    }
  }
  return sum * step * step;
}

/* 
 * Differentiate a grid with respect to x.
 *
 * grid     = grid to be differentiated (rgrid2d *).
 * gradient = differentiated grid output (rgrid2d *).
 * 
 * No return value.
 *
 */

EXPORT void rgrid2d_fd_gradient_x(const rgrid2d *grid, rgrid2d *gradient) {

  long i, j, ij;
  long ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in rgrid2d_fd_gradient_x().\n");
    return;
  }
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta,grid) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = inv_delta * (rgrid2d_value_at_index(grid, i+1, j) - rgrid2d_value_at_index(grid, i-1, j));
  }
}

/* 
 * Differentiate a grid with respect to y.
 *
 * grid     = grid to be differentiated (rgrid2d *).
 * gradient = differentiated grid output (rgrid2d *).
 * 
 * No return value.
 *
 */

EXPORT void rgrid2d_fd_gradient_y(const rgrid2d *grid, rgrid2d *gradient) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in rgrid2d_fd_gradient_y().\n");
    return;
  }
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta,grid) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = inv_delta * (rgrid2d_value_at_index(grid, i, j+1) - rgrid2d_value_at_index(grid, i, j-1));
  }
}

/*
 * Calculate gradient of a grid.
 *
 * grid       = grid to be differentiated twice (rgrid2d *).
 * gradient_x = x output grid for the operation (rgrid2d *).
 * gradient_y = y output grid for the operation (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_fd_gradient(const rgrid2d *grid, rgrid2d *gradient_x, rgrid2d *gradient_y) {

  rgrid2d_fd_gradient_x(grid, gradient_x);
  rgrid2d_fd_gradient_y(grid, gradient_y);
}

/*
 * Calculate laplacian of the grid.
 *
 * grid    = source grid (rgrid2d *).
 * laplace = output grid for the operation (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_fd_laplace(const rgrid2d *grid, rgrid2d *laplace) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_delta2 = 1.0 / (grid->step * grid->step);
  double *lvalue = laplace->value;
  
#pragma omp parallel for firstprivate(ny,nxy,lvalue,inv_delta2,grid) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = inv_delta2 * (-4.0 * rgrid2d_value_at_index(grid, i, j) + rgrid2d_value_at_index(grid, i, j+1) + rgrid2d_value_at_index(grid, i, j-1) + rgrid2d_value_at_index(grid,i+1,j) + rgrid2d_value_at_index(grid,i-1,j));
  }
}

/*
 * Calculate dot product of the gradient of the grid.
 *
 * grid          = source grid for the operation (rgrid2d *).
 * grad_dot_grad = destination grid (rgrid2d *).
 *
 * No return value.
 *
 * Note: grid and grad_dot_grad may not be the same grid.
 *
 */

EXPORT void rgrid2d_fd_gradient_dot_gradient(const rgrid2d *grid, rgrid2d *grad_dot_grad) {

  long i, j, ij, ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_2delta2 = 1.0 / (2.0 * grid->step * 2.0 * grid->step);
  double *gvalue = grad_dot_grad->value;
  
/*  grad f(x,y,z) dot grad f(x,y,z) = [ |f(+,0) - f(-,0)|^2 + |f(0,+) - f(0,-)|^2] / (2h)^2 */
#pragma omp parallel for firstprivate(ny,nxy,gvalue,inv_2delta2,grid) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    gvalue[ij] = inv_2delta2 * (sqnorm(rgrid2d_value_at_index(grid, i, j+1) - rgrid2d_value_at_index(grid, i, j-1)) + sqnorm(rgrid2d_value_at_index(grid, i+1, j) - rgrid2d_value_at_index(grid, i-1, j)));
  }
}

/*
 * Print the grid with both real and imaginary parts into file (ASCII format).
 *
 * grid = grid to be printed out (rgrid2d *).
 * out  = output file pointer (FILE *).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_print(const rgrid2d *grid, FILE *out) {

  long i, j;
  
  for( i = 0; i < grid->nx; i++ ) {
    for( j = 0; j < grid->ny; j++ ) {
      fprintf(out, "%16.8le  ", rgrid2d_value_at_index(grid, i, j));
    }
    fprintf(out, "\n");
  }
}

/*
 * Perform Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (rgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       Also no normalization is performed.
 *
 */

EXPORT void rgrid2d_fft(rgrid2d *grid) {

  if (!grid->plan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in rgrid2d_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    rgrid2d_fftw_alloc(grid);
  }  
  rgrid2d_fftw(grid);
}

/*
 * Perform normalized Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (rgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void rgrid2d_fourier_transform(rgrid2d *grid) {

  rgrid2d_fft(grid);
  cgrid2d_multiply(grid->cint, grid->step * grid->step);
}

/*
 * Perform inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output. No normalization.
 *
 */

EXPORT void rgrid2d_inverse_fft(rgrid2d *grid) {

  if (!grid->iplan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in rgrid2d_inverse_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    rgrid2d_fftw_alloc(grid);
  }
  rgrid2d_fftw_inv(grid);
}

/*
 * Perform scaled inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid2d *).
 * c    = scaling factor (i.e. the output is multiplied by this constant) (double).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void rgrid2d_scaled_inverse_fft(rgrid2d *grid, double c) {

  rgrid2d_inverse_fft(grid);
  rgrid2d_multiply(grid, c);  
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by 1 / (Nx * Ny). Note that Nx, Ny are the logical sizes of the arrays.
 * (Ni = number of grid points along axis i)
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void rgrid2d_inverse_fft_norm(rgrid2d *grid) {

  rgrid2d_scaled_inverse_fft(grid, 1.0 / (grid->nx * grid->ny));
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by 1 / (Nx*step * Ny*step). Nx, NY are the logical sizes.
 * (Ni = number of grid points along axis i, step = grid step length)
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid2d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       This is the inverse routine to be used with rgrid2d_fourier_transform().
 *
 */

EXPORT void rgrid2d_inverse_fourier_transform(rgrid2d *grid) {

  rgrid2d_scaled_inverse_fft(grid, 1.0 / (grid->nx * grid->ny * grid->step * grid->step));
}

/*
 * Convolute FFT transformed grids. To apply this on two grids (grida and gridb)
 * and place the result in gridc:
 * rgrid2d_fft(grida);
 * rgrid2d_fft(gridb);
 * rgrid2d_convolue(gridc, grida, gridb);
 * rgrid2d_inverse_fft(gridc);
 * gridc now contains the convolution of grida and gridb.
 *
 * grida = 1st grid to be convoluted (rgrid2d *).
 * gridb = 2nd grid to be convoluted (rgrid2d *).
 * gridc = output (rgrid2d *).
 *
 * No return value.
 *
 * Note: the input/output grids may be the same.
 *
 */

EXPORT void rgrid2d_fft_convolute(rgrid2d *gridc, const rgrid2d *grida, const rgrid2d *gridb) {

  long i, j, ij, nx, ny, nxy;
  double step, norm;
  double complex *avalue, *cvalue, *bvalue;
  
  /* int f(r) g(r-r') d^2r' = iF[ F[f] F[g] ] = (step / N)^2 iFFT[ FFT[f] FFT[g] ] */
  nx = gridc->cint->nx;
  ny = gridc->cint->ny;
  nxy = nx * ny;
  step = gridc->cint->step;
  
  avalue = grida->cint->value;
  bvalue = gridb->cint->value;
  cvalue = gridc->cint->value;
  
  norm = step * step / (double) (gridc->nx * gridc->ny);  
#pragma omp parallel for firstprivate(nx,ny,nxy,avalue,bvalue,cvalue,norm) private(i,j,ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    /* if odd */
    if ((i + j) & 1)
      cvalue[ij] = -norm * avalue[ij] * bvalue[ij];
    else
      cvalue[ij] = norm * avalue[ij] * bvalue[ij];
  }
}

/* Boundary condition routines */

EXPORT double rgrid2d_value_outside_constantdirichlet(const rgrid2d *grid, long i, long j) {

  return *((double *) grid->outside_params_ptr);
}

EXPORT double rgrid2d_value_outside_neumann(const rgrid2d *grid, long i, long j) {

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

EXPORT double rgrid2d_value_outside_periodic(const rgrid2d *grid, long i, long j) {

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
 * grid = grid to be accessed (rgrid2d *).
 * i    = index along x (long).
 * j    = index along y (long).
 *
 * Returns grid value at index (i, j).
 *
 */

EXPORT inline double rgrid2d_value_at_index(const rgrid2d *grid, long i, long j) {

  if (i < 0 || j < 0 || i >= grid->nx || j >= grid->ny)
    return grid->value_outside(grid, i, j);
  return grid->value[i * grid->ny + j];
}

/*
 * Access grid point at given (x,y) point (with linear interpolation).
 *
 * grid = grid to be accessed (rgrid2d *).
 * x    = x value (double).
 * y    = y value (double).
 *
 * Returns grid value at (x,y).
 *
 */

EXPORT inline double rgrid2d_value(const rgrid2d *grid, double x, double y) {

  double f00, f10, f01, f11;
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
  f00 = rgrid2d_value_at_index(grid, i, j);
  f10 = rgrid2d_value_at_index(grid, i+1, j);
  f01 = rgrid2d_value_at_index(grid, i, j+1);
  f11 = rgrid2d_value_at_index(grid, i+1, j+1);
  
  omx = 1 - x;
  omy = 1 - y;

  return omx * omy * f00 + x * omy * f10 + omx * y * f01 + x * y * f11;
}

/*
 * Copy a real grid to a complex grid (to real part).
 *
 * dest   = destination grid (cgrid2d *).
 * source = source grid (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_real_to_complex_re(cgrid2d *dest, rgrid2d *source) {
  
  long ij, nxy = source->nx * source->ny;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,dst,src) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    dst[ij] = (double complex) src[ij];
}

/*
 * Copy a real grid to a complex grid (to imaginary part).
 *
 * dest   = destination grid (cgrid2d *).
 * source = source grid (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_real_to_complex_im(cgrid2d *dest, rgrid2d *source) {
  
  long ij, nxy = source->nx * source->ny;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,dst,src) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    dst[ij] = I * (double complex) src[ij];
}

/*
 * Add a real grid to a complex grid (to real part).
 *
 * dest   = destination grid (cgrid2d *).
 * source = source grid (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_add_real_to_complex_re(cgrid2d *dest, rgrid2d *source) {
  
  long ij, nxy = source->nx * source->ny;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,dst,src) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    dst[ij] += (double complex) src[ij];
}

/*
 * Add a real grid to a complex grid (to imaginary part).
 *
 * dest   = destination grid (cgrid2d *).
 * source = source grid (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_add_real_to_complex_im(cgrid2d *dest, rgrid2d *source) {
  
  long ij, nxy = source->nx * source->ny;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,dst,src) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    dst[ij] += I * (double complex) src[ij];
}

/*
 * Product of a real grid with a complex grid.
 *
 * dest   = destination grid (cgrid2d *).
 * source = source grid (rgrid2d *).
 *
 * "dest(complex) = dest(complex) * source(real)"
 *
 * No return value.
 *
 */

EXPORT void grid2d_product_complex_with_real(cgrid2d *dest, rgrid2d *source) {
  
  long ij, nxy = source->nx * source->ny;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,dst,src) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    dst[ij] *= (double complex) src[ij];
}

/*
 * Copy imaginary part of a complex grid to a real grid.
 *
 * dest   = destination grid (rgrid2d *).
 * source = source grid (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_complex_im_to_real(rgrid2d *dest, cgrid2d *source) {
  
  long ij, nxy = source->nx * source->ny;
  double complex *src = source->value;
  double *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,dst,src) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    dst[ij] = cimag(src[ij]);
}

/*
 * Copy real part of a complex grid to a real grid.
 *
 * dest   = destination grid (rgrid2d *).
 * source = source grid (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_complex_re_to_real(rgrid2d *dest, cgrid2d *source) {
  
  long ij, nxy = source->nx * source->ny;
  double complex *src = source->value;
  double *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,dst,src) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    dst[ij] = creal(src[ij]);
}

/*
 * Extrapolate between two different grid sizes.
 *
 * dest = Destination grid (rgrid2d *; output).
 * src  = Source grid (rgrid2d *; input).
 *
 */

EXPORT void rgrid2d_extrapolate(rgrid2d *dest, rgrid2d *src) {

  long i, j, nx = dest->nx, ny = dest->ny;
  double step = dest->step, x, y;

  for (i = 0; i < nx; i++) {
    x = (i - nx/2.0) * step;
    for (j = 0; j < ny; j++) {
      y = (j - ny/2.0) * step;
      dest->value[i * ny + j] = rgrid2d_value(src, x, y);
    }
  }
}

/*
 * Get the largest value contained in a grid
 */
EXPORT double rgrid2d_max(rgrid2d *grid){
	long nxy = grid->nx * grid->ny ;
	double *val = grid->value , max_val = val[0] ;
	long i ;

#pragma omp parallel for firstprivate(nxy, val) private(i) reduction(max: max_val) default(none) schedule(runtime)
	for(i=0 ; i< nxy ; i++){
		if(val[i] > max_val)
			max_val = val[i] ;
	}
	return max_val ;
}

/*
 * Get the lowest value contained in a grid
 */
EXPORT double rgrid2d_min(rgrid2d *grid){
	long nxy = grid->nx * grid->ny ;
	double *val = grid->value , min_val = val[0] ;
	long i ;

#pragma omp parallel for firstprivate(nxy, val) private(i) reduction(min: min_val) default(none) schedule(runtime)
	for(i=0 ; i< nxy ; i++){
		if(val[i] < min_val)
			min_val = val[i] ;
	}
	return min_val ;
}
