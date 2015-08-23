  /*
 * Routines for 3D real grids.
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
 * Allocate 3D real grid.
 *
 * nx                 = number of points on the grid along x (long).
 * ny                 = number of points on the grid along y (long).
 * nz                 = number of points on the grid along z (long).
 * step               = spatial step length on the grid (double).
 * value_outside      = condition for accessing boundary points:
 *                      RGRID3D_DIRICHLET_BOUNDARY: Dirichlet boundary.
 *                      or RGRID3D_NEUMANN_BOUNDARY: Neumann boundary.
 *                      or RGRID3D_PERIODIC_BOUNDARY: Periodic boundary.
 *                      or user supplied function with pointer to grid and
 *                         grid index as parameters to provide boundary access.
 * outside_params_ptr = pointer for passing parameters for the given boundary
 *                      access function. Use 0 to with the predefined boundary
 *                      functions (void *).
 *
 * Return value: pointer to the allocated grid (rgrid3d *). Returns NULL on
 * error.
 *
 */

EXPORT rgrid3d *rgrid3d_alloc(long nx, long ny, long nz, double step, double (*value_outside)(const rgrid3d *grid, long i, long j, long k), void *outside_params_ptr) {

  rgrid3d *grid;
  
  if (!(grid = (rgrid3d *) malloc(sizeof(rgrid3d)))) {
    fprintf( stderr, "libgrid: Error in rgrid3d_alloc(). Could not allocate memory for 3d grid.\n");
   return 0;
  }
  
  if (!(grid->value = (double *) fftw_malloc(nx * ny * (nz/2 + 1) * sizeof(double complex)))) {  /* Extra space needed to hold the FFT data */
    fprintf(stderr, "libgrid: Error in rgrid3d_alloc(). Could not allocate memory for rgrid3d->value.\n");
    free(grid);
    return 0;
  }
  if(!(grid->cint = (cgrid3d *) malloc(sizeof(cgrid3d)))) {
    fprintf(stderr, "libtrid: Error in rgrid3d_alloc(). Could not allocate memory for rgrid3d->cint->value.\n");
    free(grid);
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
  /* complex structure interface (after FFT) */
  grid->cint->value = (double complex *) grid->value;
  grid->cint->step = grid->step;
  grid->cint->nx = grid->nx;
  grid->cint->ny = grid->ny;
  grid->cint->nz = grid->nz/2 + 1;
  /* value outside not set */
  
  if (value_outside)
    grid->value_outside = value_outside;
  else
    grid->value_outside = rgrid3d_value_outside_constantdirichlet;

  if (outside_params_ptr)
    grid->outside_params_ptr = outside_params_ptr;
  else {
    grid->default_outside_params = 0.0;
    grid->outside_params_ptr = &grid->default_outside_params;
  }
  
  grid->plan = grid->iplan = NULL;
  
#if __DEBUG__
  allocated_grids++;
  fprintf(stderr, "libgrid(debug): %3d 3d real grids allocated.\n", allocated_grids);
#endif
  
  rgrid3d_constant(grid, NAN);

  return grid;
}

/*
 * Set the origin of coordinates, meaning the coordinates of the grid will be:
 * 	x(i)  = (i - nx/2)* step - x0
 * 	y(j)  = (j - ny/2)* step - y0
 * 	z(k)  = (k - nz/2)* step - z0
 */
EXPORT void rgrid3d_set_origin( rgrid3d *grid , double x0, double y0, double z0){
	grid->x0 = x0 ;
	grid->y0 = y0 ;
	grid->z0 = z0 ;
}

/* Shift the origin */
EXPORT void rgrid3d_shift_origin( rgrid3d *grid , double x0, double y0, double z0){
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
EXPORT void rgrid3d_set_momentum( rgrid3d *grid , double kx0, double ky0, double kz0){
	grid->kx0 = kx0 ;
	grid->ky0 = ky0 ;
	grid->kz0 = kz0 ;
}

/*
 * Free 3D grid.
 *
 * grid = pointer to 3D grid to be freed (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_free(rgrid3d *grid) {

  if (grid) {
    if (grid->value) fftw_free(grid->value);
    if (grid->cint) free(grid->cint);
    rgrid3d_fftw_free(grid);
    free(grid);
  }
}

/* 
 * Write 3D grid on disk in binary format.
 *
 * grid = 3D grid to be written (rgrid3d *).
 * out  = file handle for the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_write(rgrid3d *grid, FILE *out) {

  fwrite(&grid->nx, sizeof(long), 1, out);
  fwrite(&grid->ny, sizeof(long), 1, out);
  fwrite(&grid->nz, sizeof(long), 1, out);
  fwrite(&grid->step, sizeof(double), 1, out);
  fwrite(grid->value, sizeof(double), grid->nx * grid->ny * grid->nz, out);
}

/* 
 * Read 3D grid from disk in binary format.
 *
 * grid = 3D grid to be read (rgrid3d *).
 * in   = file handle for reading the file (FILE * as defined in stdio.h).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_read(rgrid3d *grid, FILE *in) {

  long nx, ny, nz;
  double step;
  
  fread(&nx, sizeof(long), 1, in);
  fread(&ny, sizeof(long), 1, in);
  fread(&nz, sizeof(long), 1, in);
  fread(&step, sizeof(double), 1, in);
  
  if (nx != grid->nx || ny != grid->ny || nz != grid->nz || step != grid->step) {
    rgrid3d *tmp;

    fprintf(stderr, "libgrid: Grid in file has different size than grid in memory.\n");
    fprintf(stderr, "libgrid: Interpolating between grids.\n");
    if(!(tmp = rgrid3d_alloc(nx, ny, nz, step, grid->value_outside, NULL))) {
      fprintf(stderr, "libgrid: Error allocating grid in rgrid3d_read().\n");
      abort();
    }
    fread(tmp->value, sizeof(double), nx * ny * nz, in);
    rgrid3d_extrapolate(grid, tmp);
    rgrid3d_free(tmp);
    return;
  }
  
  fread(grid->value, sizeof(double), grid->nx * grid->ny * grid->nz, in);
}

/*
 * Copy 3D grid from one grid to another.
 *
 * copy = destination grid (rgrid3d *).
 * grid = source grid (rgrid3d *).
 *
 * No return value.
 *
 * NOTE: If FFT transformed grid needs to be copied,
 * one must use the grid's cint interface and cgrid3d_copy
 * to do so.
 *
 */

EXPORT void rgrid3d_copy(rgrid3d *copy, const rgrid3d *grid) {

  long i, nx = grid->nx, nyz = grid->ny * grid->nz, bytes = grid->ny * grid->nz * sizeof(double);
  double *gvalue = grid->value;
  double *cvalue = copy->value;
  
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
 * Shift 3D grid by given amount spatially.
 *
 * shifted = destination grid for the operation (rgrid3d *).
 * grid    = source grid for the operation (rgrid3d *).
 * x       = shift spatially by this amount in x (double).
 * y       = shift spatially by this amount in y (double).
 * z       = shift spatially by this amount in z (double).
 *
 * No return value.
 *
 * NOTE: Source and destination may be the same grid.
 *
 */

EXPORT void rgrid3d_shift(rgrid3d *shifted, const rgrid3d *grid, double x, double y, double z) {

  sShiftParametersr3d params;

  /* shift by (x,y,z) i.e. current grid center to (x,y,z) */
  params.x = x;  params.y = y;  params.z = z;  params.grid = grid;
  rgrid3d_map(shifted, shift_rgrid3d, &params);
}

/* 
 * Zero 3D grid.
 *
 * grid = grid to be zeroed (rgrid3d *).
 *
 * No return value.
 * 
 */

EXPORT void rgrid3d_zero(rgrid3d *grid) { 

  rgrid3d_constant(grid, 0.0); 
}

/* 
 * Set 3D grid to a constant value.
 *
 * grid = grid to be set (rgrid3d *).
 * c    = value (double).
 *
 * No return value.
 */

EXPORT void rgrid3d_constant(rgrid3d *grid, double c) {

   long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
   double *value = grid->value;
  
#pragma omp parallel for firstprivate(nxy,nz,value,c) private(ij,ijnz,k) default(none) schedule(runtime)
   for(ij = 0; ij < nxy; ij++) {
     ijnz = ij * nz;
     for(k = 0; k < nz; k++)
       value[ijnz + k] = c;
  }
 }

/*
 * Multiply a given grid by a function.
 *
 * grid = destination grid for the operation (rgrid3d *).
 * func = function providing the mapping (double (*)(void *, double, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and x,y,z are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_product_func(rgrid3d *grid, double (*func)(void *arg, double x, double y, double z), void *farg) {

  long i, j, k, ij, ijnz, nx = grid->nx , ny = grid->ny, nxy = nx * ny, nz = grid->nz;
  double x, y, z, step = grid->step;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nz,nxy,step,func,value,x0,y0,z0) private(i,j,ij,ijnz,k,x,y,z) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    x = (i - nx/2) * step - x0;
    y = (j - ny/2) * step - y0;
    
    for(k = 0; k < nz; k++) {
      z = (k - nz/2) * step - z0;
      value[ijnz + k] *= func(farg, x, y, z);
    }
  }
}

/*
 * Map a given function onto 3D grid.
 *
 * grid = destination grid for the operation (rgrid3d *).
 * func = function providing the mapping (double (*)(void *, double, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and x,y,z are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_map(rgrid3d *grid, double (*func)(void *arg, double x, double y, double z), void *farg) {

  long i, j, k, ij, ijnz, nx = grid->nx , ny = grid->ny, nxy = nx * ny , nz = grid->nz;
  double x,y,z, step = grid->step;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  double *value = grid->value;
  
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
 * Map a given function onto 3D grid in k-space.
 *
 * grid = destination grid for the operation (rgrid3d *).
 * func = function providing the mapping (double (*)(void *, double, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and x,y,z are the coordinates (double) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_mapk(rgrid3d *grid, double (*func)(void *arg, double kx, double ky, double kz), void *farg) {

  long i, j, k, ij, ijnz, nxy = grid->nx * grid->ny, nx = grid->nx, ny = grid->ny, nz = grid->nz;
  double kx,ky,kz, step = grid->step;
  double kx0 = grid->kx0 , ky0 = grid->ky0 , kz0 = grid->kz0 ;
  double *value = grid->value;

  if(grid->value_outside == RGRID3D_PERIODIC_BOUNDARY ){
    #pragma omp parallel for firstprivate(farg,nx,ny,nz,nxy,step,func,value,kx0,ky0,kz0) private(i,j,ij,ijnz,k,kx,ky,kz) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      ijnz = ij * nz;
      i = ij / ny;
      j = ij % ny;
      if (i < nx / 2) kx = 2.0 * M_PI * i / (nx * step) - kx0;
        else kx = 2.0 * M_PI * (i - nx) / (nx * step) - kx0;
      if (j < ny / 2) ky = 2.0 * M_PI * j / (ny * step) - ky0;
        else ky = 2.0 * M_PI * (j - ny) / (ny * step) - ky0;
      for(k = 0; k < nz; k++) {
        if (k < nz / 2) kz = 2.0 * M_PI * k / (nz * step) - kz0;
          else kz = 2.0 * M_PI * (k - nz) / (nx * step) - kz0;
        value[ijnz + k] = func(farg, kx, ky, kz);
      }
    }
  }else{ 
    #pragma omp parallel for firstprivate(farg,nx,ny,nz,nxy,step,func,value,kx0,ky0,kz0) private(i,j,ij,ijnz,k,kx,ky,kz) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      ijnz = ij * nz;
      i = ij / ny;
      j = ij % ny;
      kx = M_PI * i / (nx * step) - kx0;
      ky = M_PI * j / (ny * step) - ky0;
      for(k = 0; k < nz; k++) {
        kz = M_PI * k / (nz * step) - kz0;
        value[ijnz + k] = func(farg, kx, ky, kz);
      }
    }
  }
}


/*
 * Map a given function onto 3D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid.
 * *
 * grid = destination grid for the operation (rgrid3d *).
 * func = function providing the mapping (double (*)(void *, double, double, double)).
 *        The first argument (void *) is for external user specified data
 *        and (x, y, z) is the point (doubles) where the function is evaluated.
 * farg = pointer to user specified data (void *).
 * ns   = number of intermediate points to be used in smoothing (int).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_smooth_map(rgrid3d *grid, double (*func)(void *arg, double x, double y, double z), void *farg, int ns) {

  long i, j, k, ij, ijnz, nx = grid->nx , ny = grid->ny, nxy = nx * ny , nz = grid->nz;
  double xc, yc, zc, step = grid->step;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  double *value = grid->value;
  
#pragma omp parallel for firstprivate(farg,nx,ny,nz,nxy,ns,step,func,value,x0,y0,z0) private(i,j,k,ijnz,xc,yc,zc) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    xc = (i - nx/2) * step - x0;
    yc = (j - ny/2) * step - y0;
    for(k = 0; k < nz; k++) {
      zc = (k - nz/2) * step - z0;
      value[ijnz + k] = linearly_weighted_integralr3d(func, farg, xc, yc, zc, step, ns);
    }
  }
}


/*
 * Map a given function onto 3D grid with linear "smoothing".
 * This can be used to weight the values at grid points to produce more
 * accurate integration over the grid. Limits for intermediate steps and
 * tolerance can be given.
 *
 * grid   = destination grid for the operation (rgrid3d *).
 * func   = function providing the mapping (double (*)(void *, double, double, double)).
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

EXPORT void rgrid3d_adaptive_map(rgrid3d *grid, double (*func)(void *arg, double x, double y, double z), void *farg, int min_ns, int max_ns, double tol) {

  long i, j, k, ij, ijnz, nx = grid->nx , ny = grid->ny, nxy = nx * ny, nz = grid->nz, ns;
  double xc, yc, zc, step = grid->step;
  double tol2 = tol * tol;
  double  sum, sump;
  double x0 = grid->x0, y0 = grid->y0, z0 = grid->z0;
  double  *value = grid->value;
  
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
        sum  = linearly_weighted_integralr3d(func, farg, xc, yc, zc, step, ns);
        sump = linearly_weighted_integralr3d(func, farg, xc, yc, zc, step, ns+1);
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
 * NON-PERIODIC version of map
 * Map a given function onto 3D grid in range [0,L] instead of [-L/2,L/2].
 * Can be removed along with all the map_nonperiodic_X functions.
 * No return value.
 *
 */
EXPORT void rgrid3d_map_nonperiodic(rgrid3d *grid, double (*func)(void *arg, double x, double y, double z), void *farg) {

	rgrid3d_set_origin(grid,
			-(grid->nx/2)*grid->step ,
			-(grid->ny/2)*grid->step ,
			-(grid->nz/2)*grid->step ) ;

	rgrid3d_map(grid, func, farg) ;
}	

/*
 * NON-PERIODIC version of smooth_map
 * Map a given function onto 3D grid in range [0,L] instead of [-L/2,L/2].
 * Can be removed along with all the map_nonperiodic_X functions.
 * No return value.
 *
 */
EXPORT void rgrid3d_smooth_map_nonperiodic(rgrid3d *grid, double (*func)(void *arg, double x, double y, double z), void *farg, int ns) {

	rgrid3d_set_origin(grid,
			-(grid->nx/2)*grid->step ,
			-(grid->ny/2)*grid->step ,
			-(grid->nz/2)*grid->step ) ;

	rgrid3d_smooth_map(grid, func, farg, ns) ;
}

/*
 * NON-PERIODIC version of adaptive_map
 * Map a given function onto 3D grid in range [0,L] instead of [-L/2,L/2].
 * Can be removed along with all the map_nonperiodic_X functions.
 * No return value.
 *
 */
EXPORT void rgrid3d_adaptive_map_nonperiodic(rgrid3d *grid, double (*func)(void *arg, double x, double y, double z), void *farg, int min_ns, int max_ns, double tol) {

	rgrid3d_set_origin(grid,
			-(grid->nx/2)*grid->step ,
			-(grid->ny/2)*grid->step ,
			-(grid->nz/2)*grid->step ) ;

	rgrid3d_adaptive_map(grid, func, farg, min_ns, max_ns, tol) ;
}

	
/*
 * Add two 3D grids ("gridc = grida + gridb").
 *
 * gridc = destination grid (rgrid3d *).
 * grida = 1st of the grids to be added (rgrid3d *).
 * gridb = 2nd of the grids to be added (rgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid3d_sum(rgrid3d *gridc, const rgrid3d *grida, const rgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
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
 * gridc = destination grid (rgrid3d *).
 * grida = 1st source grid (rgrid3d *).
 * gridb = 2nd source grid (rgrid3d *).
 *
 * No return value.
 *
 * Note: both source and destination may be the same.
 *
 */

EXPORT void rgrid3d_difference(rgrid3d *gridc, const rgrid3d *grida, const rgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;

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
 * gridc = destination grid (rgrid3d *).
 * grida = 1st source grid (rgrid3d *).
 * gridb = 2nd source grid (rgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid3d_product(rgrid3d *gridc, const rgrid3d *grida, const rgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
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
 * gridb    = destination grid (rgrid3d *).
 * grida    = 1st source grid (rgrid3d *).
 * exponent = exponent to be used (double).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same.
 *       This routine uses pow() so that the exponent can be
 *       fractional but this is slow! Do not use this for integer
 *       exponents.
 *
 */

EXPORT void rgrid3d_power(rgrid3d *gridb, const rgrid3d *grida, double exponent) {

  long ij, k, ijnz, nxy = gridb->nx * gridb->ny, nz = gridb->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
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
 * gridb    = destination grid (rgrid3d *).
 * grida    = 1st source grid (rgrid3d *).
 * exponent = exponent to be used (double).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same.
 *       This routine uses pow() so that the exponent can be
 *       fractional but this is slow! Do not use this for integer
 *       exponents.
 *
 */

EXPORT void rgrid3d_abs_power(rgrid3d *gridb, const rgrid3d *grida, double exponent) {

  long ij, k, ijnz, nxy = gridb->nx * gridb->ny, nz = gridb->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,exponent) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      bvalue[ijnz + k] = pow(fabs(avalue[ijnz + k]), exponent);
  }
}


/*
 * Divide two grids ("gridc = grida / gridb").
 *
 * gridc = destination grid (rgrid3d *).
 * grida = 1st source grid (rgrid3d *).
 * gridb = 2nd source grid (rgrid3d *).
 *
 * No return value.
 *
 * Note: Source and destination grids may be the same. EPS added to avoid
 * possible NaNs.
 *
 */

EXPORT void rgrid3d_division(rgrid3d *gridc, const rgrid3d *grida, const rgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,cvalue) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = avalue[ijnz + k] / (bvalue[ijnz + k] + GRID_EPS);
  }
}

/*
 * Add a constant to a 3D grid.
 *
 * grid = grid where the constant is added (rgrid3d *).
 * c    = constant to be added (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_add(rgrid3d *grid, double c) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double *value = grid->value;
  
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
 * grid = grid to be multiplied (rgrid3d *).
 * c    = multiplier (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_multiply(rgrid3d *grid, double c) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double *value = grid->value;
  
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
 * grid = grid to be operated (rgrid3d *).
 * ca   = constant to be added (double).
 * cm   = multiplier (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_add_and_multiply(rgrid3d *grid, double ca, double cm) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double *value = grid->value;
  
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
 * grid = grid to be operated (rgrid3d *).
 * ca   = constant to be added (double).
 * cm   = multiplier (double).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_multiply_and_add(rgrid3d *grid, double cm, double ca) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double *value = grid->value;
  
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
 * gridc = destination grid for the operation (rgrid3d *).
 * d     = multiplier for the operation (double).
 * grida = source grid for the operation (rgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid3d_add_scaled(rgrid3d *gridc, double d, const rgrid3d *grida) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double *avalue = grida->value;
  double *cvalue = gridc->value;
  
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
 * gridc = destination grid (rgrid3d *).
 * d     = constant multiplier (double).
 * grida = 1st source grid (rgrid3d *).
 * gridb = 2nd source grid (rgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid3d_add_scaled_product(rgrid3d *gridc, double d, const rgrid3d *grida, const rgrid3d *gridb) {
  
  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
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
 * gridc    = destination grid (rgrid3d *).
 * grida    = source grid (rgrid3d *).
 * operator = operator (double (*)(double)).
 *            (i.e., a function mapping a given C-number to another)
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid3d_operate_one(rgrid3d *gridc, const rgrid3d *grida, double (*operator)(double a)) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double *avalue = grida->value;
  double *cvalue = gridc->value;
  
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
 * gridc    = destination grid (rgrid3d *).
 * grida    = 1s source grid (rgrid3d *).
 * gridb    = 2nd source grid (rgrid3d *).
 * operator = operator mapping grida and gridb (double (*)(double, double)).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 *
 */

EXPORT void rgrid3d_operate_two(rgrid3d *gridc, const rgrid3d *grida, const rgrid3d *gridb, double (*operator)(double, double)) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  
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
 * grid     = grid to be operated (rgrid3d *).
 * operator = operator (void (*)(double *)).
 * 
 * No return value.
 *
 */

EXPORT void rgrid3d_transform_one(rgrid3d *grid, void (*operator)(double *a)) {

  long ij, k, ijnz, nxy = grid->nx * grid->ny, nz = grid->nz;
  double *value = grid->value;
  
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
 * grida    = grid to be operated (rgrid3d *).
 * gridb    = grid to be operated (rgrid3d *).
 * operator = operator (void (*)(double *)).
 * 
 * No return value.
 *
 */

EXPORT void rgrid3d_transform_two(rgrid3d *grida, rgrid3d *gridb, void (*operator)(double *a, double *b)) {

  long ij, k, ijnz, nxy = grida->nx * grida->ny, nz = grida->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  
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
 * grid = grid to be integrated (rgrid3d *).
 *
 * Returns the integral value (double).
 *
 * NOTE: This will integrate also over the missing points due to BC
 *       such that the symmetry is preserved.
 *
 */

EXPORT double rgrid3d_integral(const rgrid3d *grid) {

  long i, j, k, nx = grid->nx, ny = grid->ny, nz = grid->nz;
  double sum = 0.0;

  // TODO: collapse(2) ?
#pragma omp parallel for firstprivate(nx,ny,nz,grid) private(i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      for (k = 0; k <= nz; k++)
	sum += rgrid3d_value_at_index(grid, i, j, k);
  return sum * grid->step * grid->step * grid->step;
}

/*
 * Integrate over a grid with limits.
 *
 * grid = grid to be integrated (rgrid3d *).
 * xl   = lower limit for x (double).
 * xu   = upper limit for x (double).
 * yl   = lower limit for y (double).
 * yu   = upper limit for y (double).
 * zl   = lower limit for z (double).
 * zu   = upper limit for z (double).
 *
 * Returns the integral value (double).
 *
 */

EXPORT double rgrid3d_integral_region(const rgrid3d *grid, double xl, double xu, double yl, double yu, double zl, double zu) {

  long iu, il, i, ju, jl, j, ku, kl, k;
  double sum;
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
	sum += rgrid3d_value_at_index(grid, i, j, k);
  return sum * step * step * step; 
}
 
/* 
 * Integrate over the grid squared (int grid^2).
 *
 * grid = grid to be integrated (rgrid3d *).
 *
 * Returns the integral (double).
 *
 */

EXPORT double rgrid3d_integral_of_square(const rgrid3d *grid) {

  long i, j, k, nx = grid->nx, ny = grid->ny, nz = grid->nz;
  double sum = 0.0;
  
  // TODO: collapse(2) ?
#pragma omp parallel for firstprivate(nx,ny,nz,grid) private(i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      for (k = 0; k <= nz; k++)
	sum += sqnorm(rgrid3d_value_at_index(grid, i, j, k));
  return sum * grid->step * grid->step * grid->step;
}

/*
 * Calculate overlap between two grids (int grida gridb).
 *
 * grida = 1st grid (rgrid3d *).
 * gridb = 2nd grid (rgrid3d *).
 *
 * Returns the value of the overlap integral (double).
 *
 */

EXPORT double rgrid3d_integral_of_product(const rgrid3d *grida, const rgrid3d *gridb) {

  long i, j, k, nx = grida->nx, ny = grida->ny, nz = grida->nz;
  double sum = 0.0;
  
  // TODO: collapse(2) ?
#pragma omp parallel for firstprivate(nx,ny,nz,grida,gridb) private(i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      for (k = 0; k <= nz; k++)
	sum += rgrid3d_value_at_index(grida, i, j, k) * rgrid3d_value_at_index(gridb, i, j, k);
  return sum * grida->step * grida->step * grida->step;
}

/*
 * Calculate the expectation value of a grid over a grid.
 * (int gridb grida gridb = int grida gridb^2).
 *
 * grida = grid giving the probability (gridb^2) (rgrid3d *).
 * gridb = grid to be averaged (rgrid3d *).
 *
 * Returns the average value (double *).
 *
 */

EXPORT double rgrid3d_grid_expectation_value(const rgrid3d *grida, const rgrid3d *gridb) {

  long i, j, k, nx = grida->nx, ny = grida->ny, nz = grida->nz;
  double sum = 0.0;
  
  // TODO: collapse(2) ?
#pragma omp parallel for firstprivate(nx,ny,nz,grida,gridb) private(i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++)
    for (j = 0; j <= ny; j++)
      for (k = 0; k <= nz; k++)
	sum += sqnorm(rgrid3d_value_at_index(grida, i, j, k)) * rgrid3d_value_at_index(gridb, i, j, k);
  return sum * grida->step * grida->step * grida->step;
}
 
/*
 * Calculate the expectation value of a function over a grid.
 * (int grida func grida = int func grida^2).
 *
 * func  = function to be averaged (double (*)(void *, double, double, double, double)).
 *         The arguments are: optional arg, grida(x,y,z), x, y, z.
 * grida = grid giving the probability (grida^2) (rgrid3d *).
 *
 * Returns the average value (double).
 *
 */
 
EXPORT double rgrid3d_grid_expectation_value_func(void *arg, double (*func)(void *arg, double val, double x, double y, double z), const rgrid3d *grida) {
   
  long i, j, k, nx = grida->nx, ny = grida->ny, nz = grida->nz;
  double sum = 0.0, tmp, step = grida->step, x0 = grida->x0, y0 = grida->y0, z0 = grida->z0, x, y, z;
  
  // TODO: collapse(2) ? move x = ... inside the 2nd loop?
#pragma omp parallel for firstprivate(nx,ny,nz,grida,x0,y0,z0,step,func,arg) private(x,y,z,i,j,k,tmp) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    x = (i - nx/2) * step - x0;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step - y0;
      for (k = 0; k <= nz; k++) {
	z = (k - nz/2) * step - z0;
	tmp = rgrid3d_value_at_index(grida, i, j, k);
	sum += sqnorm(tmp) * func(arg, tmp, x, y, z);
      }
    }
  }
  return sum * step * step * step;
}

/* 
 * Integrate over the grid multiplied by weighting function (int grid w(x)).
 *
 * grid   = grid to be integrated over (rgrid3d *).
 * weight = function defining the weight (double (*)(double, double, double)). The arguments are (x,y,z) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double rgrid3d_weighted_integral(const rgrid3d *grid, double (*weight)(void *farg, double x, double y, double z), void *farg) {

  long i, j, k, nx = grid->nx, ny = grid->ny, nz = grid->nz;
  double sum = 0.0, step = grid->step, x0 = grid->x0, y0 = grid->y0, z0 = grid->z0, x, y, z;
  
  // TODO: collapse(2) ? move x = ... inside the 2nd loop?
#pragma omp parallel for firstprivate(nx,ny,nz,grid,x0,y0,z0,step,weight,farg) private(x,y,z,i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    x = (i - nx/2) * step - x0;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step - y0;
      for (k = 0; k <= nz; k++) {
	z = (k - nz/2) * step - z0;
	sum += weight(farg, x, y, z) * rgrid3d_value_at_index(grid, i, j, k);
      }
    }
  }
  return sum * step * step * step;
}

/* 
 * Integrate over square of the grid multiplied by weighting function (int grid^2 w(x)).
 *
 * grid   = grid to be integrated over (rgrid3d *).
 * weight = function defining the weight (double (*)(double, double, double)).
 *          The arguments are (x,y,z) coordinates.
 * farg   = argument to the weight function (void *).
 *
 * Returns the value of the integral (double).
 *
 */

EXPORT double rgrid3d_weighted_integral_of_square(const rgrid3d *grid, double (*weight)(void *farg, double x, double y, double z), void *farg) {

  long i, j, k, nx = grid->nx, ny = grid->ny, nz = grid->nz;
  double sum = 0.0, step = grid->step, x0 = grid->x0, y0 = grid->y0, z0 = grid->z0, x, y, z;
  
  // TODO: collapse(2) ? move x = ... inside the 2nd loop?
#pragma omp parallel for firstprivate(nx,ny,nz,grid,x0,y0,z0,step,weight,farg) private(x,y,z,i,j,k) reduction(+:sum) default(none) schedule(runtime)
  for (i = 0; i <= nx; i++) {
    x = (i - nx/2) * step - x0;
    for (j = 0; j <= ny; j++) {
      y = (j - ny/2) * step - y0;
      for (k = 0; k <= nz; k++) {
	z = (k - nz/2) * step - z0;
	sum += weight(farg, x, y, z) * sqnorm(rgrid3d_value_at_index(grid, i, j, k));
      }
    }
  }
  return sum * step * step * step;
}

/* 
 * Differentiate a grid with respect to x (central difference).
 *
 * grid     = grid to be differentiated (rgrid3d *).
 * gradient = differentiated grid output (rgrid3d *).
 * 
 * No return value.
 *
 */

EXPORT void rgrid3d_fd_gradient_x(const rgrid3d *grid, rgrid3d *gradient) {

  long i, j, k, ij, ijnz, ny = grid->ny, nz = grid->nz, nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in rgrid3d_fd_gradient_x().\n");
    return;
  }

#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = inv_delta * (rgrid3d_value_at_index(grid, i+1, j, k) - rgrid3d_value_at_index(grid, i-1, j, k));
  }
}

/* 
 * Differentiate a grid with respect to y.
 *
 * grid     = grid to be differentiated (rgrid3d *).
 * gradient = differentiated grid output (rgrid3d *).
 * 
 * No return value.
 *
 */

EXPORT void rgrid3d_fd_gradient_y(const rgrid3d *grid, rgrid3d *gradient) {

  long i, j, k, ij, ijnz, ny = grid->ny, nz = grid->nz, nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in rgrid3d_fd_gradient_x().\n");
    return;
  }

#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = inv_delta * (rgrid3d_value_at_index(grid, i, j+1, k) - rgrid3d_value_at_index(grid, i, j-1, k));
  }
}

/* 
 * Differentiate a grid with respect to z.
 *
 * grid     = grid to be differentiated (rgrid3d *).
 * gradient = differentiated grid output (rgrid3d *).
 * 
 * No return value.
 *
 */

EXPORT void rgrid3d_fd_gradient_z(const rgrid3d *grid, rgrid3d *gradient) {

  long i, j, k, ij, ijnz, ny = grid->ny, nz = grid->nz, nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step);
  double *lvalue = gradient->value;
  
  if(grid == gradient) {
    fprintf(stderr, "libgrid: source and destination must be different in rgrid3d_fd_gradient_x().\n");
    return;
  }

#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = inv_delta * (rgrid3d_value_at_index(grid, i, j, k+1) - rgrid3d_value_at_index(grid, i, j, k-1));
  }
}
 
/*
 * Calculate gradient of a grid.
 *
 * grid       = grid to be differentiated twice (rgrid3d *).
 * gradient_x = x output grid for the operation (rgrid3d *).
 * gradient_y = y output grid for the operation (rgrid3d *).
 * gradient_z = z output grid for the operation (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_fd_gradient(const rgrid3d *grid, rgrid3d *gradient_x, rgrid3d *gradient_y, rgrid3d *gradient_z) {

  rgrid3d_fd_gradient_x(grid, gradient_x);
  rgrid3d_fd_gradient_y(grid, gradient_y);
  rgrid3d_fd_gradient_z(grid, gradient_z);
}

/*
 * Calculate laplacian of the grid.
 *
 * grid    = source grid (rgrid3d *).
 * laplace = output grid for the operation (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_fd_laplace(const rgrid3d *grid, rgrid3d *laplace) {

  long i, j, k, ij, ijnz;
  long ny = grid->ny, nz = grid->nz;
  long nxy = grid->nx * grid->ny;
  double inv_delta2 = 1.0 / (grid->step * grid->step);
  double *lvalue = laplace->value;
  
#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta2,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = inv_delta2 * (-6.0 * rgrid3d_value_at_index(grid, i, j, k) + rgrid3d_value_at_index(grid, i, j, k+1)
				       + rgrid3d_value_at_index(grid, i, j, k-1) + rgrid3d_value_at_index(grid, i, j+1, k) 
				       + rgrid3d_value_at_index(grid, i, j-1, k) + rgrid3d_value_at_index(grid,i+1,j,k) 
				       + rgrid3d_value_at_index(grid,i-1,j,k));
  }
}

/*
 * Calculate dot product of the gradient of the grid.
 *
 * grid          = source grid for the operation (rgrid3d *).
 * grad_dot_grad = destination grid (rgrid3d *).
 *
 * No return value.
 *
 * Note: grid and grad_dot_grad may not be the same grid.
 *
 */

EXPORT void rgrid3d_fd_gradient_dot_gradient(const rgrid3d *grid, rgrid3d *grad_dot_grad) {

  long i, j, k, ij, ijnz;
  long ny = grid->ny, nz = grid->nz;
  long nxy = grid->nx * grid->ny;
  double inv_2delta2 = 1.0 / (2.0*grid->step * 2.0*grid->step);
  double *gvalue = grad_dot_grad->value;
  
/*  grad f(x,y,z) dot grad f(x,y,z) = [ |f(+,0,0) - f(-,0,0)|^2 + |f(0,+,0) - f(0,-,0)|^2 + |f(0,0,+) - f(0,0,-)|^2 ] / (2h)^2 */
#pragma omp parallel for firstprivate(ny,nz,nxy,gvalue,inv_2delta2,grid) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      gvalue[ijnz + k] = inv_2delta2 * 
	(sqnorm(rgrid3d_value_at_index(grid, i, j, k+1) - rgrid3d_value_at_index(grid, i, j, k-1)) + sqnorm(rgrid3d_value_at_index(grid, i, j+1, k) - rgrid3d_value_at_index(grid, i, j-1, k))
	 + sqnorm(rgrid3d_value_at_index(grid, i+1, j, k) - rgrid3d_value_at_index(grid, i-1, j, k)));
  }
}

/*
 * Print the grid with both real and imaginary parts into file (ASCII format).
 *
 * grid = grid to be printed out (rgrid3d *).
 * out  = output file pointer (FILE *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_print(const rgrid3d *grid, FILE *out) {

  long i, j, k;

  for(i = 0; i < grid->nx; i++) {
    for(j = 0; j < grid->ny; j++) {
      for(k = 0; k < grid->nz; k++) {
        fprintf(out, "%16.8le   ", rgrid3d_value_at_index(grid, i, j, k));
	  }
      fprintf(out, "\n");
    }
    fprintf(out, "\n");
  }
}

/*
 * Perform Fast Fourier Transformation of a grid.
 *
 * grid = grid to be Fourier transformed (input/output) (rgrid3d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       Also no normalization is performed.
 *
 */

EXPORT void rgrid3d_fft(rgrid3d *grid) {

  if (!grid->plan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in rgrid3d_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    rgrid3d_fftw_alloc(grid);
  }
  rgrid3d_fftw(grid);
}

/*
 * Perform inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid3d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *       No normalization.
 *
 */

EXPORT void rgrid3d_inverse_fft(rgrid3d *grid) {

  if (!grid->iplan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in rgrid3d_inverse_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    rgrid3d_fftw_alloc(grid);
  }
  rgrid3d_fftw_inv(grid);
}
 
 
/*
 * Perform scaled inverse Fast Fourier Transformation of a grid.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid3d *).
 * c    = scaling factor (i.e. the output is multiplied by this constant) (double).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */
 
EXPORT void rgrid3d_scaled_inverse_fft(rgrid3d *grid, double c) {
   
  rgrid3d_inverse_fft(grid);
  rgrid3d_multiply(grid, c);  
}

/*
 * Perform inverse Fast Fourier Transformation of a grid scaled by FFT norm.
 *
 * grid = grid to be inverse Fourier transformed (input/output) (rgrid3d *).
 *
 * No return value.
 *
 * Note: The input grid is overwritten with the output.
 *
 */

EXPORT void rgrid3d_inverse_fft_norm(rgrid3d *grid) {

  rgrid3d_scaled_inverse_fft(grid, grid->fft_norm);
}
  

EXPORT void rgrid3d_fft_convolute(rgrid3d *gridc, const rgrid3d *grida, const rgrid3d *gridb) {
	if( gridc->value_outside == RGRID3D_PERIODIC_BOUNDARY ){
		rgrid3d_fft_periodic_convolute(gridc, grida, gridb) ;
	}else{
		rgrid3d_fft_nonperiodic_convolute(gridc, grida, gridb) ;
	}
}
/*
 * Convolute FFT transformed grids. To apply this on two grids (grida and gridb)
 * and place the result in gridc:
 * rgrid3d_fft(grida);
 * rgrid3d_fft(gridb);
 * rgrid3d_convolue(gridc, grida, gridb);
 * rgrid3d_inverse_fft(gridc);
 * gridc now contains the convolution of grida and gridb.
 *
 * grida = 1st grid to be convoluted (rgrid3d *).
 * gridb = 2nd grid to be convoluted (rgrid3d *).
 * gridc = output (rgrid3d *).
 *
 * No return value.
 *
 * Note: the input/output grids may be the same.
 *
 */

EXPORT void rgrid3d_fft_periodic_convolute(rgrid3d *gridc, const rgrid3d *grida, const rgrid3d *gridb) {

  long i, j, k, ij, ijnz, nx, ny, nz, nxy;
  double step = gridc->step, norm = step * step * step * grida->fft_norm;
  double complex *avalue, *bvalue, *cvalue;

  /* int f(r) g(r-r') d^3r' = iF[ F[f] F[g] ] = (step / N)^3 iFFT[ FFT[f] FFT[g] ] */
  
  nx = gridc->cint->nx;
  ny = gridc->cint->ny;
  nz = gridc->cint->nz;
  nxy = nx * ny;
  avalue = grida->cint->value;
  bvalue = gridb->cint->value;
  cvalue = gridc->cint->value;
#pragma omp parallel for firstprivate(nx,ny,nz,nxy,avalue,bvalue,cvalue,norm) private(i,j,ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++) {
      /* if odd */
      if ((i + j + k) & 1)
	cvalue[ijnz + k] = -norm * avalue[ijnz + k] * bvalue[ijnz + k];
      else
	cvalue[ijnz + k] = norm * avalue[ijnz + k] * bvalue[ijnz + k];
    }
  }
}

/* 
 * Convolution for nonperiodic grids. ("gridc = norm * grida * gridb").
 *
 * gridc = destination grid (rgrid3d *).
 * grida = 1st source grid (rgrid3d *).
 * gridb = 2nd source grid (rgrid3d *).
 *
 * No return value.
 *
 * Note: source and destination grids may be the same.
 * Note: function identical to rgrid3d_product except for a normalization factor
 *
 */

EXPORT void rgrid3d_fft_nonperiodic_convolute(rgrid3d *gridc, const rgrid3d *grida, const rgrid3d *gridb) {

  long ij, k, ijnz, nxy = gridc->nx * gridc->ny, nz = gridc->nz;
  double *avalue = grida->value;
  double *bvalue = gridb->value;
  double *cvalue = gridc->value;
  double step = gridc->step;
  if (!gridc->plan) {
    if (grid_threads() == 0) {
      fprintf(stderr, "libgrid: Error in rgrid3d_fft(). Function grid_threads_init() must be called before this function in parallel programs.\n");
      abort();
    }
    rgrid3d_fftw_alloc(gridc);
  }
  double norm = step * step * step * gridc->fft_norm;

#pragma omp parallel for firstprivate(nxy,nz,avalue,bvalue,cvalue,norm) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = norm * avalue[ijnz + k] * bvalue[ijnz + k];
  }
}


/* Boundary condition routines */

EXPORT double rgrid3d_value_outside_constantdirichlet(const rgrid3d *grid, long i, long j, long k) {

  return *((double *) grid->outside_params_ptr);
}

/*
 * The symmetry point are i=0 and i=nx-1 for consistency with the FFT plan FFTW_REDFT00. 
 * If one wants to use REDFT01 the symmetry points are i=-0.5 and i=nx-0.5 
 */
EXPORT double rgrid3d_value_outside_neumann(const rgrid3d *grid, long i, long j, long k) {

  long nx = grid->nx, ny = grid->ny, nz = grid->nz, nd;

  nd = nx * 2;
  if (i < 0) i  = -i ;
  if (i >= nd) i %= nd;
  if (i >= nx) i  = nd - i ;
  
  nd = ny * 2;
  if (j < 0) j  = -j ;
  if (j >= nd) j %= nd;
  if (j >= ny) j  = nd - j ;
  
  nd = nz * 2;
  if (k < 0) k  = -k ;
  if (k >= nd) k %= nd;
  if (k >= nz) k  = nd - k ;

  return grid->value[(i*ny + j)*nz + k];
}

EXPORT double rgrid3d_value_outside_periodic(const rgrid3d *grid, long i, long j, long k) {

  long nx = grid->nx, ny = grid->ny, nz = grid->nz;
  
  i %= nx;
  if (i < 0) i = nx + i;
  j %= ny;
  if (j < 0) j = ny + j;
  k %= nz;
  if (k < 0) k = nz + k;
  
  return grid->value[(i*ny + j)*nz + k];  
}

EXPORT double rgrid3d_value_outside_vortex(const rgrid3d *grid, long i, long j, long k) {

  long nx = grid->nx, ny = grid->ny, nz = grid->nz;

  i %= nx;
  if (i < 0) i = nx + i;
  j %= ny;
  if (j < 0) j = ny + j;
  k %= nz;
  if (k < 0) k = nz + k;
  
  return grid->value[(i*ny + j)*nz + k];  
}

/* End boundary condition routines */

/*
 * Access grid point at given index.
 *
 * grid = grid to be accessed (rgrid3d *).
 * i    = index along x (long).
 * j    = index along y (long).
 * k    = index along z (long).
 *
 * Returns grid value at index (i, j, k).
 *
 */

EXPORT inline double rgrid3d_value_at_index(const rgrid3d *grid, long i, long j, long k) {

  if (i < 0 || j < 0 || k < 0 || i >= grid->nx || j >= grid->ny || k >= grid->nz)
    return grid->value_outside(grid, i, j, k);
  return grid->value[(i*grid->ny + j)*grid->nz + k];
}

/*
 * Access grid point at given (x,y,z) point (with linear interpolation).
 *
 * grid = grid to be accessed (rgrid3d *).
 * x    = x value (double).
 * y    = y value (double).
 * z    = z value (double).
 *
 * Returns grid value at (x,y,z).
 *
 * (the most positive elements are missing - rolled over to negative; 
 *  periodic boundaries)
 *
 */

EXPORT inline double rgrid3d_value(const rgrid3d *grid, double x, double y, double z) {

  double f000, f100, f010, f001, f110, f101, f011, f111;
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
  f000 = rgrid3d_value_at_index(grid, i, j, k);
  f100 = rgrid3d_value_at_index(grid, i+1, j, k);
  f010 = rgrid3d_value_at_index(grid, i, j+1, k);
  f001 = rgrid3d_value_at_index(grid, i, j, k+1);
  f110 = rgrid3d_value_at_index(grid, i+1, j+1, k);
  f101 = rgrid3d_value_at_index(grid, i+1, j, k+1);
  f011 = rgrid3d_value_at_index(grid, i, j+1, k+1);
  f111 = rgrid3d_value_at_index(grid, i+1, j+1, k+1);
  
  omx = 1.0 - x;
  omy = 1.0 - y;
  omz = 1.0 - z;

  return omx * omy * omz * f000 + x * omy * omz * f100 + omx * y * omz * f010 + omx * omy * z * f001
    + x * y * omz * f110 + x * omy * z * f101 + omx * y * z * f011 + x * y * z * f111;
}
 
/*
 * Copy a real grid to a complex grid (to real part).
 *
 * dest   = destination grid (cgrid3d *).
 * source = source grid (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_real_to_complex_re(cgrid3d *dest, rgrid3d *source) {
  
  long ij, k, nz = source->nz, nxy = source->nx * source->ny, ijnz;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->nz = source->nz;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,nz,dst,src) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      dst[ijnz + k] = (double complex) src[ijnz + k];
  }
}

/*
 * Copy a real grid to a complex grid (to imaginary part).
 *
 * dest   = destination grid (cgrid3d *).
 * source = source grid (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_real_to_complex_im(cgrid3d *dest, rgrid3d *source) {
  
  long ij, k, nz = source->nz, nxy = source->nx * source->ny, ijnz;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->nz = source->nz;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,nz,dst,src) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      dst[ijnz + k] = I * (double complex) src[ijnz + k];
  }
}

/*
 * Add a real grid to a complex grid (to real part).
 *
 * dest   = destination grid (cgrid3d *).
 * source = source grid (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_add_real_to_complex_re(cgrid3d *dest, rgrid3d *source) {
  
  long ij, k, nz = source->nz, nxy = source->nx * source->ny, ijnz;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->nz = source->nz;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,nz,dst,src) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      dst[ijnz + k] += (double complex) src[ijnz + k];
  }
}

/*
 * Add a real grid to a complex grid (to real part).
 *
 * dest   = destination grid (cgrid3d *).
 * source = source grid (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_add_real_to_complex_im(cgrid3d *dest, rgrid3d *source) {
  
  long ij, k, nz = source->nz, nxy = source->nx * source->ny, ijnz;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->nz = source->nz;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,nz,dst,src) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      dst[ijnz + k] += I * (double complex) src[ijnz + k];
  }
}


/*
 * Product of a real grid with a complex grid.
 *
 * dest   = destination grid (cgrid3d *).
 * source = source grid (rgrid3d *).
 *
 * "dest(complex) = dest(complex) * source(real)"
 *
 * No return value.
 *
 */

EXPORT void grid3d_product_complex_with_real(cgrid3d *dest, rgrid3d *source) {
  
  long ij, k, nz = source->nz, nxy = source->nx * source->ny, ijnz;
  double *src = source->value;
  double complex *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->nz = source->nz;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,nz,dst,src) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      dst[ijnz + k] *= (double complex) src[ijnz + k];
  }
}


/*
 * Copy imaginary part of a complex grid to a real grid.
 *
 * dest   = destination grid (rgrid3d *).
 * source = source grid (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_complex_im_to_real(rgrid3d *dest, cgrid3d *source) {
  
  long ij, k, nz = source->nz, nxy = source->nx * source->ny, ijnz;
  double complex *src = source->value;
  double *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->nz = source->nz;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,nz,dst,src) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      dst[ijnz + k] = cimag(src[ijnz + k]);
  }
}

/*
 * Copy real part of a complex grid to a real grid.
 *
 * dest   = destination grid (rgrid3d *).
 * source = source grid (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_complex_re_to_real(rgrid3d *dest, cgrid3d *source) {
  
  long ij, k, nz = source->nz, nxy = source->nx * source->ny, ijnz;
  double complex *src = source->value;
  double *dst = dest->value;
  
  dest->nx = source->nx;
  dest->ny = source->ny;
  dest->nz = source->nz;
  dest->step = source->step;
  
#pragma omp parallel for firstprivate(nxy,nz,dst,src) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      dst[ijnz + k] = creal(src[ijnz + k]);
  }
}

/*
 * Extrapolate between two different grid sizes.
 *
 * dest = Destination grid (rgrid3d *; output).
 * src  = Source grid (rgrid3d *; input).
 *
 */

EXPORT void rgrid3d_extrapolate(rgrid3d *dest, rgrid3d *src) {

  long i, j, k, nx = dest->nx, ny = dest->ny, nz = dest->nz;
  double x0 = dest->x0, y0 = dest->y0, z0 = dest->z0;
  double step = dest->step, x, y, z;

  for (i = 0; i < nx; i++) {
    x = (i - nx/2) * step - x0;
    for (j = 0; j < ny; j++) {
      y = (j - ny/2) * step - y0;
      for (k = 0; k < nz; k++) {
	z = (k - nz/2) * step - z0;
	dest->value[i * ny * nz + j * nz + k] = rgrid3d_value(src, x, y, z);
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
 *
 */

double rgrid3d_value_rotate_z(void *arg, double x, double y, double z) {

  /* Unpack the values in arg */ 
  rgrid3d *grid = ((rotation *) arg)->rgrid;
  double sth = ((rotation *) arg)->sinth, cth = ((rotation *) arg)->costh, xp, yp;

  xp = -y * sth + x * cth; 
  yp =  y * cth + x * sth;

  return rgrid3d_value(grid, xp, yp, z);
}

/*
 * Rotate a grid by a given angle around the z-axis.
 *  cgrid3d *in : pointer with original grid.
 *  cgrid3d *out : pointer with rotated grid.
 *  double th: angle (radians) of the rotation.
 *
 *  The grid in and out CANNOT be the same.
 */
EXPORT void rgrid3d_rotate_z(rgrid3d *out, rgrid3d *in, double th) {

  rotation *r;
  if (in == out){
	  fprintf(stderr,"libgrid: in and out grids in rgrid3d_rotate_z must be different\n") ;
	  abort() ;
  }

  r = malloc(sizeof(rotation));
  r->rgrid = in;
  r->sinth = sin(-th);  // same direction of rotation as -wLz
  r->costh = cos(th);
  
  rgrid3d_map(out, rgrid3d_value_rotate_z, (void *) r);
  free(r);
}

/*
 * Get the largest value contained in a grid
 */

EXPORT double rgrid3d_max(rgrid3d *grid) {

  long nxyz = grid->nx * grid->ny * grid->nz;
  double *val = grid->value, max_val = val[0];
  long i;
  
#pragma omp parallel for firstprivate(nxyz, val) private(i) reduction(max: max_val) default(none) schedule(runtime)
  for(i=0 ; i< nxyz ; i++)
    if(val[i] > max_val) max_val = val[i];

  return max_val;
}

/*
 * Get the lowest value contained in a grid
 */

EXPORT double rgrid3d_min(rgrid3d *grid) {

  long nxyz = grid->nx * grid->ny * grid->nz;
  double *val = grid->value, min_val = val[0];
  long i;
  
#pragma omp parallel for firstprivate(nxyz, val) private(i) reduction(min: min_val) default(none) schedule(runtime)
  for(i=0 ; i< nxyz ; i++)
    if(val[i] < min_val) min_val = val[i];
  return min_val;
}

/*
 * Add random noise to grid.
 *
 * grid  = Grid where the noise will be added (cgrid3d *).
 * scale = Scaling for random numbers [-1,+1[ (double).
 *
 */

EXPORT void rgrid3d_random(rgrid3d *grid, double scale) {

  static int been_here = 0;
  long i;

  if(!been_here) {
    srand48(time(0));
    been_here = 1;
  }
  for (i = 0; i < grid->nx * grid->ny * grid->nz; i++)
    grid->value[i] += scale * 2.0 * (drand48() - 0.5);
}
