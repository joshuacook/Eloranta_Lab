/*
 * Interface to FFTW. Real version.
 *
 * No user callable routines.
 *
 * TODO: Somehow this only works for 3-D arrays where all dimensions
 * are multiples of two! (?) For example, 256, 256, 256 will work
 * fine but 255, 255, 255 does not. This is for periodic.
 *
 */

#include "grid.h"


static double *temp = NULL;
static long temp_len = 0;
static double complex *temp2 = NULL;

/*
 * These functions select between periodic and non-periodic functions.
 * The periodic FFT does r2c transform and allocate different
 * pointers for the in/out (use of cint interface).
 * The non-periodic FFT does r2r in place.
 *
 */

EXPORT void rgrid3d_fftw_alloc(rgrid3d *grid) {
	if(grid->value_outside == RGRID3D_PERIODIC_BOUNDARY)
		rgrid3d_fftw_periodic_alloc(grid) ;
	else
		rgrid3d_fftw_nonperiodic_alloc(grid) ;
}

EXPORT void rgrid3d_fftw(rgrid3d *grid) {
	if(grid->value_outside == RGRID3D_PERIODIC_BOUNDARY)
		rgrid3d_fftw_periodic(grid) ;
	else
		fftw_execute(grid->plan) ;
}

EXPORT void rgrid3d_fftw_inv(rgrid3d *grid) {
	if(grid->value_outside == RGRID3D_PERIODIC_BOUNDARY)
		rgrid3d_fftw_periodic_inv(grid) ;
	else
		fftw_execute(grid->iplan) ;
}




/*
 * Allocate FFT buffers for a given grid.
 * This must be called before FFT can be carried out
 * (called automatically, users don't need to call this)
 * 
 * grid = 3D grid for which the FFT allocation is to be done (rgrid3d *).
 *
 * Notes: 
 *  - This uses c2r transformation since it is much faster than r2r.
 *
 * No return value.
 *
 */
EXPORT void rgrid3d_fftw_periodic_alloc(rgrid3d *grid) {

  long s;
  double *plan_temp;

  s = 2 * grid->nx * grid->ny * (grid->nz/2 + 1);

  if(!(plan_temp = fftw_malloc(sizeof(double) * s))) {
    fprintf(stderr, "libgrid: Out of memory in rgrid3d_fft().\n");
    return;
  }

  memcpy(plan_temp, grid->value, sizeof(double) * s);
 
  if(s > temp_len) {
    if(temp) free(temp);
    if(!(temp = fftw_malloc(sizeof(double) * s))) {
      fprintf(stderr, "libgrid: Out of memory in rgrid3d_fft().\n");
      return;
    }
    temp_len = s;
    temp2 = (double complex *) temp;
  }
  fftw_plan_with_nthreads(grid_threads());
  grid->plan = 
    fftw_plan_dft_r2c_3d(grid->nx, grid->ny, grid->nz, grid->value, temp2, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
  fftw_plan_with_nthreads(grid_threads());
  grid->iplan = 
    fftw_plan_dft_c2r_3d(grid->nx, grid->ny, grid->nz, grid->cint->value, temp, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
  grid->fft_norm = 1.0 / (grid->nx * grid->ny * grid->nz);
  memcpy(grid->value, plan_temp, sizeof(double) * s);
  free(plan_temp);
}

/* Non-periodic version of alloc */
/* The use of REDFT00 is prefered over REDFT01 because 
 * REDFT00 starts at x=0 and REDFT01 at x=0.5*step. This makes
 * REDFT00 a better choice to reproduce the
 * coarse-graining and Lennard-Jones kernels. 
 * With REDFT00 their integral (or their k=0 FT) is closer
 * to 1 and -b respectively.
 *
 * Note the normalization depends on the plan used:
 *  REDFT00 -> 1/(2*nx-1)
 *  REDFT01 -> 1/(2*nx)
 *
 */
EXPORT void rgrid3d_fftw_nonperiodic_alloc(rgrid3d *grid) {
  long s = grid-> nx * grid->ny * grid-> nz ;
  double *temp ;
  temp = malloc( sizeof(double)*s ) ;
  memcpy(temp, grid->value, sizeof(double)*s ) ;

  fftw_plan_with_nthreads(grid_threads());
  grid->plan = fftw_plan_r2r_3d(grid->nx, grid->ny, grid->nz,
		  		grid->value, grid->value,
				FFTW_REDFT10 , FFTW_REDFT10 , FFTW_REDFT10 ,
				GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
  fftw_plan_with_nthreads(grid_threads());
  grid->iplan = fftw_plan_r2r_3d(grid->nx, grid->ny, grid->nz,
		  		grid->value, grid->value,
				FFTW_REDFT01 , FFTW_REDFT01 , FFTW_REDFT01 ,
				GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
  grid->fft_norm = 1.0 / (2 * (grid->nx ) * 2 * (grid->ny ) * 2 * (grid->nz ));
  
  memcpy(grid->value, temp, sizeof(double)*s );
  free(temp);
}



/*
 * Free FFT buffers. Used only internally.
 *
 * grid = 3D grid for which the FFT buffers are to be freed (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_fftw_free(rgrid3d *grid) {

  if(grid->plan) fftw_destroy_plan(grid->plan);
  if(grid->iplan) fftw_destroy_plan(grid->iplan);
}

/*
 * Forward FFT using FFTW. Used only internally.
 *
 * grid = 3D grid to be transformed (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_fftw_periodic(rgrid3d *grid) {

  long i, nx = grid->cint->nx, nyz = grid->cint->ny * grid->cint->nz, bytes = sizeof(double complex) * nyz;
  double complex *value = grid->cint->value;

  fftw_execute(grid->plan);
  //memcpy(grid->cint->value, temp2, sizeof(double complex) * grid->cint->nx * grid->cint->ny * grid->cint->nz);
#pragma omp parallel for firstprivate(nx, nyz, bytes, value, temp2) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    memcpy(&value[i * nyz], &temp2[i * nyz], bytes);
}

/*
 * Backward FFT using FFTW. Used only internally.
 *
 * grid = 3D grid to be transformed (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid3d_fftw_periodic_inv(rgrid3d *grid) {

  long i, nx = grid->nx, nyz = grid->ny * grid->nz, bytes = sizeof(double) * nyz;
  double *value = grid->value;    
  
  fftw_execute(grid->iplan);
  
  // memcpy(grid->value, temp, sizeof(double) * grid->nx * grid->ny * grid->nz);
#pragma omp parallel for firstprivate(nx, nyz, bytes, value, temp) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    memcpy(&value[i * nyz], &temp[i * nyz], bytes);
}
