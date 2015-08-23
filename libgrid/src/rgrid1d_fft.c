/*
 * Interface to FFTW (real).
 * No user callable routines.
 *
 */

#include "grid.h"

static double *temp = NULL;
static long temp_len = 0;
static double complex *temp2 = NULL;

/*
 * Allocate FFT buffers for a given grid.
 * This must be called before FFT can be carried out
 * (called automatically, users don't need to call this)
 * 
 * grid = 1D grid for which the FFT allocation is to be done (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_fftw_alloc(rgrid1d *grid) {

  double *plan_temp;
  long s;

  s = 2 * (grid->nx/2 + 1);
  if(s > temp_len) {
    if(temp) free(temp);
    if(!(temp = fftw_malloc(sizeof(double) * s))) {
      fprintf(stderr, "libgrid: Out of memory in rgrid3d_fft().\n");
      return;
    }
    temp_len = s;
    temp2 = (double complex *) temp;
  }
  if(!(plan_temp = fftw_malloc(sizeof(double) * s))) {
    fprintf(stderr, "libgrid: Out of memory in rgrid3d_fft().\n");
    return;
  }

  memcpy(plan_temp, grid->value, sizeof(double) * s);
 
  fftw_plan_with_nthreads(grid_threads());

  grid->plan = 
    fftw_plan_dft_r2c_1d(grid->nx, grid->value, temp2, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);

  grid->iplan = 
    fftw_plan_dft_c2r_1d(grid->nx, grid->cint->value, temp, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);

  memcpy(grid->value, plan_temp, sizeof(double) * s);

  fftw_free(plan_temp);
}

/*
 * Free FFT buffers. Used only internally.
 *
 * grid = 1D grid for which the FFT buffers are to be freed (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_fftw_free(rgrid1d *grid) {

  if(grid->plan) fftw_destroy_plan(grid->plan);
  if(grid->iplan) fftw_destroy_plan(grid->iplan);
}

/*
 * Forward FFT using FFTW. Used only internally.
 *
 * grid = 1D grid to be transformed (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_fftw(rgrid1d *grid) {

  fftw_execute(grid->plan);
  memcpy(grid->cint->value, temp2, sizeof(double complex) * grid->cint->nx);
}

/*
 * Backward FFT using FFTW. Used only internally.
 *
 * grid = 1D grid to be transformed (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void rgrid1d_fftw_inv(rgrid1d *grid) {
  
  fftw_execute(grid->iplan);
  memcpy(grid->value, temp, sizeof(double) * grid->nx);
}
