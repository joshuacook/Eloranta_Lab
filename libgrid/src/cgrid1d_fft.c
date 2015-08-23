/*
 * Interface to FFTW.
 * No user callable routines.
 *
 */

#include "grid.h"

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

EXPORT void cgrid1d_fftw_alloc(cgrid1d *grid) {

  double complex *temp;
  long size = sizeof(double complex) * grid->nx;

  temp = fftw_malloc(size);
  memcpy(temp, grid->value, size);
  fftw_plan_with_nthreads(grid_threads());
  grid->plan = fftw_plan_dft_1d(grid->nx, grid->value, grid->value, FFTW_FORWARD, GRID_FFTW_PLAN);
  fftw_plan_with_nthreads(grid_threads());
  grid->iplan = fftw_plan_dft_1d(grid->nx, grid->value, grid->value, FFTW_BACKWARD, GRID_FFTW_PLAN);
  memcpy(grid->value, temp, size);
  fftw_free(temp);
}

/*
 * Free FFT buffers. Used only internally.
 *
 * grid = 1D grid for which the FFT buffers are to be freed (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_fftw_free(cgrid1d *grid) {

  if(grid->plan) fftw_destroy_plan(grid->plan);
  if(grid->iplan) fftw_destroy_plan(grid->iplan);
}

/*
 * Forward FFT using FFTW. Used only internally.
 *
 * grid = 1D grid to be transformed (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_fftw(cgrid1d *grid) {

  fftw_execute(grid->plan);
}

/*
 * Backward FFT using FFTW. Used only internally.
 *
 * grid = 1D grid to be transformed (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid1d_fftw_inv(cgrid1d *grid) {
  
  fftw_execute(grid->iplan);
}
