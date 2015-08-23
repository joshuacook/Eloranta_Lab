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
 * grid = 2D grid for which the FFT allocation is to be done (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_fftw_alloc(cgrid2d *grid) {

  double complex *temp;
  long size = sizeof(double complex) * grid->nx * grid->ny;

  temp = fftw_malloc(size);
  memcpy(temp, grid->value, size);
  fftw_plan_with_nthreads(grid_threads());
  grid->plan = fftw_plan_dft_2d(grid->nx, grid->ny,
				grid->value, grid->value, FFTW_FORWARD, GRID_FFTW_PLAN);
  fftw_plan_with_nthreads(grid_threads());
  grid->iplan = fftw_plan_dft_2d(grid->nx, grid->ny,
				 grid->value, grid->value, FFTW_BACKWARD, GRID_FFTW_PLAN);
  memcpy(grid->value, temp, size);
  fftw_free(temp);
}

/*
 * Free FFT buffers. Used only internally.
 *
 * grid = 2D grid for which the FFT buffers are to be freed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_fftw_free(cgrid2d *grid) {

  if(grid->plan) fftw_destroy_plan(grid->plan);
  if(grid->iplan) fftw_destroy_plan(grid->iplan);
}

/*
 * Forward FFT using FFTW. Used only internally.
 *
 * grid = 2D grid to be transformed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_fftw(cgrid2d *grid) {

  fftw_execute(grid->plan);
}

/*
 * Backward FFT using FFTW. Used only internally.
 *
 * grid = 2D grid to be transformed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid2d_fftw_inv(cgrid2d *grid) {

  fftw_execute(grid->iplan);
}
