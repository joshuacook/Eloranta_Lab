/*
 * Interface to FFTW.
 *
 * No user callable routines.
 *
 */

#include "grid.h"

/*
 * Allocate FFT buffers for a given grid.
 * This must be called before FFT can be carried out
 * (called automatically, users don't need to call this)
 * 
 * grid = 3D grid for which the FFT allocation is to be done (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_fftw_alloc(cgrid3d *grid) {

  double complex *temp;
  long size = sizeof(double complex) * grid->nx * grid->ny * grid->nz;
  int n[3] = {grid->nx, grid->ny, grid->nz};
  fftw_r2r_kind rfk[3], rbk[3], ifk[3], ibk[3];

  temp = fftw_malloc(size);
  memcpy(temp, grid->value, size);
  /* NOTE: To see which boundary condition applies, see if grid->implan (or grid->implan) is NULL */
  /* If grid->implan is NULL, this is VORTEX boundary; if not, it is standard PERIODIC boundary */
  if(grid->value_outside == CGRID3D_PERIODIC_BOUNDARY) {
    fftw_plan_with_nthreads(grid_threads());
    grid->plan = fftw_plan_dft_3d(grid->nx, grid->ny, grid->nz,
				  grid->value, grid->value, FFTW_FORWARD, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
    fftw_plan_with_nthreads(grid_threads());
    grid->iplan = fftw_plan_dft_3d(grid->nx, grid->ny, grid->nz,
				   grid->value, grid->value, FFTW_BACKWARD, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
    grid->implan = grid->iimplan = NULL;
    grid->fft_norm = 1.0 / (grid->nx * grid->ny * grid->nz);
    memcpy(grid->value, temp, size);
    fftw_free(temp);
    return;
  } else if (grid->value_outside == CGRID3D_NEUMANN_BOUNDARY) {
    rfk[0] = FFTW_REDFT10; rfk[1] = FFTW_REDFT10; rfk[2] = FFTW_REDFT10;
    rbk[0] = FFTW_REDFT01; rbk[1] = FFTW_REDFT01; rbk[2] = FFTW_REDFT01;
    ifk[0] = FFTW_REDFT10; ifk[1] = FFTW_REDFT10; ifk[2] = FFTW_REDFT10;
    ibk[0] = FFTW_REDFT01; ibk[1] = FFTW_REDFT01; ibk[2] = FFTW_REDFT01;
  } else if (grid->value_outside == CGRID3D_VORTEX_X_BOUNDARY) {
    rfk[0] = FFTW_REDFT10; rfk[1] = FFTW_RODFT10; rfk[2] = FFTW_REDFT10;
    rbk[0] = FFTW_REDFT01; rbk[1] = FFTW_RODFT01; rbk[2] = FFTW_REDFT01;
    ifk[0] = FFTW_REDFT10; ifk[1] = FFTW_REDFT10; ifk[2] = FFTW_RODFT10;
    ibk[0] = FFTW_REDFT01; ibk[1] = FFTW_REDFT01; ibk[2] = FFTW_RODFT01;
  } else if (grid->value_outside == CGRID3D_VORTEX_Y_BOUNDARY) {
    rfk[0] = FFTW_REDFT10; rfk[1] = FFTW_REDFT10; rfk[2] = FFTW_RODFT10;
    rbk[0] = FFTW_REDFT01; rbk[1] = FFTW_REDFT01; rbk[2] = FFTW_RODFT01;
    ifk[0] = FFTW_RODFT10; ifk[1] = FFTW_REDFT10; ifk[2] = FFTW_REDFT10;
    ibk[0] = FFTW_RODFT01; ibk[1] = FFTW_REDFT01; ibk[2] = FFTW_REDFT01;
  } else if (grid->value_outside == CGRID3D_VORTEX_Z_BOUNDARY) {
    rfk[0] = FFTW_RODFT10; rfk[1] = FFTW_REDFT10; rfk[2] = FFTW_REDFT10;
    rbk[0] = FFTW_RODFT01; rbk[1] = FFTW_REDFT01; rbk[2] = FFTW_REDFT01;
    ifk[0] = FFTW_REDFT10; ifk[1] = FFTW_RODFT10; ifk[2] = FFTW_REDFT10;
    ibk[0] = FFTW_REDFT01; ibk[1] = FFTW_RODFT01; ibk[2] = FFTW_REDFT01;
  } else {
    fprintf(stderr, "libgrid: Incompatible boundary condition for FFT.\n");
    exit(1);
  }
  fftw_plan_with_nthreads(grid_threads());
  grid->plan = fftw_plan_many_r2r(3, n, 1, (double *) grid->value, n, 2, 0, (double *) grid->value, n, 2, 0, rfk, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
  fftw_plan_with_nthreads(grid_threads());
  grid->iplan = fftw_plan_many_r2r(3, n, 1, (double *) grid->value, n, 2, 0, (double *) grid->value, n, 2, 0, rbk, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
  fftw_plan_with_nthreads(grid_threads());
  grid->implan = fftw_plan_many_r2r(3, n, 1, ((double *) grid->value) + 1, n, 2, 0, ((double *) grid->value) + 1, n, 2, 0, ifk, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
  fftw_plan_with_nthreads(grid_threads());
  grid->iimplan = fftw_plan_many_r2r(3, n, 1, ((double *) grid->value) + 1, n, 2, 0, ((double *) grid->value) + 1, n, 2, 0, ibk, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);
  grid->fft_norm = 1.0 / (2.0 * grid->nx * 2.0 * grid->ny * 2.0 * grid->nz);
  memcpy(grid->value, temp, size);
  fftw_free(temp);
}

/*
 * Free FFT buffers. Used only internally.
 *
 * grid = 3D grid for which the FFT buffers are to be freed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_fftw_free(cgrid3d *grid) {

  if(grid->plan) fftw_destroy_plan(grid->plan);
  if(grid->iplan) fftw_destroy_plan(grid->iplan);
  if(grid->implan) fftw_destroy_plan(grid->implan);
  if(grid->iimplan) fftw_destroy_plan(grid->iimplan);
  grid->plan = grid->iplan = grid->implan = grid->iimplan = NULL;
}

/*
 * Forward FFT using FFTW. Used only internally.
 *
 * grid = 3D grid to be transformed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_fftw(cgrid3d *grid) {

  fftw_execute(grid->plan);
  if(grid->implan) fftw_execute(grid->implan);
}

/*
 * Backward FFT using FFTW. Used only internally.
 *
 * grid = 3D grid to be transformed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void cgrid3d_fftw_inv(cgrid3d *grid) {

  fftw_execute(grid->iplan);
  if(grid->iimplan) fftw_execute(grid->iimplan);
}
