#include <stdio.h>
#include <math.h>
#include <grid/grid.h>
#include <grid/au.h>
#include <dft/dft.h>
#include <dft/ot.h>

#define NX 256
#define NY 256
#define NZ 256
#define STEP 1.0

main() {
  
  rgrid3d *gauss;
  wf3d *wf;
  double inv_width = 1.0 / 1.0;

  grid_threads_init(16);
  gauss = rgrid3d_alloc(NX, NY, NZ, STEP, RGRID3D_VORTEX_BOUNDARY, 0);
  wf = grid3d_wf_alloc(NX, NY, NZ, STEP, 100.0, WF3D_VORTEX_Z_BOUNDARY, WF3D_2ND_ORDER_PROPAGATOR);
  
  rgrid3d_map(gauss, dft_common_gaussian, &inv_width);
  dft_driver_write_density(gauss, "ri");
  rgrid3d_power(gauss, gauss, 0.5);
  grid3d_real_to_complex_re(wf->grid, gauss);
  grid3d_wf_propagate_kinetic_fft(wf, 1.0 / GRID_AUTOFS);
  grid3d_wf_density(wf, gauss);
  dft_driver_write_density(gauss, "co");
}

