#include <stdio.h>
#include <math.h>
#include <grid/grid.h>
#include <dft/dft.h>
#include <dft/ot.h>

#define NZ 1024
#define NR 1024
#define STEP 0.1

main() {
  
  rgrid2d *gauss;
  double inv_width = 1.0 / 10.0;

  grid_threads_init(1);
  gauss = rgrid2d_alloc(NZ, NR, STEP, RGRID2D_NEUMANN_BOUNDARY, 0);
  
  rgrid2d_map_cyl(gauss, dft_common_gaussian_2d, &inv_width);
  rgrid2d_fft_cylindrical(gauss);
  rgrid2d_inverse_fft_cylindrical(gauss);
  dft_driver_write_density_2d(gauss, "o");
}

