/*
 * Test program for wf1d routines.
 *
 * Propagate wavepacket in 3D.
 *
 * Try for example:
 * ./test_wf3d 1 128 0.1 200 8.0 8.0 8.0 0.25 0.25 0.25 -2.0 -2.0 -2.0 > test.dat
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grid/grid.h>
#include <omp.h>

double complex wavepacket(void *arg, double x, double y, double z);
double complex harmonic(void *arg, double x, double y, double z);
double complex dip(void *arg, double x, double y, double z);

typedef struct wparams_struct {
  double kx, ky, kz;
  double wx, wy, wz;
  double xc, yc, zc;
} wparams;

typedef struct pparams_struct {
  double kx, ky, kz;
  double delta;
} pparams;

double complex external_potential(void *arg, double x, double y, double z) {

  return harmonic(arg, x, y, z);
}

int main(int argc, char *argv[]) {

  long i, j, k, l, n, iterations, threads;
  double x, y, z, step, lx, threshold;
  double time_step;
  double complex time;
  
  wf3d *gwf = NULL;
  cgrid3d *potential = NULL;
  cgrid3d *workspace = NULL;
  cgrid3d *workspace2 = NULL;
  cgrid3d *sq_grad_pot = NULL;
  rgrid3d *rworkspace = NULL;
  
  pparams potential_params;
  wparams wp_params;
  
  /* parameters */
  if (argc != 14) {
    fprintf(stderr, "Usage: %s <threads> <points/axis> <time_step> <iterations> <kx> <ky> <kz> <wx> <wy> <wz> <xc> <yc> <zc>\n", argv[0]);
    return -1;
  }
  
  threads = atoi(argv[1]);
  n = atoi(argv[2]);
  time_step = atof(argv[3]);
  iterations = atol(argv[4]);
  wp_params.kx = atof(argv[5]);
  wp_params.ky = atof(argv[6]);
  wp_params.kz = atof(argv[7]);
  wp_params.wx = atof(argv[8]);
  wp_params.wy = atof(argv[9]);
  wp_params.wz = atof(argv[10]);
  wp_params.xc = atof(argv[11]);
  wp_params.yc = atof(argv[12]);
  wp_params.zc = atof(argv[13]);
  if(wp_params.wx == 0.0 || wp_params.wy == 0.0 || wp_params.wz == 0.0) {
    fprintf(stderr, "Width cannot be zero.\n");
    exit(1);
  }
  if(wp_params.kx == 0.0 || wp_params.ky == 0.0 || wp_params.kz == 0.0) {
    fprintf(stderr, "force constant cannot be zero.\n");
    exit(1);
  }
  
  step = 0.4 / (n / 16.0);
  lx = n * step;
  
  fprintf(stderr, "Grid (%ldX%ldX%ld)\n", n, n, n);
  
  /* potential parameters */
  potential_params.kx = lx * 2.0;
  potential_params.ky = lx * 2.0;
  potential_params.kz = lx * 2.0;
  potential_params.delta = 1;
  
  grid_threads_init(threads);
  
  /* allocate memory (mass = 1.0) */
  gwf = grid3d_wf_alloc(n, n, n, step, 1.0, WF3D_PERIODIC_BOUNDARY, WF3D_2ND_ORDER_PROPAGATOR);
  //gwf = grid3d_wf_alloc(n, n, n, step, 1.0, WF3D_NEUMANN_BOUNDARY, WF3D_2ND_ORDER_PROPAGATOR);
  potential = cgrid3d_alloc(n, n, n, step, CGRID3D_PERIODIC_BOUNDARY, 0);
  workspace = cgrid3d_alloc(n, n, n, step, CGRID3D_PERIODIC_BOUNDARY, 0);
  workspace2 = cgrid3d_alloc(n, n, n, step, CGRID3D_PERIODIC_BOUNDARY, 0);
  rworkspace = rgrid3d_alloc(n, n, n, step, RGRID3D_PERIODIC_BOUNDARY, 0);
  sq_grad_pot = cgrid3d_alloc(n, n, n, step, CGRID3D_PERIODIC_BOUNDARY, 0);
  
  /* initialize wave function */
  grid3d_wf_map(gwf, wavepacket, &wp_params);
  grid3d_wf_normalize(gwf);
  
  /* map potential */
  cgrid3d_smooth_map(potential, external_potential, &potential_params, 1);
  
  /* solve */
  time = time_step;
  for(l = 0; l < iterations; l++) {
    grid3d_wf_square_of_potential_gradient(sq_grad_pot, potential, workspace, workspace2);
    grid3d_wf_propagate(gwf, potential, sq_grad_pot, time, workspace);
#if 1
    grid3d_wf_density(gwf, rworkspace);
    for (i = 0; i < n; i++) {
      x = (i - n/2) * step;
      printf("%le %le %le\n", time_step * (double) l, x, rgrid3d_value_at_index(rworkspace, i, n/2, n/2));
    }
    printf("\n");
    fflush(stdout);
#endif
  }

  /* release resources */
  grid3d_wf_free(gwf);
  cgrid3d_free(workspace);
  cgrid3d_free(workspace2);
  rgrid3d_free(rworkspace);
  cgrid3d_free(sq_grad_pot);
  cgrid3d_free(potential);
  
  return 0;
}

double complex wavepacket(void *arg, double x, double y, double z) {

  double kx = ((wparams *) arg)->kx;
  double ky = ((wparams *) arg)->ky;
  double kz = ((wparams *) arg)->kz;
  double wx = ((wparams *) arg)->wx;
  double wy = ((wparams *) arg)->wy;
  double wz = ((wparams *) arg)->wz;
  double xc = ((wparams *) arg)->xc;
  double yc = ((wparams *) arg)->yc;
  double zc = ((wparams *) arg)->zc;
  
  x -= xc;
  y -= yc;
  z -= zc;

  return cexp(-x/wx*x/wx + I * kx * x
	      -y/wy*y/wy + I * ky * y
	      -z/wz*z/wz + I * kz * z);
}

double complex harmonic(void *arg, double x, double y, double z) {

  pparams params = *((pparams *) arg);

  return 0.5 * params.delta * (
			       params.kx * params.kx * x * x
			       + params.ky * params.ky * y * y
			       + params.kz * params.kz * z * z
			       );
    
}

double complex dip(void *arg, double x, double y, double z) {

  double pot;
  pparams params = *((pparams *) arg);

  x = 2.0 * M_PI * x / (2.0 * params.kx); 
  y = 2.0 * M_PI * y / (2.0 * params.ky); 
  z = 2.0 * M_PI * z / (2.0 * params.kz); 

  pot = params.delta * 
    (1.0 - cos(x) * cos(x))
    * (1.0 - cos(y) * cos(y))
    * (1.0 - cos(z) * cos(z));

  return pot;
}
