/*
 * Test program for wf2d routines.
 *
 * Propagate wavepacket in 2D.
 *
 * Try for example:
 * ./test_wf2d 1 128 0.1 200 8.0 8.0 0.25 0.25 -2.0 -2.0 > test.dat
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grid/grid.h>
#include <omp.h>

double complex wavepacket(void *arg, double x, double y);
double complex harmonic(void *arg, double x, double y);
double complex dip(void *arg, double x, double y);

typedef struct wparams_struct {
  double kx, ky;
  double wx, wy;
  double xc, yc;
} wparams;

typedef struct pparams_struct {
  double kx, ky;
  double delta;
} pparams;

double complex external_potential(void *arg, double x, double y) {

  return harmonic(arg, x, y);
}

int main(int argc, char *argv[]) {

  long i, j, k, l, n, iterations, threads;
  double x, y, z, step, lx, threshold;
  double time_step;
  double complex time;
  
  wf2d *gwf = NULL;
  cgrid2d *potential = NULL;
  cgrid2d *cworkspace = NULL;
  rgrid2d *workspace = NULL;
  cgrid2d *sq_grad_pot = NULL;
  
  pparams potential_params;
  wparams wp_params;
  
  /* parameters */
  if (argc != 11) {
    fprintf(stderr, "Usage: %s <threads> <points/axis> <time_step> <iterations> <kx> <ky> <wx> <wy> <xc> <yc>\n", argv[0]);
    return -1;
  }
  
  threads = atoi(argv[1]);
  n = atoi(argv[2]);
  time_step = atof(argv[3]);
  iterations = atol(argv[4]);
  wp_params.kx = atof(argv[5]);
  wp_params.ky = atof(argv[6]);
  wp_params.wx = atof(argv[7]);
  wp_params.wy = atof(argv[8]);
  wp_params.xc = atof(argv[9]);
  wp_params.yc = atof(argv[10]);
  if(wp_params.wx == 0.0 || wp_params.wy == 0.0) {
    fprintf(stderr, "Width cannot be zero.\n");
    exit(1);
  }
  if(wp_params.kx == 0.0 || wp_params.ky == 0.0) {
    fprintf(stderr, "force constant cannot be zero.\n");
    exit(1);
  }
  
  step = 0.4 / (n / 16.0);
  lx = n * step;
  
  fprintf(stderr, "Grid (%ldX%ld)\n", n, n);
  
  /* potential parameters */
  potential_params.kx = lx * 2.0;
  potential_params.ky = lx * 2.0;
  potential_params.delta = 1;
  
  grid_threads_init(threads);
  
  /* allocate memory (mass = 1.0) */
  //gwf = grid2d_wf_alloc(n, n, step, 1.0, WF2D_PERIODIC_BOUNDARY, WF2D_2ND_ORDER_PROPAGATOR);
  gwf = grid2d_wf_alloc(n, n, step, 1.0, WF2D_NEUMANN_BOUNDARY, WF2D_2ND_ORDER_PROPAGATOR);
  potential = cgrid2d_alloc(n, n, step, CGRID2D_PERIODIC_BOUNDARY, 0);
  workspace = rgrid2d_alloc(n, n, step, RGRID2D_PERIODIC_BOUNDARY, 0);
  cworkspace = cgrid2d_alloc(n, n, step, CGRID2D_PERIODIC_BOUNDARY, 0);
  sq_grad_pot = cgrid2d_alloc(n, n, step, CGRID2D_PERIODIC_BOUNDARY, 0);
  
  /* initialize wave function */
  grid2d_wf_map(gwf, wavepacket, &wp_params);
  grid2d_wf_normalize(gwf);
  
  /* map potential */
  cgrid2d_smooth_map(potential, external_potential, &potential_params, 1);
  
  /* solve */
  time = time_step;
  for(l = 0; l < iterations; l++) {
    grid2d_wf_square_of_potential_gradient(sq_grad_pot, potential, cworkspace);
    grid2d_wf_propagate(gwf, potential, sq_grad_pot, time, cworkspace);
#if 1
    grid2d_wf_density(gwf, workspace);
    for (i = 0; i < n; i++) {
      x = (i - n/2) * step;
      printf("%le %le %le\n", time_step * (double) l, x, workspace->value[i]);
    }
    printf("\n");
    fflush(stdout);
#endif
  }

  /* release resources */
  grid2d_wf_free(gwf);
  cgrid2d_free(cworkspace);
  rgrid2d_free(workspace);
  cgrid2d_free(sq_grad_pot);
  cgrid2d_free(potential);
  
  return 0;
}

double complex wavepacket(void *arg, double x, double y) {

  double kx = ((wparams *) arg)->kx;
  double ky = ((wparams *) arg)->ky;
  double wx = ((wparams *) arg)->wx;
  double wy = ((wparams *) arg)->wy;
  double xc = ((wparams *) arg)->xc;
  double yc = ((wparams *) arg)->yc;
  
  x -= xc;
  y -= yc;

  return cexp(-x/wx*x/wx + I * kx * x
	      -y/wy*y/wy + I * ky * y);
}

double complex harmonic(void *arg, double x, double y) {

  pparams params = *((pparams *) arg);

  return 0.5 * params.delta * (
			       params.kx * params.kx * x * x
			       + params.ky * params.ky * y * y
			       );
    
}

double complex dip(void *arg, double x, double y) {

  double pot;
  pparams params = *((pparams *) arg);

  x = 2.0 * M_PI * x / (2.0 * params.kx); 
  y = 2.0 * M_PI * y / (2.0 * params.ky); 

  pot = params.delta * 
    (1.0 - cos(x) * cos(x))
    * (1.0 - cos(y) * cos(y));

  return pot;
}
