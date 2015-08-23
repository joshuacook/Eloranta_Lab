/*
 * Test program for wf1d routines.
 *
 * Propagate wavepacket in 1D.
 * 
 * Try for example:
 * ./test_wf1d 1 5000 0.01 200 8.0 0.25 -2.0 > test.dat
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grid/grid.h>
#include <omp.h>

double complex wavepacket(void *arg, double x);
double complex harmonic(void *arg, double x);
double complex dip(void *arg, double x);

typedef struct wparams_struct {
  double kx;
  double wx;
  double xc;
} wparams;

typedef struct pparams_struct {
  double kx;
  double delta;
} pparams;

double complex external_potential(void *arg, double x) {

  return harmonic(arg, x);
}

int main(int argc, char *argv[]) {

  long i, j, k, l, n, iterations, threads;
  double x, step, lx, threshold;
  double time_step;
  double complex time;
  
  wf1d *gwf = NULL;
  cgrid1d *potential = NULL;
  cgrid1d *cworkspace = NULL;
  rgrid1d *workspace = NULL;
  cgrid1d *sq_grad_pot = NULL;
  
  pparams potential_params;
  wparams wp_params;
  
  /* parameters */
  if (argc != 8) {
    fprintf(stderr, "Usage: %s <threads> <points/axis> <time_step> <iterations> <kx> <wx> <xc>\n", argv[0]);
    return -1;
  }
  
  threads = atoi(argv[1]);
  n = atoi(argv[2]);
  time_step = atof(argv[3]);
  iterations = atol(argv[4]);
  wp_params.kx = atof(argv[5]);
  wp_params.wx = atof(argv[6]);
  wp_params.xc = atof(argv[7]);
  if(wp_params.wx == 0.0) {
    fprintf(stderr, "Width cannot be zero.\n");
    exit(1);
  }
  if(wp_params.kx == 0.0) {
    fprintf(stderr, "force constant cannot be zero.\n");
    exit(1);
  }
  
  step = 0.4 / (n / 16.0);
  lx = n * step;
  
  fprintf(stderr, "Grid (%ld)\n", n);
  
  /* potential parameters */
  potential_params.kx = lx * 2.0;
  potential_params.delta = 1;
  
  grid_threads_init(threads);
  
  /* allocate memory (mass = 1.0) */
  gwf = grid1d_wf_alloc(n, step, 1.0, WF1D_PERIODIC_BOUNDARY, WF1D_2ND_ORDER_PROPAGATOR);
  potential = cgrid1d_alloc(n, step, CGRID1D_PERIODIC_BOUNDARY, 0);
  cworkspace = cgrid1d_alloc(n, step, CGRID1D_PERIODIC_BOUNDARY, 0);
  workspace = rgrid1d_alloc(n, step, RGRID1D_PERIODIC_BOUNDARY, 0);
  sq_grad_pot = cgrid1d_alloc(n, step, CGRID1D_PERIODIC_BOUNDARY, 0);
  
  /* initialize wave function */
  grid1d_wf_map(gwf, wavepacket, &wp_params);
  grid1d_wf_normalize(gwf);
  
  /* map potential */
  cgrid1d_smooth_map(potential, external_potential, &potential_params, 1);
  
  /* solve */
  time = time_step;
  for(l = 0; l < iterations; l++) {
    grid1d_wf_square_of_potential_gradient(sq_grad_pot, potential);
    grid1d_wf_propagate(gwf, potential, sq_grad_pot, time, cworkspace);
#if 1
    grid1d_wf_density(gwf, workspace);
    for (i = 0; i < n; i++) {
      x = (i - n/2) * step;
      printf("%le %le %le\n", time_step * (double) l, x, workspace->value[i]);
    }
    printf("\n");
    fflush(stdout);
#endif
  }

  /* release resources */
  grid1d_wf_free(gwf);
  rgrid1d_free(workspace);
  cgrid1d_free(cworkspace);
  cgrid1d_free(sq_grad_pot);
  cgrid1d_free(potential);
  
  return 0;
}

double complex wavepacket(void *arg, double x) {

  double kx = ((wparams *) arg)->kx;
  double wx = ((wparams *) arg)->wx;
  double xc = ((wparams *) arg)->xc;
  
  x -= xc;

  return cexp(-x/wx*x/wx + I * kx * x);
}

double complex harmonic(void *arg, double x) {

  pparams params = *((pparams *) arg);

  return 0.5 * params.delta * params.kx * params.kx * x * x;
}

double complex dip(void *arg, double x) {

  double pot;
  pparams params = *((pparams *) arg);
  x = 2.0 * M_PI * x / (2.0 * params.kx); 
  pot = params.delta * (1.0 - cos(x) * cos(x));
  return pot;
}
