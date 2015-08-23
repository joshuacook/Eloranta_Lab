/*
 * Test program for itp_linear routines.
 *
 * Try for example:
 * ./test_wf1d_itp_linear 1 32 5 1 1.0 1e-4 1000 > test.dat
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grid/grid.h>
#include <omp.h>

double complex random_wf(void *arg, double x);

typedef struct pparams_struct {
  double kx;
  double delta;
} pparams;

double complex harmonic(void *arg, double x);
double complex dip(void *arg, double x);

double complex external_potential(void *arg, double x) {

  return dip(arg, x);
}

int main(int argc, char *argv[]) {

  long i, j, k, l, n, iterations, states, virtuals, threads;
  int riterations;
  double step, lx, threshold, rtau;
  double tau, erms;
  double complex time;
  wf1d **gwf = 0;
  cgrid1d *potential = 0;
  cgrid1d *workspace = 0;
  pparams potential_params;
  
  /* parameters */
  if (argc < 8) {
    fprintf(stderr, "Usage: %s <threads> <points/axis> <states> <virtuals> <time_step> <threshold> <iterations>\n", argv[0]);
    return 0;
  }
  
  threads = atoi(argv[1]);
  n = atoi(argv[2]);
  states = atoi(argv[3]);
  virtuals = atoi(argv[4]);
  tau = atof(argv[5]);
  threshold = atof(argv[6]);
  iterations = atol(argv[7]);
  
  step = 0.4 / (n / 16.0);
  lx = n * step;
  
  fprintf(stderr, "Grid (%ld)\n", n);
  
  /* potential parameters */
  potential_params.kx = lx;
  potential_params.delta = 10;
    
  grid_threads_init(threads);
  
  /* allocate memory */
  gwf = (wf1d **) malloc(states * sizeof(wf1d *));
  for(i = 0; i < states; i++)
    gwf[i] = grid1d_wf_alloc(n, step, 1.0, WF1D_PERIODIC_BOUNDARY, WF1D_2ND_ORDER_PROPAGATOR);
  potential = cgrid1d_alloc(n, step, CGRID1D_PERIODIC_BOUNDARY, 0);
  
  /* initialize wave function */
  for(i = 0; i < states; i++) {
    grid1d_wf_map(gwf[i], random_wf, 0);
    grid1d_wf_normalize(gwf[i]);
  }
  
  /* map potential */
  cgrid1d_smooth_map(potential, external_potential, &potential_params, 1);
  
  /* solve */
  erms = grid1d_itp_linear(gwf, states, virtuals, potential, tau, threshold, iterations, &rtau, &riterations);
  fprintf(stderr, "RMS of error = %le\n", erms); 
  
  /* print energies */
  workspace = cgrid1d_alloc(n, step, CGRID1D_PERIODIC_BOUNDARY, 0);
  
#if 1
  for(i = 0; i < states; i++)
    fprintf(stderr, " %24.16lf ", grid1d_wf_energy(gwf[i], potential, workspace));
  fprintf(stderr, "\n");
#if 0
  /* print wave function */
  for(i = 0; i < states; i++)
    grid1d_wf_print(gwf[i], stdout);
#endif
#endif
    
  /* release resources */
  for(i = 0; i < states; i++)
    grid1d_wf_free(gwf[i]);
  cgrid1d_free(workspace);
  cgrid1d_free(potential);
  
  free(gwf);
  
  return 0;
}


double complex random_wf(void *arg, double x) {

  return drand48();
}

double complex harmonic(void *arg, double x) {

  pparams params = *((pparams *) arg);

  return 0.5 * params.delta * (params.kx*params.kx * x*x);
}

double complex dip(void *arg, double x) {

  double pot;
  pparams params = *((pparams *) arg);

  x = 2.0 * M_PI * x / (2.0 * params.kx); 
  pot = params.delta * (1.0 - cos(x) * cos(x));
  return pot;
}
