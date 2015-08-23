/*
 * Test program for itp_linear routines.
 *
 * Try for example:
 * ./test_wf2d_itp_linear 1 32 5 1 1.0 1e-4 1000 > test.dat
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grid/grid.h>

double complex random_wf(void *arg, double x, double y);

typedef struct pparams_struct {
  double kx, ky;
  double delta;
} pparams;

double complex harmonic(void *arg, double x, double y);
double complex dip(void *arg, double x, double y);

double complex external_potential(void *arg, double x, double y) {

  return dip(arg,x,y);
}

int main(int argc, char *argv[]) {

  long i, j, k, l, n, iterations, states, virtuals, threads;
  int riterations;
  double step, lx, threshold, rtau;
  double tau, erms;
  double complex time;
  wf2d **gwf = 0;
  cgrid2d *potential = 0;
  cgrid2d *workspace = 0;
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
  
  fprintf(stderr, "Grid (%ldx%ld)\n", n, n);
  
  /* potential parameters */
  potential_params.kx = lx;
  potential_params.ky = lx;
  potential_params.delta = 10;
    
  grid_threads_init(threads);
  
  /* allocate memory */
  gwf = (wf2d **) malloc(states * sizeof(wf2d *));
  for(i = 0; i < states; i++)
    gwf[i] = grid2d_wf_alloc(n, n, step, 1.0, WF2D_PERIODIC_BOUNDARY, WF2D_2ND_ORDER_PROPAGATOR);
  potential = cgrid2d_alloc(n, n, step, CGRID2D_PERIODIC_BOUNDARY, 0);
  
  /* initialize wave function */
  for(i = 0; i < states; i++) {
    grid2d_wf_map(gwf[i], random_wf, 0);
    grid2d_wf_normalize(gwf[i]);
  }
  
  /* map potential */
  cgrid2d_smooth_map(potential, external_potential, &potential_params, 1);
  
  /* solve */
  erms = grid2d_itp_linear(gwf, states, virtuals, potential, tau, threshold, iterations, &rtau, &riterations);
  fprintf(stderr, "RMS of error = %le\n", erms); 
  
  /* print energies */
  workspace = cgrid2d_alloc(n, n, step, CGRID2D_PERIODIC_BOUNDARY, 0);
  
#if 1
  for(i = 0; i < states; i++)
    fprintf(stderr, " %24.16lf ", grid2d_wf_energy(gwf[i], potential, workspace));
  fprintf(stderr, "\n");
#if 0
  /* print wave function */
  for(i = 0; i < states; i++)
    grid2d_wf_print(gwf[i], stdout);
#endif
#endif
    
  /* release resources */
  for(i = 0; i < states; i++)
    grid2d_wf_free(gwf[i]);
  cgrid2d_free(workspace);
  cgrid2d_free(potential);
  
  free(gwf);
  
  return 0;
}


double complex random_wf(void *arg, double x, double y) {

  return drand48();
}

double complex harmonic(void *arg, double x, double y) {

  pparams params = *((pparams *) arg);

  return 0.5 * params.delta * (params.kx*params.kx * x*x + 
			       params.ky*params.ky * y*y);
}

double complex dip(void *arg, double x, double y) {

  double pot;
  pparams params = *((pparams *) arg);

  x = 2.0 * M_PI * x / (2.0 * params.kx); 
  y = 2.0 * M_PI * y / (2.0 * params.ky); 
  pot = params.delta * (1.0 - cos(x) * cos(x) * cos(y) * cos(y));
  return pot;
}
