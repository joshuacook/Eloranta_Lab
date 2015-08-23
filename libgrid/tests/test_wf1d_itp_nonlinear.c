/*
 * Test program for itp_nonlinear routines.
 *
 * Try for example:
 * ./test_wf1d_itp_nonlinear 1 32 1 10000 1e-4 1000 > test.dat
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grid/grid.h>
#include <grid/au.h>

double complex random_wf(void *arg, double x);

typedef struct lpparams_struct {
  double kx;
  double delta;
} lpparams;

typedef struct gpparams_struct {
  double mu, lambda, particles;
} gpparams;

typedef struct pparams_struct {
  lpparams linear;
  gpparams gp;
} pparams;

double complex harmonic(void *arg, double x);
double complex dip(void *arg, double x);
double complex external_potential(void *arg, double x) {

  return dip(arg, x);
}

void calculate_potentials(cgrid1d **potential, void *arg, wf1d **gwf, int states);

int main(int argc, char *argv[]) {

  long i, j, k, l, n, iterations, states = 1, virtuals = 0, threads, particles;
  int riterations;
  double step, lx, threshold;
  double tau, erms, rtau;
  double complex time;
  wf1d **gwf = 0;
  cgrid1d *workspace = 0, **potential = 0;
  pparams potential_params;
  
  /* parameters */
  if (argc < 6) {
    fprintf(stderr, "Usage: %s <threads> <points/axis> <particles> <time_step> <threshold> <iterations>\n", argv[0]);
    return 0;
  }
  
  threads = atoi(argv[1]);
  n = atoi(argv[2]);
  particles = atoi(argv[3]);
  tau = atof(argv[4]);
  threshold = atof(argv[5]);
  iterations = atol(argv[6]);
  
  step = 1.6 / (n / 16.0);
  lx = n * step;
  
  fprintf(stderr, "Grid (%ld)\n", n);
  
  /* potential parameters */
  potential_params.linear.kx = lx;
  potential_params.linear.delta = 100.0 / GRID_AUTOK;
  potential_params.gp.mu = 7.0 / GRID_AUTOK;
  potential_params.gp.lambda = potential_params.gp.mu / (GRID_AUTOFS * (GRID_AUTOANG*GRID_AUTOANG*GRID_AUTOANG)); /* mu / rho_0 */
  potential_params.gp.mu = 0.0 * 7.0 / GRID_AUTOK;
  potential_params.gp.particles = particles;
  
  grid_threads_init(threads);
  
  /* allocate memory */
  gwf = (wf1d **) malloc(states * sizeof(wf1d *));
  potential = (cgrid1d **) malloc(states * sizeof(cgrid1d *));
  for(i = 0; i < states; i++) {
    gwf[i] = grid1d_wf_alloc(n, step, 4.0026 / GRID_AUTOAMU /*He*/, WF1D_PERIODIC_BOUNDARY, WF1D_2ND_ORDER_PROPAGATOR);
    potential[i] = cgrid1d_alloc(n, step, CGRID1D_PERIODIC_BOUNDARY, 0);
  }
  workspace = cgrid1d_alloc(n, step, CGRID1D_PERIODIC_BOUNDARY, 0);
  
  /* initialize wave function */
  for(i = 0; i < states; i++) {
    grid1d_wf_map(gwf[i], random_wf, 0);
    grid1d_wf_normalize(gwf[i]);
  }
  
  /* solve */
  erms = grid1d_itp_nonlinear(gwf, states, virtuals, calculate_potentials, &potential_params, tau, threshold, iterations, &rtau, &riterations);
  fprintf(stderr, "RMS of error = %le\n", erms);   

  calculate_potentials(potential, &potential_params, gwf, states);
  for (i = 0; i < states; i++)
    fprintf(stderr, " %24.16lf ", grid1d_wf_energy(gwf[i], potential[i], workspace));
  fprintf(stderr, "\n");

#if 0
  /* print wfs */  
  for(i = 0; i < states; i++)
    grid1d_wf_print(gwf[i], stdout);  
#endif
  
  /* release resources */
  for(i = 0; i < states; i++)
    grid1d_wf_free(gwf[i]);
  cgrid1d_free(workspace);
  
  free(gwf);
  
  return 0;
}

double complex random_wf(void *arg, double x) {

  return drand48();
}

double complex harmonic(void *arg, double x) {

  lpparams params = *((lpparams *) arg);

  return 0.5 * params.delta * (params.kx*params.kx * x*x);
}

double complex dip(void *arg, double x) {

  double pot;
  lpparams params = *((lpparams *) arg);

  x = 2.0 * M_PI * x / (2.0 * params.kx); 
  pot = params.delta * (1.0 - cos(x)*cos(x));
  return pot;
}

void calculate_potentials(cgrid1d **potential, void *arg, wf1d **gwf, int states) {

  int i;
  pparams *params = (pparams *) arg; 
  cgrid1d *tmp = cgrid1d_alloc(potential[0]->nx, potential[0]->step, CGRID1D_PERIODIC_BOUNDARY, 0);
  
  for(i = 0; i < states; i++) {
    cgrid1d_conjugate_product(potential[i], gwf[i]->grid, gwf[i]->grid);
    cgrid1d_multiply_and_add(potential[i], 
			    params->gp.particles * params->gp.lambda, 
			    - params->gp.mu);
    cgrid1d_map(tmp, external_potential, &params->linear);
    cgrid1d_sum(potential[i], potential[i], tmp);
  }
  
  cgrid1d_free(tmp);
}
