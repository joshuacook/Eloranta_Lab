/*
 * Test program for itp_nonlinear routines.
 *
 * Try for example:
 * ./test_wf3d_itp_nonlinear 1 32 1 10000 1e-4 1000 > test.dat
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <grid/grid.h>
#include <grid/au.h>

double complex random_wf(void *arg, double x, double y, double z);

typedef struct lpparams_struct {
  double kx, ky, kz;
  double delta;
} lpparams;

typedef struct gpparams_struct {
  double mu, lambda, particles;
} gpparams;

typedef struct pparams_struct {
  lpparams linear;
  gpparams gp;
} pparams;

double complex harmonic(void *arg, double x, double y, double z);
double complex dip(void *arg, double x, double y, double z);
double complex external_potential(void *arg, double x, double y, double z) {

  return dip(arg, x, y, z);
}

void calculate_potentials(cgrid3d **potential, void *arg, wf3d **gwf, int states);

int main(int argc, char *argv[]) {

  long i, j, k, l, n, iterations, states = 1, virtuals = 0, threads, particles;
  int riterations;
  double step, lx, threshold;
  double tau, erms, rtau;
  double complex time;
  wf3d **gwf = 0;
  cgrid3d *workspace = 0, **potential = 0;
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
  
  fprintf(stderr, "Grid (%ldx%ldx%ld)\n", n, n, n);
  
  /* potential parameters */
  potential_params.linear.kx = lx;
  potential_params.linear.ky = lx;
  potential_params.linear.kz = lx;
  potential_params.linear.delta = 100.0 / GRID_AUTOK;
  potential_params.gp.mu = 7.0 / GRID_AUTOK;
  potential_params.gp.lambda = potential_params.gp.mu / (GRID_AUTOFS * (GRID_AUTOANG*GRID_AUTOANG*GRID_AUTOANG)); /* mu / rho_0 */
  potential_params.gp.mu = 0.0 * 7.0 / GRID_AUTOK;
  potential_params.gp.particles = particles;
  
  grid_threads_init(threads);
  
  /* allocate memory */
  gwf = (wf3d **) malloc(states * sizeof(wf3d *));
  potential = (cgrid3d **) malloc(states * sizeof(cgrid3d *));
  for(i = 0; i < states; i++) {
    gwf[i] = grid3d_wf_alloc(n, n, n, step, 4.0026 / GRID_AUTOAMU /*He*/, WF3D_PERIODIC_BOUNDARY, WF3D_2ND_ORDER_PROPAGATOR);
    potential[i] = cgrid3d_alloc(n, n, n, step, CGRID3D_PERIODIC_BOUNDARY, 0);
  }
  workspace = cgrid3d_alloc(n, n, n, step, CGRID3D_PERIODIC_BOUNDARY, 0);
  
  /* initialize wave function */
  for(i = 0; i < states; i++) {
    grid3d_wf_map(gwf[i], random_wf, 0);
    grid3d_wf_normalize(gwf[i]);
  }
  
  /* solve */
  erms = grid3d_itp_nonlinear(gwf, states, virtuals, calculate_potentials, &potential_params, tau, threshold, iterations, &rtau, &riterations);
  fprintf(stderr, "RMS of error = %le\n", erms);   

  calculate_potentials(potential, &potential_params, gwf, states);
  for (i = 0; i < states; i++)
    fprintf(stderr, " %24.16lf ", grid3d_wf_energy(gwf[i], potential[i], workspace));
  fprintf(stderr, "\n");

#if 0
  /* print wfs */  
  for(i = 0; i < states; i++)
    grid3d_wf_print(gwf[i], stdout);  
#endif
  
  /* release resources */
  for(i = 0; i < states; i++)
    grid3d_wf_free(gwf[i]);
  cgrid3d_free(workspace);
  
  free(gwf);
  
  return 0;
}

double complex random_wf(void *arg, double x, double y, double z) {

  return drand48();
}

double complex harmonic(void *arg, double x, double y, double z) {

  lpparams params = *((lpparams *) arg);

  return 0.5 * params.delta * (params.kx*params.kx * x*x + 
			       params.ky*params.ky * y*y + 
			       params.kz*params.kz * z*z);
}

double complex dip(void *arg, double x, double y, double z) {

  double pot;
  lpparams params = *((lpparams *) arg);

  x = 2.0 * M_PI * x / (2.0 * params.kx); 
  y = 2.0 * M_PI * y / (2.0 * params.ky); 
  z = 2.0 * M_PI * z / (2.0 * params.kz); 
  pot = params.delta * (1.0 - cos(x)*cos(x) * cos(y)*cos(y) * cos(z)*cos(z));
  return pot;
}

void calculate_potentials(cgrid3d **potential, void *arg, wf3d **gwf, int states) {

  int i;
  pparams *params = (pparams *) arg; 
  cgrid3d *tmp = cgrid3d_alloc(potential[0]->nx, potential[0]->ny, potential[0]->nz, potential[0]->step, CGRID3D_PERIODIC_BOUNDARY, 0);
  
  for(i = 0; i < states; i++) {
    cgrid3d_conjugate_product(potential[i], gwf[i]->grid, gwf[i]->grid);
    cgrid3d_multiply_and_add(potential[i], 
			    params->gp.particles * params->gp.lambda, 
			    - params->gp.mu);
    cgrid3d_map(tmp, external_potential, &params->linear);
    cgrid3d_sum(potential[i], potential[i], tmp);
  }
  
  cgrid3d_free(tmp);
}
