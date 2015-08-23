/*
 * Routines for handling 1D wavefunctions.
 *
 */

#include "grid.h"
#include "private.h"
#include "private1d.h"

/*
 * Allocate a 1D wavefunction.
 *
 * nx         = number of spatial grid points along x (int).
 * step       = spatial step size (double).
 * mass       = mass of the particle corresponding to this wavefunction (double).
 * boundary   = boundary condition (int):
 *              WF1D_DIRICHLET_BOUNDARY = Dirichlet boundary condition.
 *              WF1D_NEUMANN_BOUNDARY   = Neumann boundary condition.
 *              WF1D_PERIODIC_BOUNDARY  = Periodic boundary condition.
 * propagator = which time propagator to use for this wavefunction:
 *              WF1D_2ND_ORDER_PROPAGATOR = 2nd order in time.
 *              WF1D_4TH_ORDER_PROPAGATOR = 4th order in time.
 *
 * Return value is a pointer to the allocated wavefunction.
 * This routine returns NULL if allocation fails.
 *
 */

EXPORT wf1d *grid1d_wf_alloc(long nx, double step, double mass, int boundary, int propagator) {

  wf1d *wf;
  double complex (*value_outside)(const struct cgrid1d_struct *grid, long i);
  
  if (boundary != WF1D_DIRICHLET_BOUNDARY
       && boundary != WF1D_NEUMANN_BOUNDARY 
       && boundary != WF1D_PERIODIC_BOUNDARY) {
    fprintf(stderr, "libgrid: Error in grid1d_wf_alloc(). Unknown boundary condition (index = %d).\n", boundary);
    return 0;
  }
  
  if (propagator != WF1D_2ND_ORDER_PROPAGATOR
       && propagator != WF1D_4TH_ORDER_PROPAGATOR) {
    fprintf(stderr, "libgrid: Error in grid1d_wf_alloc(). Unknown propagator (index = %d).\n", propagator);
    return 0;
  }
  
  if ((boundary == WF1D_DIRICHLET_BOUNDARY || boundary == WF1D_NEUMANN_BOUNDARY)
       && propagator == WF1D_4TH_ORDER_PROPAGATOR) {
    fprintf(stderr, "libgrid: Error in grid1d_wf_alloc(). Invalid boundary condition - propagator combination. 4th order propagator can be used only with periodic boundary conditions.\n");
    return 0;
  }
  
  wf = (wf1d *) malloc(sizeof(wf1d));
  if (!wf) {
    fprintf(stderr, "libgrid: Error in grid1d_wf_alloc(). Could not allocate memory for wf1d.\n");
    return 0;
  }
  
  value_outside = NULL;
  if (boundary == WF1D_DIRICHLET_BOUNDARY)
    value_outside = cgrid1d_value_outside_constantdirichlet;
  else if (boundary == WF1D_NEUMANN_BOUNDARY)
    value_outside = cgrid1d_value_outside_neumann;
  else if (boundary == WF1D_PERIODIC_BOUNDARY)
    value_outside = cgrid1d_value_outside_periodic;
  
  wf->grid = cgrid1d_alloc(nx, step, value_outside, 0);
  
  if (!wf->grid) {
    fprintf(stderr, "libgrid: Error in grid1d_wf_alloc(). Could not allocate memory for wf1d->grid.\n");
    free(wf);
    return 0;
  }
  
  wf->mass = mass;
  wf->norm = 1.0;
  wf->boundary = boundary;
  wf->propagator = propagator;
  
  return wf;
}

/*
 * Free 1D wavefunction.
 *
 * gwf = wavefunction to be freed (wf1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_free(wf1d *gwf) {

  if (gwf) {
    if (gwf->grid) cgrid1d_free(gwf->grid);
    free(gwf);
  }
}

/* 
 * Calculate absorbing boundary potential.
 *
 * Excitations entering the "region" (see below) will be damped out
 * such that no back reflections occur from the boundary of the finite grid.
 * This is achieved by adding an imarinary potential which has the magnitude
 * depending on the difference between the present density (|psi_cur|^2) and 
 * the initial (desired) density (|psi_ini|^2). Where the difference between
 * the two is greatest, the largest damping occurs there.
 *
 * potential   = absorbing potential (cgrid1d *; output; overwritten).
 * density     = current density (rgrid1d *; input).
 * rho0        = desired density in the boundary region (double, input).
 * region      = function that will multiply the absorbing potential (double (*)(double); input).
 *               0.0 = no absorption, \approx 0.1 = full absorption (values this
 *               large often makes the time propagation unstable).
 *               Even a linear function starting from the absorbing region edge (value 0) 
 *               and then gradually increasing up to 0.1 will often work well.
 * workspace   = temporary space needed for the operation (rgrid1d *).
 *
 * The calculated potential should be added to the potential that will be
 * propagated. Note that it is important that this is also included in the
 * possible predict-correct cycle as that improves the numericla stability.
 *
 */

EXPORT void grid1d_wf_absorb(cgrid1d *potential, rgrid1d *density, double rho0, double (*region)(void *, double), rgrid1d *workspace) {

  rgrid1d_copy(workspace, density);
  rgrid1d_add(workspace, -rho0);
  rgrid1d_product_func(workspace, region, NULL);
  rgrid1d_multiply(workspace, -1.0);
  grid1d_add_real_to_complex_im(potential, workspace);
}

/*
 * Calculate the probability flux.
 *
 * gwf        = wavefunction for the operation (wf1d *).
 * flux_x     = output grid containing the flux (rgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_probability_flux(const wf1d *gwf, rgrid1d *flux_x) {

  /*
   * J(r) = -i (hbar/2m) (psi^* grad psi - psi grad psi^*)
   *      = (hbar/m) Im[psi^* grad psi] 
   */
  
#if 0
  // old code
  cgrid1d_fd_gradient_x(gwf->grid, workspace);
  cgrid1d_conjugate_product(workspace, gwf->grid, workspace);
  grid1d_complex_im_to_real(flux_x, workspace);
  rgrid1d_multiply(flux_x, HBAR / gwf->mass);
#endif
  // new code without additional workspace
  long i;
  cgrid1d *grid = gwf->grid;
  long nx = grid->nx;
  double inv_delta = 1.0 / (2.0 * grid->step), mass = gwf->mass;
  double *lvalue = flux_x->value;
  double complex (*value_at)(const cgrid1d *grid, long i) = cgrid1d_value_at_index;

#pragma omp parallel for firstprivate(value_at,nx,lvalue,inv_delta,grid,mass) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    lvalue[i] = cimag(conj(value_at(grid, i)) * inv_delta * (value_at(grid, i+1) - value_at(grid, i-1))) * (HBAR / mass);
}

/*
 * Calculate the momentum.
 *
 * gwf        = wavefunction for the operation (wf1d *).
 * momentum_x = output grid containing the momentum (cgrid1d *).
 * workspace  = additional storage needed for the operation (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_momentum(const wf1d *gwf, cgrid1d *momentum_x, cgrid1d *workspace) {

  cgrid1d_copy(workspace, gwf->grid);
  cgrid1d_fft(workspace);
  cgrid1d_fft_gradient(workspace, momentum_x);
  cgrid1d_inverse_fft(momentum_x);
  cgrid1d_multiply(momentum_x, -I * HBAR / (2.0 * gwf->mass));
}

/*
 * Calculate energy for the wavefunction.
 *
 * gwf       = wavefunction for the energy calculation (wf1d *).
 * potential = grid containing the potential (cgrid1d *).
 * workspace = additional storage needed (cgrid1d *).
 *
 * Returns the energy (double).
 *
 */

EXPORT double grid1d_wf_energy(const wf1d *gwf, const cgrid1d *potential, cgrid1d *workspace) {

  if (gwf->boundary == WF1D_DIRICHLET_BOUNDARY || gwf->boundary == WF1D_NEUMANN_BOUNDARY)
    return grid1d_wf_energy_cn(gwf, gwf, potential, workspace);
  else if (gwf->boundary == WF1D_PERIODIC_BOUNDARY)
    return grid1d_wf_energy_fft(gwf, potential, workspace);
  else
    abort();
}

/*
 * Auxiliary routine for calculating the energy (Crank-Nicolson).
 * Users should rather call grid1d_wf_energy().
 *
 * gwfa      = (left) wavefunction for the energy calculation (wf1d *).
 * gwfb      = (right) wavefunction for the energy calculation (wf1d *).
 *             Normally gwfa = gwfb.
 * potential = Potential grid (cgrid3d *).
 * workspace = Additional workspace needed for the operation (cgrid1d *).
 * 
 * Returns the energy.
 *
 */

EXPORT double grid1d_wf_energy_cn(const wf1d *gwfa, const wf1d *gwfb, const cgrid1d *potential, cgrid1d *workspace) {
  
  long i;

  /* (-2m/hbar^2) T psi */
  cgrid1d_fd_laplace(gwfb->grid, workspace);
  cgrid1d_multiply(workspace, -HBAR * HBAR / (2.0 * gwfb->mass));

  /* V psi */
  if(potential)
    for(i = 0; i < gwfb->grid->nx; i++)
      workspace->value[i] += potential->value[i] * gwfb->grid->value[i];
  
  /* int psi^* (T + V) psi d^3r */
  return creal(cgrid1d_integral_of_conjugate_product(gwfa->grid, workspace));
}

/*
 * Auxiliary routine for calculating the energy (FFT).
 * Users should rather call grid1d_wf_energy().
 *
 * gwf       = wavefunction for the energy calculation (wf1d *).
 * potential = Potential grid.
 * workspace = Additional workspace needed for the operation (rgrid1d *).
 * 
 * Returns the energy.
 *
 */

EXPORT double grid1d_wf_energy_fft(const wf1d *gwf, const cgrid1d *potential, cgrid1d *workspace) {

  double en;

  /* delta (- k^2) fft[f(x)] / N */
  cgrid1d_copy(workspace, gwf->grid);
  cgrid1d_fft(workspace);
  
  en = -HBAR*HBAR / (2.0 * gwf->mass) * cgrid1d_fft_laplace_expectation_value(workspace, workspace);
  if(potential) en += creal(cgrid1d_grid_expectation_value(gwf->grid, potential));
  return en;
}

/*
 * Auxiliary routine for calculating kinetic energy (FFT).
 * This is used by grid1d_wf_energy_fft().
 * 
 * gwf       = wavefunction for the kinetic energy calculation (wf1d *).
 * workspace = additional workspace required for the operation (cgrid1d *).
 *
 * Returns the kinetic energy.
 *
 */

EXPORT double grid1d_wf_kinetic_energy_fft(const wf1d *gwf, cgrid1d *workspace) {

  /* delta (- k^2) FFT[f(x)] / N */
  cgrid1d_copy(workspace, gwf->grid);
  cgrid1d_fft(workspace);
  
  return -HBAR * HBAR / (2.0 * gwf->mass) * cgrid1d_fft_laplace_expectation_value(workspace, workspace);
}

/*
 * Auxiliary routine for calculating potential energy.
 * 
 * gwf       = wavefunction for the kinetic energy calculation (wf1d *).
 * workspace = additional workspace required for the operation (cgrid1d *).
 *
 * Returns the potential energy.
 *
 */

EXPORT double grid1d_wf_potential_energy(const wf1d *gwf, const cgrid1d *potential) {

  return creal(cgrid1d_grid_expectation_value(gwf->grid, potential));
}

/*
 * Calculate energy (E) and the error (dE).
 *
 * gwf       = wavefunction for the operation (wf1d *).
 * potential = grid containing the potential (cgrid1d *).
 * workspace = additional workspace required for the operation (cgrid1d *).
 * error     = error estimate for energy (double *).
 * 
 * Returns the energy.
 *
 */

EXPORT double grid1d_wf_energy_and_error(const wf1d *gwf, const cgrid1d *potential, cgrid1d *workspace, double *error) { 
  
  double energy;

  /* 
   * energy and its error
   * dE^2 = int dE^2 |psi|^2 dtau = ... = int E_local^2 |psi|^2 dtau - E_avg^2
   * dE = E_local - E_avg
   * E_local psi = H psi
   *
   */

  /* T psi */
  if (gwf->boundary == WF1D_DIRICHLET_BOUNDARY || gwf->boundary == WF1D_NEUMANN_BOUNDARY) {
    cgrid1d_fd_laplace(gwf->grid, workspace);
    cgrid1d_multiply(workspace, -HBAR * HBAR / (2.0 * gwf->mass));
  }
  else if (gwf->boundary == WF1D_PERIODIC_BOUNDARY) {
    cgrid1d_copy(workspace, gwf->grid);
    cgrid1d_fft(workspace);
    cgrid1d_fft_laplace(workspace, workspace);
    cgrid1d_scaled_inverse_fft(workspace, -HBAR * HBAR / (2.0 * gwf->mass));
  }
  else {
    fprintf(stderr, "libgrid: Error in grid1d_wf_energy_and_error(). Invalid boundary condition, index = %d\n", gwf->boundary);
    abort();
  }
  
  /* H psi */
  cgrid1d_add_scaled_product(workspace, 1.0, potential, gwf->grid);
  
  /* int E_local^2 |psi|^2 dtau */
  *error = cgrid1d_integral_of_square(workspace);
  
  /* int E_local |psi|^2 dtau */
  cgrid1d_conjugate_product(workspace, gwf->grid, workspace);
  energy = creal(cgrid1d_integral(workspace));
  
  /* sqrt(int E_local^2 |psi|^2 dtau - (int E_local |psi|^2 dtau)^2) */
  *error = sqrt(*error - energy * energy );
  
  return energy;
}

/*
 * Propagate wavefunction in time subject to given potential.
 *
 * gwf         = wavefunction to be propagated (wf1d *).
 * potential   = grid containing the potential (cgrid1d *).
 * sq_grad_pot = grid containing square of potential gradient (cgrid1d *).
 * time        = time step (double complex). Note this may be either real or imaginary.
 * workspace   = additional workspace needed for the operation (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_propagate(wf1d *gwf, const cgrid1d *potential, const cgrid1d *sq_grad_pot, double complex time, cgrid1d *workspace) {  
  
  double complex half_time = 0.5 * time;
  double complex one_sixth_time = time / 6.0;
  double complex two_thirds_time = 2.0 * time / 3.0;
  
  if ((gwf->boundary == WF1D_DIRICHLET_BOUNDARY || gwf->boundary == WF1D_NEUMANN_BOUNDARY)
       && gwf->propagator == WF1D_2ND_ORDER_PROPAGATOR) {
    grid1d_wf_propagate_potential(gwf, potential, half_time);
    grid1d_wf_propagate_kinetic_cn(gwf, time, workspace);
    grid1d_wf_propagate_potential(gwf, potential, half_time);
  } else if (gwf->boundary == WF1D_PERIODIC_BOUNDARY 
            && gwf->propagator == WF1D_2ND_ORDER_PROPAGATOR) {
    grid1d_wf_propagate_potential(gwf, potential, half_time);
    grid1d_wf_propagate_kinetic_fft(gwf, time);
    grid1d_wf_propagate_potential(gwf, potential, half_time);
  } else if (gwf->boundary == WF1D_PERIODIC_BOUNDARY 
            && gwf->propagator == WF1D_4TH_ORDER_PROPAGATOR) {
    grid1d_wf_propagate_potential(gwf, potential, one_sixth_time);
    grid1d_wf_propagate_kinetic_fft(gwf, half_time);
    cgrid1d_copy(workspace, potential);
    cgrid1d_add_scaled(workspace, (1.0 / 48.0 * HBAR * HBAR / gwf->mass) * sqnorm(time), sq_grad_pot);	
    grid1d_wf_propagate_potential(gwf, workspace, two_thirds_time);    
    grid1d_wf_propagate_kinetic_fft(gwf, half_time);
    grid1d_wf_propagate_potential(gwf, potential, one_sixth_time);
  } else {
    fprintf(stderr, "libgrid: Error in grid1d_wf_propagate(). Unknown propagator - boundary value combination (propagator index = %d, boundary index = %d).\n", gwf->propagator, gwf->boundary);
    abort();
  }
}

/*
 * Auxiliary routine to propagate kinetic energy using FFT.
 * Users should call grid1d_wf_propaget().
 *
 * gwf  = wavefunction to be propagated (wf1d *).
 * time = time step (double complex).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_propagate_kinetic_fft(wf1d *gwf, double complex time) {

  long i, nx;
  double kx, step, norm, mass;
  double complex *value = gwf->grid->value;
  
  nx = gwf->grid->nx;
  step = gwf->grid->step;
  mass = gwf->mass;

  cgrid1d_fft(gwf->grid);
  
  /* f(x) = iFFT[FFT[f(x)]] / N */
  norm = 1.0 / nx;
  
#pragma omp parallel for firstprivate(norm,nx,step,value,time,mass) private(i,kx) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {

    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (i < nx / 2)
      kx = 2.0 * M_PI * i / (nx * step);
    else 
      kx = 2.0 * M_PI * (i - nx) / (nx * step);

    /* psi(t+dt) = psi(t) exp( - i (hbar^2 * k^2 / 2m) dt / hbar ) */	  
    value[i] *= norm * cexp(-I * time * HBAR * HBAR * kx * kx / (2.0 * mass * HBAR));
  }
  
  cgrid1d_inverse_fft(gwf->grid);
}

/*
 * Auxiliary routine to propagate kinetic energy (Crank-Nicolson).
 * Users should call grid1d_wf_propaget().
 *
 * gwf       = wavefunction to be propagated (wf1d *).
 * time      = time step (double complex).
 * workspace = additional workspace required for the operation (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_propagate_kinetic_cn(wf1d *gwf, double complex time, cgrid1d *workspace) {

  if (gwf->boundary == WF1D_DIRICHLET_BOUNDARY)
    grid1d_wf_propagate_kinetic_cn_dbc(gwf, time, workspace);
  else if (gwf->boundary == WF1D_NEUMANN_BOUNDARY)
    grid1d_wf_propagate_kinetic_cn_nbc(gwf, time, workspace);
  else if (gwf->boundary == WF1D_PERIODIC_BOUNDARY)
    grid1d_wf_propagate_kinetic_cn_pbc(gwf, time, workspace);
  else {
    fprintf(stderr, "libgrid: Error in grid1d_wf_propagate_cn(). Unknown boundary condition (index = %d).\n", gwf->boundary);
    abort();
  }
}

/*
 * Auxiliary routine to propagate potential energy.
 * Users should call grid1d_wf_propaget().
 *
 * gwf       = wavefunction to be propagated (wf1d *).
 * potential = grid containing the potential (cgrid1d *).
 * time      = time step (double complex).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_propagate_potential(wf1d *gwf, const cgrid1d *potential, double complex time) {

  long i, nx = gwf->grid->nx;
  double complex c, *psi = gwf->grid->value, *pot = potential->value;
  
  c = -I * time / HBAR;
  
#pragma omp parallel for firstprivate(nx,psi,pot,c) private(i) default(none) schedule(runtime)
  for (i = 0; i < nx; i++)
    /* psi(t+dt) = exp( - i V dt / hbar ) psi(t) */
    psi[i] *= cexp(c*pot[i]);
}

/*
 * Project "gwfb" out from "gwfa".
 *
 * gwfa = input wavefunction (wf1d *).
 * gwfb = this will be projected out from gwfa (wf1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_project_out(wf1d *gwfa, const wf1d *gwfb) {

  double complex overlap = grid1d_wf_overlap(gwfa, gwfb);

  cgrid1d_add_scaled(gwfa->grid, -overlap, gwfb->grid);
}

/*
 * "traditional" diagonalization of Hamiltonian.
 *
 * gwf    = an array of wavefunctions (wf1d **).
 * states = number of states (int).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_diagonalize(wf1d **gwf, int states) {

  int i, j;
  double *eigenvalue = (double *) malloc(states * sizeof(double));
  double complex *overlap = (double complex *) malloc(states * states * sizeof(double complex));
  wf1d *gwf_tmp;
  
  if (states == 1) {
    grid1d_wf_normalize(gwf[0]);
    return;
  }
  
  /* overlap matrix */
  for( i = 0; i < states; i++ ) {
    for( j = 0; j <= i; j++ ) {	  
      /* fortran (column major) matrix order, i is row (minor) index, j is column (major) index */
      overlap[i + j * states] = grid1d_wf_overlap(gwf[i], gwf[j]);
      overlap[j + i * states] = conj(overlap[i + j * states]);
    }
  }
  
  /* diagonalize */
  grid_hermitian_eigenvalue_problem(eigenvalue, overlap, states);
  
  /* phi_i = 1 / sqrt(m_i) C_ij psi_j, C (row major) matrix order ???is it??? */
  for(i = 0; i < states; i++)
    for(j = 0; j < states; j++)
      overlap[i * states + j] /= sqrt(eigenvalue[i]);
  
  grid1d_wf_linear_transform(gwf, overlap, states);
  
  /* invert order */
  for(i = 0; i < states/2; i++) {
    gwf_tmp = gwf[i];
    gwf[i] = gwf[states-i-1];
    gwf[states-i-1] = gwf_tmp;	
  }
  
  /* free memory */
  free(eigenvalue);
  free(overlap);
}

/*
 * Linear transform a set of wavefunctions.
 *
 * gwf       = an array of wavefunctions (wf1d **).
 * transform = transformation matrix (double complex *).
 * states    = number of states (int).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_linear_transform(wf1d **gwf, double complex *transform, int states) {

  int p, q, offset;
  long i, nx;
  double complex **value, *tmp;
  
  nx = gwf[0]->grid->nx;
  
  /* + 16 to prevent "write locks" */
  tmp = (double complex *) malloc(omp_get_max_threads() * (states + 16) * sizeof(double complex));
  
  value = (double complex **) malloc(states * sizeof(double complex *));
  for(p = 0; p < states; p++)
    value[p] = gwf[p]->grid->value;
  
#pragma omp parallel for firstprivate(nx,value,states,transform,tmp) private(i,p,q,offset) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {
    offset = (states + 16) * omp_get_thread_num();
    for(p = 0; p < states; p++)
      tmp[offset + p] = 0.0;
    
    for(p = 0; p < states; p++)
      for(q = 0; q < states; q++)
	tmp[offset + p] += transform[p * states + q] * value[q][i];
    
    for(p = 0; p < states; p++)
      value[p][i] = tmp[offset + p];    
  }

  free(value);
  free(tmp);
}

/*
 * Calcuate square of potential gradient.
 *
 * sq_grad_pot = output grid (cgrid1d *).
 * potential   = potental input grid (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_square_of_potential_gradient(cgrid1d *sq_grad_pot, const cgrid1d *potential) {

  cgrid1d_copy(sq_grad_pot, potential);
  cgrid1d_fft(sq_grad_pot);
  cgrid1d_fft_gradient(sq_grad_pot, sq_grad_pot);
  
  cgrid1d_inverse_fft(sq_grad_pot);
  
  cgrid1d_conjugate_product(sq_grad_pot, sq_grad_pot, sq_grad_pot);
}

/*
 * Auxiliary routine for propagating kinetic energy using Dirichlet boundaries (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf1d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_propagate_kinetic_cn_dbc(wf1d *gwf, double complex time, cgrid1d *workspace) {

  /*
   * exp( -i Tx dt / hbar ) 
   *   
   */
  
  grid1d_wf_propagate_kinetic_x_cn_dbc(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along x using Dirichlet boundaries (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf1d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_propagate_kinetic_x_cn_dbc(wf1d *gwf, double complex time, cgrid1d *workspace) {

  double complex c, ca, boundary, *ad, *ld, *b, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long i, nx = gwf->grid->nx;

  /*
   * (1 + .5 i T dt) psi(t+dt) = (1 - .5 i T dt) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar^2 dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * HBAR * time);
  
  /* create Cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;
  /* in 1d we ran out of workspace... fix later */
  if(!(ad = (double complex *) malloc(sizeof(double complex) * nx))) {
    fprintf(stderr, "libgrid: Out of memory in grid1d_wf_propagate_kinetic_x_cn_dbc().\n");
    exit(1);
  }
  
  /* create matrix */
#pragma omp parallel for firstprivate(ad, ca, nx) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    ad[i] = ca;

  /* create decomposition */
  grid_cholesky_decomposition(ad, nx, 0);

  /* matrix is replaced by decomposition */
  ld = ad;
  
  b = wrk;
  
  /* create right-hand side vector */
#pragma omp parallel for firstprivate(c, b, nx, psi) private(i) default(none) schedule(runtime)  
  for(i = 1; i < nx-1; i++)
    /* (C - dx^2 laplace) psi(t) */
    b[i] = c * psi[i] + 2.0 * psi[i] - (psi[i+1] + psi[i-1]);
  
  /* dirichlet boundaries */
  boundary = 0.0;
  i = 0;
  b[i] = c * psi[i] + 2.0 * psi[i] - (psi[i+1] + boundary);
  i = nx-1;
  b[i] = c * psi[i] + 2.0 * psi[i] - (boundary + psi[i-1]);
  
  /* solve */
  grid_cholesky_substitute(ld, b, nx, 0);
  
  /* copy */
#pragma omp parallel for firstprivate(psi, b, nx) private(i) default(none) schedule(runtime) 
  for(i = 0; i < nx; i++)
    psi[i] = b[i];

  free(ad); /* goes with the malloc above */
}

/*
 * Auxiliary routine for propagating kinetic energy using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf1d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT void grid1d_wf_propagate_kinetic_cn_nbc(wf1d *gwf, double complex time, cgrid1d *workspace) {

  /*
   * exp( -i Tx dt / hbar ) 
   *   
   */
  
  grid1d_wf_propagate_kinetic_x_cn_dbc(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along x using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf1d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid1d *).
 *
 * No return value.
 *
 * TODO: could add openmp pragmas but 1D problems are usually small anyway...
 *
 */

EXPORT void grid1d_wf_propagate_kinetic_x_cn_nbc(wf1d *gwf, double complex time, cgrid1d *workspace) {

  double complex c, ca, *psi = gwf->grid->value;
  double step = gwf->grid->step;
  long i, nx = gwf->grid->nx;
  double complex *du = workspace->value;
  double complex *dl;
  double complex *d;
  double complex *b;

  /* in 1D we don't have enough space in workspace, we need to allocate more memory */
  if(!(dl = (double complex *) malloc(sizeof(double complex) * nx))) {
    fprintf(stderr, "libgrid: Out of memory in grid1d_wf_propagate_kinetic_x_cn_nbc().\n");
    exit(1);
  }
  if(!(d = (double complex *) malloc(sizeof(double complex) * nx))) {
    fprintf(stderr, "libgrid: Out of memory in grid1d_wf_propagate_kinetic_x_cn_nbc().\n");
    exit(1);
  }
  if(!(b = (double complex *) malloc(sizeof(double complex) * nx))) {
    fprintf(stderr, "libgrid: Out of memory in grid1d_wf_propagate_kinetic_x_cn_nbc().\n");
    exit(1);
  }
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* create cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

  fprintf(stderr, "libgrid: Using Thomas algorithm.\n");   /* LU not implemented.. would need even more workspace */

  for(i = 0; i < nx; i++) {
      
    if(i == nx-1) dl[nx-1] = 2.0;             /* Neumann outer boundary */
    else if(i > 0) dl[i] = 1.0;
    
    if(i == 0) du[0] = 2.0;                   /* Neumann inner boundary */
    else if(i < nx-1) du[i] = 1.0;
    
    d[i] = ca;
  }

  /* create right-hand side vector */
  for(i = 1; i < nx - 1; i++) {
    /* (C - dx^2 laplace) psi(t) */
    b[i] = c * psi[i] - (psi[i + 1] - 2.0 * psi[i] + psi[i - 1]);
  }
 
  /* neumann boundaries */
  i = 0; 
  b[i] = c * psi[i] - (2.0 * psi[i + 1] - 2.0 * psi[i]);
  
  i = nx - 1;
  b[i] = c * psi[i] - (2.0 * psi[i - 1] - 2.0 * psi[i]);
  
  grid_solve_tridiagonal_system(nx, dl, d, du, b, dl);
  
  for(i = 0; i < nx; i++)
    psi[i] = dl[i];
  free(dl);
  free(d);
  free(b);
}

/*
 * Auxiliary routine for propagating kinetic energy using periodic boundary condition (Crank-Nicolson).
 * 
 * gwf       = wavefunction to be propagated (wf1d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid1d *).
 *
 * No return value.
 * 
 */

EXPORT void grid1d_wf_propagate_kinetic_cn_pbc(wf1d *gwf, double complex time, cgrid1d *workspace) {

  /*
   * exp(-i Tx dt / hbar) 
   *
   */

  grid1d_wf_propagate_kinetic_x_cn_pbc(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along x using periodic boundary condition (Crank-Nicolson).
 * 
 * gwf       = wavefunction to be propagated (wf1d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (grid1d *).
 *
 * No return value.
 * 
 */

EXPORT void grid1d_wf_propagate_kinetic_x_cn_pbc(wf1d *gwf, double complex time, cgrid1d *workspace) {

  double complex c, ca, boundary, *ad, *ld, *b, *w, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long i, nx = gwf->grid->nx;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / hbar dt)
   */
  c = (I * 4.0 * gwf->mass * step * step / HBAR) / time;
  
  /* create Cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;
  /* TODO: running out of space in 1D ... */
  if(!(ad = (double complex *) malloc(sizeof(double complex) * nx))) {
    fprintf(stderr, "libgrid: Out of memory in grid1d_wf_propagate_kinetic_x_cn_pbc().\n");
    exit(1);
  }
  
  /* create matrix A' */
#pragma omp parallel for firstprivate(ad, ca, nx) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    ad[i] = ca;  
  /* periodic boundaries */
  ad[0]    -= 1.0;
  ad[nx-1] -= 1.0;
  
  /* create Cholesky decomposition */
  grid_cholesky_decomposition(ad, nx, 0);

  /* matrix is replaced by decomposition */
  ld = ad;
  
  /* periodic boundaries */
  /* TODO: running out of space in 1D ... */
  if(!(w = (double complex *) malloc(sizeof(double complex) * nx))) {
    fprintf(stderr, "libgrid: Out of memory in grid1d_wf_propagate_kinetic_x_cn_pbc().\n");
    exit(1);
  }
#pragma omp parallel for firstprivate(w, nx) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) 
    w[i] = 0.0;
  w[0] = w[nx-1] = 1.0;

  /* solve A' z = u */
  grid_cholesky_substitute(ld, w, nx, 0);
  
  /* use diffent workspace for different threads */
  b = wrk;
      
  /* create right-hand side vector */
#pragma omp parallel for firstprivate(c, b, nx, psi) private(i) default(none) schedule(runtime)  
  for(i = 1; i < nx-1; i++)
    /* (C - dx^2 laplace) psi(t) */
    b[i] = c * psi[i] + 2.0 * psi[i] - (psi[i+1] + psi[i-1]);
  
  /* periodic boundaries */
  boundary = psi[nx-1];
  i = 0;
  b[i] = c * psi[i] + 2.0 * psi[i] - (psi[i+1] + boundary);
  
  boundary = psi[0];
  i = nx-1;
  b[i] = c * psi[i] + 2.0 * psi[i] - (boundary + psi[i-1]);
  
  /* solve A' x' = b */
  grid_cholesky_substitute(ld, b, nx, 0);
      
  /* periodic boundaries */
  /* shermann-morrison form NR in C
   *             v . x'                x'_0 + x'_(n-1)
   * x = x' - ----------- w = x' - -------------------- w
   *           1 + v . w             1 + w_0 + w_(n-1)
   *
   * where v = (1, 0, 0, 0, ..., 0, 1) 
   */	  
  ca = (b[0] + b[nx-1]) / (1.0 + w[0] + w[nx-1]);
#pragma omp parallel for firstprivate(b, nx, ca, w) private(i) default(none) schedule(runtime) 
  for(i = 0; i < nx; i++)
    b[i] -= ca * w[i];
  
  /* copy back */
#pragma omp parallel for firstprivate(psi, b, nx) private(i) default(none) schedule(runtime) 
  for(i = 0; i < nx; i++)
    psi[i] = b[i];

  free(ad);
  free(w);
}

/*
 * Produce density grid from a given wavefunction.
 *
 * gwf     = wavefunction (wf1d *).
 * density = output density grid (cgrid1d *).
 *
 * No return value.
 *
 */

EXPORT inline void grid1d_wf_density(const wf1d *gwf, rgrid1d *density) {

  long i, nx = gwf->grid->nx;
  double complex *avalue = gwf->grid->value;
  double *cvalue = density->value;
  
#pragma omp parallel for firstprivate(nx,avalue,cvalue) private(i) default(none) schedule(runtime)
  for(i = 0; i < nx; i++)
    cvalue[i] = conj(avalue[i]) * avalue[i];
}

/*
 * Zero wavefunction.
 *
 * gwf = wavefunction to be zeroed (wf3d *).
 *
 * No return value.
 *
 */

EXPORT inline void grid1d_wf_zero(wf1d *gwf) { 

  cgrid1d_zero(gwf->grid); 
}

/*
 * Set wavefunction to some constant value.
 *
 * gwf = wavefunction to be set (wf1d *).
 * c   = value (double complex).
 *
 * No return value.
 *
 */

EXPORT inline void grid1d_wf_constant(wf1d *gwf, double complex c) { 

  cgrid1d_constant(gwf->grid, c); 
}

/*
 * Map a given function on a wavefunction.
 *
 * gwf  = wavefunction where function will be mapped to (wf1d *).
 * func = function providing the mapping (double complex (*)(void *, double)).
 * farg = optional argument for passing parameters to func (void *).
 *
 * No return value.
 *
 */

EXPORT inline void grid1d_wf_map(wf1d *gwf, double complex (*func)(void *arg, double x), void *farg) { 

  cgrid1d_map(gwf->grid, func, farg); 
}

/*
 * Calculate the norm of the given wavefunction.
 *
 * gwf = wavefunction for the calculation (wf1d *).
 *
 * Returns the norm (double).
 *
 */

EXPORT inline double grid1d_wf_norm(const wf1d *gwf) { 

  return cgrid1d_integral_of_square(gwf->grid); 
}

/*
 * Normalize wavefunction (to the value given in gwf->norm).
 *
 * gwf = wavefunction to be normalized (wf1d *).
 *
 * Returns the normalization constant (double).
 *
 */

EXPORT inline double grid1d_wf_normalize(wf1d *gwf) { 

  double norm = grid1d_wf_norm(gwf);

  cgrid1d_multiply(gwf->grid, sqrt(gwf->norm / norm));
  return norm; 
}

/*
 * Calculate overlap between two wavefunctions.
 *
 * gwfa = 1st wavefunction (wf1d *).
 * gwfb = 2nd wavefunction (wf1d *).
 *
 * Returns the overlap (double complex).
 *
 */

EXPORT inline double complex grid1d_wf_overlap(const wf1d *gwfa, const wf1d *gwfb) { 

  return cgrid1d_integral_of_conjugate_product(gwfa->grid, gwfb->grid); 
}

/*
 * Output wavefunction.
 *
 * gwf = wavefunction to be printed (wf1d *).
 * out = output file pointer (FILE *).
 *
 * No return value.
 *
 */

EXPORT inline void grid1d_wf_print(const wf1d *gwf, FILE *out) { 

  cgrid1d_print(gwf->grid, out); 
}
