/*
 * Routines for handling 2D wavefunctions.
 *
 */

#include "grid.h"
#include "private.h"
#include "private2d.h"

/*
 * Allocate a 2D wavefunction.
 *
 * nx         = number of spatial grid points along x (int).
 * ny         = number of spatial grid points along y (int).
 * step       = spatial step size (double).
 * mass       = mass of the particle corresponding to this wavefunction (double).
 * boundary   = boundary condition (int):
 *              WF2D_DIRICHLET_BOUNDARY = Dirichlet boundary condition.
 *              WF2D_NEUMANN_BOUNDARY   = Neumann boundary condition.
 *              WF2D_PERIODIC_BOUNDARY  = Periodic boundary condition.
 * propagator = which time propagator to use for this wavefunction:
 *              WF2D_2ND_ORDER_PROPAGATOR = 2nd order in time.
 *              WF2D_4TH_ORDER_PROPAGATOR = 4th order in time.
 *
 * Return value is a pointer to the allocated wavefunction.
 * This routine returns NULL if allocation fails.
 *
 */

EXPORT wf2d *grid2d_wf_alloc(long nx, long ny, double step, double mass, int boundary, int propagator) {

  wf2d *wf;
  double complex (*value_outside)(const struct cgrid2d_struct *grid, long i, long j);
  
  if(boundary != WF2D_DIRICHLET_BOUNDARY 
       && boundary != WF2D_NEUMANN_BOUNDARY 
       && boundary != WF2D_PERIODIC_BOUNDARY) {
    fprintf(stderr, "libgrid: Error in grid2d_wf_alloc(). Unknown boundary condition (index = %d).\n", boundary);
    return 0;
  }
  
  if(propagator != WF2D_2ND_ORDER_PROPAGATOR 
       && propagator != WF2D_4TH_ORDER_PROPAGATOR) {
    fprintf(stderr, "libgrid: Error in grid2d_wf_alloc(). Unknown propagator (index = %d).\n", propagator);
    return 0;
  }
  
  if ((boundary == WF2D_DIRICHLET_BOUNDARY || boundary == WF2D_NEUMANN_BOUNDARY)
      && propagator == WF2D_4TH_ORDER_PROPAGATOR) {
    fprintf(stderr, "libgrid: Error in grid2d_wf_alloc(). Invalid boundary condition - propagator combination. 4th order propagator can be used only with periodic boundary conditions.\n");
    return 0;
  }
  
  wf = (wf2d *) malloc(sizeof(wf2d));
  if (!wf) {
    fprintf(stderr, "libgrid: Error in grid2d_wf_alloc(). Could not allocate memory for wf2d.\n");
    return 0;
  }
  
  value_outside = NULL;
  if (boundary == WF2D_DIRICHLET_BOUNDARY)
    value_outside = cgrid2d_value_outside_constantdirichlet;
  else if (boundary == WF2D_NEUMANN_BOUNDARY)
    value_outside = cgrid2d_value_outside_neumann;
  else if (boundary == WF2D_PERIODIC_BOUNDARY)
    value_outside = cgrid2d_value_outside_periodic;
  
  wf->grid = cgrid2d_alloc(nx, ny, step, value_outside, 0);
  
  if (!wf->grid) {
    fprintf(stderr, "libgrid: Error in grid2d_wf_alloc(). Could not allocate memory for wf2d->grid.\n");
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
 * Free 2D wavefunction.
 *
 * gwf = wavefunction to be freed (wf2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_free(wf2d *gwf) {

  if (gwf) {
    if (gwf->grid) cgrid2d_free(gwf->grid);
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
 * potential   = absorbing potential (cgrid2d *; output; overwritten).
 * density     = current density (rgrid2d *; input).
 * rho0        = desired density in the boundary region (double, input).
 * region      = function that will multiply the absorbing potential (double (*)(double, double); input).
 *               0.0 = no absorption, \approx 0.1 = full absorption (values this
 *               large often makes the time propagation unstable).
 *               Even a linear function starting from the absorbing region edge (value 0) 
 *               and then gradually increasing up to 0.1 will often work well.
 * workspace   = temporary space needed for the operation (rgrid2d *).
 *
 * The calculated potential should be added to the potential that will be
 * propagated. Note that it is important that this is also included in the
 * possible predict-correct cycle as that improves the numericla stability.
 *
 */

EXPORT void grid2d_wf_absorb(cgrid2d *potential, rgrid2d *density, double rho0, double (*region)(void *, double, double), rgrid2d *workspace) {

  rgrid2d_copy(workspace, density);
  rgrid2d_add(workspace, -rho0);
  rgrid2d_product_func(workspace, region, NULL);
  rgrid2d_multiply(workspace, -1.0);
  grid2d_add_real_to_complex_im(potential, workspace);
}

/*
 * Calculate the x component of the probability flux.
 * 
 * gwf        = wavefunction for the operation (wf2d *).
 * flux_x     = x output frid containing the flux (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_probability_flux_x(const wf2d *gwf, rgrid2d *flux_x) {

#if 0
  // old code
  cgrid2d_fd_gradient_x(gwf->grid, workspace);
  cgrid2d_conjugate_product(workspace, gwf->grid, workspace);
  grid2d_complex_im_to_real(flux_x, workspace);
  rgrid2d_multiply(flux_x, HBAR / gwf->mass);
#endif
  // new code without additional workspace
  long i, j, ij;
  cgrid2d *grid = gwf->grid;
  long ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step), mass = gwf->mass;
  double *lvalue = flux_x->value;
  double complex (*value_at)(const cgrid2d *grid, long i, long j) = cgrid2d_value_at_index;
  
#pragma omp parallel for firstprivate(value_at,ny,nxy,lvalue,inv_delta,grid,mass) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = cimag(conj(value_at(grid, i, j)) * inv_delta * (value_at(grid, i+1, j) - value_at(grid, i-1, j))) * (HBAR / mass);
  }
}

/*
 * Calculate the y component of the probability flux.
 * 
 * gwf        = wavefunction for the operation (wf2d *).
 * flux_y     = y output frid containing the flux (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_probability_flux_y(const wf2d *gwf, rgrid2d *flux_y) {

#if 0
  // old code
  cgrid2d_fd_gradient_x(gwf->grid, workspace);
  cgrid2d_conjugate_product(workspace, gwf->grid, workspace);
  grid2d_complex_im_to_real(flux_y, workspace);
  rgrid2d_multiply(flux_y, HBAR / gwf->mass);
#endif
  // new code without additional workspace
  long i, j, ij;
  cgrid2d *grid = gwf->grid;
  long ny = grid->ny;
  long nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step), mass = gwf->mass;
  double *lvalue = flux_y->value;
  double complex (*value_at)(const cgrid2d *grid, long i, long j) = cgrid2d_value_at_index;
  
  /* main loop */
#pragma omp parallel for firstprivate(value_at,ny,nxy,lvalue,inv_delta,grid,mass) private(ij,i,j) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    lvalue[ij] = cimag(conj(value_at(grid, i, j)) * inv_delta * (value_at(grid, i, j+1) - value_at(grid, i, j-1))) * (HBAR / mass);
  }
}

/*
 * Calculate the probability flux.
 *
 * gwf        = wavefunction for the operation (wf2d *).
 * flux_x     = x output grid containing the flux (rgrid2d *).
 * flux_y     = y output grid containing the flux (rgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_probability_flux(const wf2d *gwf, rgrid2d *flux_x, rgrid2d *flux_y) {
  
  /*
   * J(r) = -i (hbar/2m) (psi^* grad psi - psi grad psi^*)
   *      = (hbar/m) Im[ psi^* grad psi ] 
   */
  grid2d_wf_probability_flux_x(gwf, flux_x);
  grid2d_wf_probability_flux_y(gwf, flux_y);
}

/*
 * Calculate the x component of the momentum.
 * 
 * gwf        = wavefunction for the operation (wf2d *).
 * momentum_x = output grid containing the momentum (cgrid2d *).
 * workspace  = additional storage needed for the operation (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_momentum_x(const wf2d *gwf, cgrid2d *momentum_x, cgrid2d *workspace) {

  cgrid2d_copy(workspace, gwf->grid);
  cgrid2d_fft(workspace);
  cgrid2d_fft_gradient_x(workspace, momentum_x);
  cgrid2d_inverse_fft(momentum_x);
  cgrid2d_multiply(momentum_x, -I*HBAR / (2.0 * gwf->mass));
}

/*
 * Calculate the y component of the momentum.
 * 
 * gwf        = wavefunction for the operation (wf2d *).
 * momentum_y = output frid containing the momentum (cgrid2d *).
 * workspace  = additional storage needed for the operation (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_momentum_y(const wf2d *gwf, cgrid2d *momentum_y, cgrid2d *workspace) {

  cgrid2d_copy(workspace, gwf->grid);
  cgrid2d_fft(workspace);
  cgrid2d_fft_gradient_y(workspace, momentum_y);
  cgrid2d_inverse_fft(momentum_y);
  cgrid2d_multiply(momentum_y, -I*HBAR / (2.0 * gwf->mass));
}

/*
 * Calculate the momentum.
 *
 * gwf        = wavefunction for the operation (wf2d *).
 * momentum_x = x output grid containing the momentum (cgrid2d *).
 * momentum_y = y output grid containing the momentum (cgrid2d *).
 * workspace  = additional storage needed for the operation (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_momentum(const wf2d *gwf, cgrid2d *momentum_x, cgrid2d *momentum_y, cgrid2d *workspace) {
  
  grid2d_wf_momentum_x(gwf, momentum_x, workspace);
  grid2d_wf_momentum_y(gwf, momentum_y, workspace);
}

/*
 * Calculate energy for the wavefunction.
 *
 * gwf       = wavefunction for the energy calculation (wf2d *).
 * potential = grid containing the potential (cgrid2d *).
 * workspace = additional storage needed (cgrid2d *).
 *
 * Returns the energy (double).
 *
 */

EXPORT double grid2d_wf_energy(const wf2d *gwf, const cgrid2d *potential, cgrid2d *workspace) {

  if (gwf->boundary == WF2D_DIRICHLET_BOUNDARY || gwf->boundary == WF2D_NEUMANN_BOUNDARY)
    return grid2d_wf_energy_cn(gwf, gwf, potential, workspace);
  else if (gwf->boundary == WF2D_PERIODIC_BOUNDARY)
    return grid2d_wf_energy_fft(gwf, potential, workspace);
  else
    abort();
}

/*
 * Auxiliary routine for calculating the energy (Crank-Nicolson).
 * Users should rather call grid2d_wf_energy().
 *
 * gwfa      = (left) wavefunction for the energy calculation (wf2d *).
 * gwfb      = (right) wavefunction for the energy calculation (wf2d *).
 *             Normally gwfa = gwfb.
 * potential = Potential grid (cgrid2d *).
 * workspace = Additional workspace needed for the operation (cgrid2d *).
 * 
 * Returns the energy.
 *
 */

EXPORT double grid2d_wf_energy_cn(const wf2d *gwfa, const wf2d *gwfb, const cgrid2d *potential, cgrid2d *workspace) {
  
  long i;

  /* (-2m/hbar^2) T psi */
  cgrid2d_fd_laplace(gwfb->grid, workspace);
  cgrid2d_multiply(workspace, -HBAR * HBAR / (2.0 * gwfb->mass));

  /* V psi */
  if(potential) 
    for(i = 0; i < gwfb->grid->nx * gwfb->grid->ny; i++)
      workspace->value[i] += potential->value[i] * gwfb->grid->value[i];
  
  /* int psi^* (T + V) psi d^3r */
  return creal(cgrid2d_integral_of_conjugate_product(gwfa->grid, workspace));
}

/*
 * Auxiliary routine for calculating the energy (FFT).
 * Users should rather call grid2d_wf_energy().
 *
 * gwf       = wavefunction for the energy calculation (wf2d *).
 * potential = Potential grid (cgrid3d *).
 * workspace = Additional workspace needed for the operation (cgrid2d *).
 * 
 * Returns the energy.
 *
 */

EXPORT double grid2d_wf_energy_fft(const wf2d *gwf, const cgrid2d *potential, cgrid2d *workspace) {

  double en;

  /* delta (- k^2) fft[f(x)] / N */
  cgrid2d_copy(workspace, gwf->grid);
  cgrid2d_fft(workspace);
  
  en = -HBAR*HBAR / (2.0 * gwf->mass) * cgrid2d_fft_laplace_expectation_value(workspace, workspace);
  if(potential) en += creal(cgrid2d_grid_expectation_value(gwf->grid, potential));
  return en;
}

/*
 * Auxiliary routine for calculating kinetic energy (FFT).
 * This is used by grid2d_wf_energy_fft().
 * 
 * gwf       = wavefunction for the kinetic energy calculation (wf2d *).
 * workspace = additional workspace required for the operation (cgrid2d *).
 *
 * Returns the kinetic energy.
 *
 */

EXPORT double grid2d_wf_kinetic_energy_fft(const wf2d *gwf, cgrid2d *workspace) {

  /* delta (- k^2) fft[f(x)] / N */
  cgrid2d_copy(workspace, gwf->grid);
  cgrid2d_fft(workspace);
  
  return -HBAR*HBAR / (2.0 * gwf->mass) * cgrid2d_fft_laplace_expectation_value(workspace, workspace);
}

/*
 * Auxiliary routine for calculating potential energy.
 * 
 * gwf       = wavefunction for the kinetic energy calculation (wf2d *).
 * workspace = additional workspace required for the operation (cgrid2d *).
 *
 * Returns the potential energy.
 *
 */

EXPORT double grid2d_wf_potential_energy(const wf2d *gwf, const cgrid2d *potential) {

  return creal(cgrid2d_grid_expectation_value(gwf->grid, potential));
}

/*
 * Calculate energy (E) and the error (dE).
 *
 * gwf       = wavefunction for the operation (wf2d *).
 * potential = grid containing the potential (cgrid2d *).
 * workspace = additional workspace required for the operation (cgrid2d *).
 * error     = error estimate for energy (double *).
 * 
 * Returns the energy.
 *
 */

EXPORT double grid2d_wf_energy_and_error(const wf2d *gwf, const cgrid2d *potential, cgrid2d *workspace, double *error) { 

  double energy;
  
/* 
 * energy and its error
 * dE^2 = int dE^2 |psi|^2 dtau = ... = int E_local^2 |psi|^2 dtau - E_avg^2
 * dE = E_local - E_avg
 * E_local psi = H psi
 */

  /* T psi */
  if (gwf->boundary == WF2D_DIRICHLET_BOUNDARY || gwf->boundary == WF2D_NEUMANN_BOUNDARY) {
    cgrid2d_fd_laplace(gwf->grid, workspace);
    cgrid2d_multiply(workspace, -HBAR*HBAR / (2.0 * gwf->mass));
  }
  else if (gwf->boundary == WF2D_PERIODIC_BOUNDARY) {
    cgrid2d_copy(workspace, gwf->grid);
    cgrid2d_fft(workspace);
    cgrid2d_fft_laplace(workspace, workspace);
    cgrid2d_scaled_inverse_fft(workspace, -HBAR*HBAR / (2.0 * gwf->mass));
  }
  else {
    fprintf(stderr, "libgrid: Error in grid2d_wf_energy_and_error(). Invalid boundary condition, index = %d\n", gwf->boundary);
    abort();
  }
  
  /* H psi */
  cgrid2d_add_scaled_product(workspace, 1.0, potential, gwf->grid);
  
  /* int E_local^2 |psi|^2 dtau */
  *error = cgrid2d_integral_of_square(workspace);
  
  /* int E_local |psi|^2 dtau */
  cgrid2d_conjugate_product(workspace, gwf->grid, workspace);
  energy = creal(cgrid2d_integral(workspace));
  
  /* sqrt( int E_local^2 |psi|^2 dtau - ( int E_local |psi|^2 dtau )^2 ) */
  *error = sqrt(*error - energy * energy);
  
  return energy;
}

/*
 * Propagate wavefunction in time subject to given potential.
 *
 * gwf         = wavefunction to be propagated (wf2d *).
 * potential   = grid containing the potential (cgrid2d *).
 * sq_grad_pot = grid containing square of potential gradient (cgrid2d *).
 * time        = time step (double complex). Note this may be either real or imaginary.
 * workspace   = additional workspace needed for the operation (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate(wf2d *gwf, const cgrid2d *potential, const cgrid2d *sq_grad_pot, double complex time, cgrid2d *workspace) {  
  
  double complex half_time = 0.5 * time;
  double complex one_sixth_time = time / 6.0;
  double complex two_thirds_time = 2.0 * time / 3.0;
  
  if ((gwf->boundary == WF2D_DIRICHLET_BOUNDARY || gwf->boundary == WF2D_NEUMANN_BOUNDARY)
       && gwf->propagator == WF2D_2ND_ORDER_PROPAGATOR) {    
    grid2d_wf_propagate_potential(gwf, potential, half_time);
    grid2d_wf_propagate_kinetic_cn(gwf, time, workspace);
    grid2d_wf_propagate_potential(gwf, potential, half_time);
  } else if (gwf->boundary == WF2D_PERIODIC_BOUNDARY
	     && gwf->propagator == WF2D_2ND_ORDER_PROPAGATOR) {
    grid2d_wf_propagate_potential(gwf, potential, half_time);
    grid2d_wf_propagate_kinetic_fft(gwf, time);
    grid2d_wf_propagate_potential(gwf, potential, half_time);
  } else if (gwf->boundary == WF2D_PERIODIC_BOUNDARY 
            && gwf->propagator == WF2D_4TH_ORDER_PROPAGATOR) {    
    grid2d_wf_propagate_potential(gwf, potential, one_sixth_time);
    grid2d_wf_propagate_kinetic_fft(gwf, half_time);
    
    cgrid2d_copy(workspace, potential);
    cgrid2d_add_scaled(workspace, (1/48.0 * HBAR * HBAR / gwf->mass) * sqnorm(time), sq_grad_pot);	
    grid2d_wf_propagate_potential(gwf, workspace, two_thirds_time);
    
    grid2d_wf_propagate_kinetic_fft(gwf, half_time);
    grid2d_wf_propagate_potential(gwf, potential, one_sixth_time);
  } else {
    fprintf(stderr, "libgrid: Error in grid2d_wf_propagate(). Unknown propagator - boundary value combination (propagator index = %d, boundary index = %d).\n", gwf->propagator, gwf->boundary);
    abort();
  }
}

/*
 * Auxiliary routine to propagate kinetic energy using FFT.
 * Users should call grid2d_wf_propaget().
 *
 * gwf  = wavefunction to be propagated (wf2d *).
 * time = time step (double complex).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_fft(wf2d *gwf, double complex time) {

  long i, j, ij, nx, ny, nxy;
  double kx, ky, step, norm, mass;
  double complex *value = gwf->grid->value;
  
  cgrid2d_fft(gwf->grid);
  
  nx = gwf->grid->nx;
  ny = gwf->grid->ny;
  nxy = nx * ny;
  step = gwf->grid->step;
  mass = gwf->mass;
  
  /* f(x) = ifft[fft[f(x)]] / N */
  norm = (1.0 / nx) * (1.0 / ny);
  
  /*cTimer_start(&timer);*/
  
#pragma omp parallel for firstprivate(norm,nx,ny,nxy,step,value,time,mass) private(i,j,ij,kx,ky) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    i = ij / ny;
    j = ij % ny;
    
    /* 
     * k = 2 pi n / L 
     * if k < n/2, k = k
     * else k = -k
     */
    if (i < nx / 2)
      kx = 2.0 * M_PI * i / (nx * step);
    else 
      kx = 2.0 * M_PI * (i - nx) / (nx * step);
    
    if (j < ny / 2)
      ky = 2.0 * M_PI * j / (ny * step);
    else 
      ky = 2.0 * M_PI * (j - ny) / (ny * step);
    
    /* psi(t+dt) = psi(t) exp( - i (hbar^2 * k^2 / 2m) dt / hbar ) */	  
    value[ij] *= norm * cexp(-I * time * HBAR * HBAR * (kx*kx + ky*ky) / (2 * mass * HBAR));
  }
  
  cgrid2d_inverse_fft(gwf->grid);
}

/*
 * Auxiliary routine to propagate kinetic energy (Crank-Nicolson).
 * Users should call grid2d_wf_propaget().
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional workspace required for the operation (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_cn(wf2d *gwf, double complex time, cgrid2d *workspace) {

  if (gwf->boundary == WF2D_DIRICHLET_BOUNDARY)
    grid2d_wf_propagate_kinetic_cn_dbc(gwf, time, workspace);
  else if (gwf->boundary == WF2D_NEUMANN_BOUNDARY)
    grid2d_wf_propagate_kinetic_cn_nbc(gwf, time, workspace);
  else if (gwf->boundary == WF2D_PERIODIC_BOUNDARY)
    grid2d_wf_propagate_kinetic_cn_pbc(gwf, time, workspace);
  else {
    fprintf(stderr, "libgrid: Error in grid2d_wf_propagate_cn(). Unknown boundary condition (index = %d).\n", gwf->boundary);
    abort();
  }
}

/*
 * Auxiliary routine to propagate potential energy.
 * Users should call grid2d_wf_propaget().
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * potential = grid containing the potential (cgrid2d *).
 * time      = time step (double complex).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_potential(wf2d *gwf, const cgrid2d *potential, double complex time) {

  long ij, nxy = gwf->grid->nx * gwf->grid->ny;
  double complex c, *psi = gwf->grid->value, *pot = potential->value;
  
  /* c = - i dt / hbar */
  c = -I * time / HBAR;
  
  #pragma omp parallel for firstprivate(nxy,psi,pot,c) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    /* psi(t+dt) = exp(- i V dt / hbar) psi(t) */
    psi[ij] *= cexp(c * pot[ij]);
}

/*
 * Project "gwfb" out from "gwfa".
 *
 * gwfa = input wavefunction (wf2d *).
 * gwfb = this will be projected out from gwfa (wf2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_project_out(wf2d *gwfa, const wf2d *gwfb) {

  double complex overlap = grid2d_wf_overlap(gwfa, gwfb);

  cgrid2d_add_scaled(gwfa->grid, -overlap, gwfb->grid);
}

/*
 * "traditional" diagonalization of Hamiltonian.
 *
 * gwf    = an array of wavefunctions (wf2d **).
 * states = number of states (int).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_diagonalize(wf2d **gwf, int states) {

  int i,j;
  double *eigenvalue = (double *) malloc(states * sizeof(double));
  double complex *overlap = (double complex *) malloc(states * states * sizeof(double complex));
  wf2d *gwf_tmp;
  
  if (states == 1) {
    grid2d_wf_normalize(gwf[0]);
    return;
  }
  
  /* overlap matrix */
  for(i = 0; i < states; i++) {
    for(j = 0; j <= i; j++) {
      /* fortran (column major) matrix order, i is row (minor) index, j is column (major) index */
      overlap[i + j * states] = grid2d_wf_overlap(gwf[i], gwf[j]);
      overlap[j + i * states] = conj(overlap[i + j * states]);
    }
  }
  
  /* diagonalize */
  grid_hermitian_eigenvalue_problem(eigenvalue, overlap, states);
  
  /* phi_i = 1 / sqrt(m_i) C_ij psi_j, C (row major) matrix order ???is it??? */
  for(i = 0; i < states; i++)
    for(j = 0; j < states; j++)
      overlap[i * states + j] /= sqrt(eigenvalue[i]);
  
  grid2d_wf_linear_transform(gwf, overlap, states);
  
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
 * gwf       = an array of wavefunctions (wf2d **).
 * transform = transformation matrix (double complex *).
 * states    = number of states (int).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_linear_transform(wf2d **gwf, double complex *transform, int states) {

  int p, q, offset;
  long ij, nxy;
  double complex **value, *tmp;
  
  nxy = gwf[0]->grid->nx * gwf[0]->grid->ny;
  
  /* + 16 to prevent "write locks" */
  tmp = (double complex *) malloc(omp_get_max_threads() * (states + 16) * sizeof(double complex));
  
  value = (double complex **) malloc(states * sizeof(double complex *));

  for(p = 0; p < states; p++)
    value[p] = gwf[p]->grid->value;
  
#pragma omp parallel for firstprivate(nxy,value,states,transform,tmp) private(ij,p,q,offset) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    
    offset = (states + 16) * omp_get_thread_num();
    
    for(p = 0; p < states; p++)
      tmp[offset + p] = 0.0;
    
    for(p = 0; p < states; p++)
      for(q = 0; q < states; q++)
	tmp[offset + p] += transform[p * states + q] * value[q][ij];
    
    for(p = 0; p < states; p++)
      value[p][ij] = tmp[offset + p];
  }
  
  free(value);
  free(tmp);
}

/*
 * Calcuate square of potential gradient.
 *
 * sq_grad_pot = output grid (cgrid2d *).
 * potential   = potental input grid (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_square_of_potential_gradient(cgrid2d *sq_grad_pot, const cgrid2d *potential, cgrid2d *workspace) {

  cgrid2d_copy(sq_grad_pot, potential);
  cgrid2d_fft(sq_grad_pot);
  cgrid2d_fft_gradient(sq_grad_pot, sq_grad_pot, workspace);
  
  cgrid2d_inverse_fft(sq_grad_pot);
  cgrid2d_inverse_fft(workspace);
  
  cgrid2d_conjugate_product(sq_grad_pot, sq_grad_pot, sq_grad_pot);
  cgrid2d_conjugate_product(workspace, workspace, workspace);
  
  cgrid2d_sum(sq_grad_pot, sq_grad_pot, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy using Dirichlet boundaries (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_cn_dbc(wf2d *gwf, double complex time, cgrid2d *workspace) {

  /*
   * exp(-i (Tx + Ty) dt / hbar) 
   *   = exp(-i Tx dt / hbar) exp(-i Ty dt / hbar) + O(dt^2)
   *   
   */
  
  grid2d_wf_propagate_kinetic_x_cn_dbc(gwf, time, workspace);
  grid2d_wf_propagate_kinetic_y_cn_dbc(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along x using Dirichlet boundaries (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_x_cn_dbc(wf2d *gwf, double complex time, cgrid2d *workspace) {

  double complex c, ca, boundary, *ad, *ld, *b, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / hbar dt)
   */
  c = I * 4.0 * gwf->mass * step * step / HBAR / time;

  /* create cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;
  ad = &wrk[omp_get_max_threads() * nx];
  
  /* create matrix */
  for(i = 0; i < nx; i++)
    ad[i] = ca;  
  
  /* create decomposition */
  grid_cholesky_decomposition(ad, nx, 0);
  /* matrix is replaced by decomposition */
  ld = ad;

  /* solve problem for each y */
#pragma omp parallel for firstprivate(c,psi,nx,ny,ld,wrk) private(i,j,ind,tid,b,boundary) default(none) schedule(runtime)
  for(j = 0; j < ny; j++) {
    tid = omp_get_thread_num();
      
    /* use diffent workspace for different threads  */
    b = &wrk[tid * nx];
      
    /* create right-hand side vector */
    for(i = 1; i < nx - 1; i++) {
      ind = i * ny + j;
      /* (C - dx^2 laplace) psi(t) */
      b[i] = c * psi[ind] - (psi[ind + ny] - 2.0 * psi[ind] + psi[ind - ny]);
    }
      
    /* dirichlet boundaries */
    boundary = 0.0;
    i = 0; 
    ind = i * ny + j;
    b[i] = c * psi[ind] - (psi[ind + ny] - 2.0 * psi[ind] + boundary);
    i = nx - 1;
    ind = i * ny + j;
    b[i] = c * psi[ind] - (boundary - 2.0 * psi[ind] + psi[ind - ny]);
    
    /* solve */
    grid_cholesky_substitute(ld, b, nx, 0);
    
    /* copy */
    for(i = 0; i < nx; i++)
      psi[i * ny + j] = b[i];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along y using Dirichlet boundaries (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_y_cn_dbc(wf2d *gwf, double complex time, cgrid2d *workspace) {

  double complex c, ca, boundary, *ad, *ld, *b, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind, iny;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / hbar dt)
   */
  c = I * (4.0 * gwf->mass * step * step / HBAR) / time;
  
  /* create cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;
  ad = &wrk[omp_get_max_threads() * ny];
  
  /* create matrix */
  for(j = 0; j < ny; j++)
    ad[j] = ca;  
  
  /* create decomposition */
  grid_cholesky_decomposition(ad, ny, 0);
  /* matrix is replaced by decomposition */
  ld = ad;
  
  /* solve problem for each x */
#pragma omp parallel for firstprivate(c,psi,nx,ny,ld,wrk) private(i,j,ind,iny,tid,b,boundary) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {
    iny = i * ny;

    tid = omp_get_thread_num();
      
    /* use diffent workspace for different threads  */
    b = &wrk[tid * ny];
      
    /* create right-hand side vector */
    for(j = 1; j < ny - 1; j++) {
      ind = iny + j;
      /* (C - dx^2 laplace) psi(t) */
      b[j] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + psi[ind - 1]);
    }
      
    /* dirichlet boundaries */
    boundary = 0.0;
    j = 0;
    ind = iny + j;
    b[j] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + boundary);
    j = ny - 1;
    ind = iny + j;
    b[j] = c * psi[ind] - (boundary - 2.0 * psi[ind] + psi[ind - 1]);
    
    /* solve */
    grid_cholesky_substitute(ld, b, ny, 0);
    
    /* copy */
    for(j = 0; j < ny; j++)
      psi[iny + j] = b[j];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 *
 * FIXME: This has still issues with boundary condition (i.e., private.h solver is still messed up).
 *        One could switch to lapack just like in *_cyl.c
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_cn_nbc(wf2d *gwf, double complex time, cgrid2d *workspace) {

  /*
   * exp(-i (Tx + Ty) dt / hbar) 
   *   = exp(-i Tx dt / hbar) exp(-i Ty dt / hbar) + O(dt^2)
   *   
   */
  
  grid2d_wf_propagate_kinetic_x_cn_nbc(gwf, time, workspace);
  grid2d_wf_propagate_kinetic_y_cn_nbc(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along x using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_x_cn_nbc(wf2d *gwf, double complex time, cgrid2d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  double complex *du, *dl, *d, *b;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* create cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nx,ny,psi,wrk,c) private(tid,i,j,dl,d,du,b,ind) default(none) schedule(runtime)
  for (j = 0; j < ny; j++) {  /* for each y */
    tid = omp_get_thread_num();
    dl = &wrk[nx * (4 * tid + 0)];
    d  = &wrk[nx * (4 * tid + 1)];
    du = &wrk[nx * (4 * tid + 2)];
    b  = &wrk[nx * (4 * tid + 3)];

    for(i = 0; i < nx; i++) {
      
      if(i == nx-1) dl[nx-1] = 2.0;              /* Neumann outer boundary */
      else if(i > 0) dl[i] = 1.0;
      
      if(i == 0) du[0] = 2.0;                   /* Neumann inner boundary */
      else if(i < nx-1) du[i] = 1.0;
      
      d[i] = ca;
    }

    /* create right-hand side vector */
    for(i = 1; i < nx - 1; i++) {
      ind = i * ny + j;
      /* (C - dx^2 laplace) psi(t) */
      b[i] = c * psi[ind] - (psi[ind + ny] - 2.0 * psi[ind] + psi[ind - ny]);
    }
 
    /* neumann boundaries */
    i = 0; 
    ind = i * ny + j;
    b[i] = c * psi[ind] - (2.0 * psi[ind + ny] - 2.0 * psi[ind]);

    i = nx - 1;
    ind = i * ny + j;
    b[i] = c * psi[ind] - (2.0 * psi[ind - ny] - 2.0 * psi[ind]);
    
    grid_solve_tridiagonal_system(nx, dl, d, du, b, dl);

    for(i = 0; i < nx; i++)
      psi[i * ny + j] = dl[i];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along y using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_y_cn_nbc(wf2d *gwf, double complex time, cgrid2d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind, iny;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  double complex *dl, *d, *du, *b;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = 4 i m dx^2 / (hbar dt)
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);
  
  /* create cholesky decompostions for kinetic x */
  
  /* (C + dr^2 laplace + (1/r) dr^2 grad) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nx,ny,psi,wrk,c,step) private(tid,i,j,dl,d,du,b,ind,iny) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) { /* for each x */
    tid = omp_get_thread_num();
    dl = &wrk[nx * (4 * tid + 0)];
    d  = &wrk[nx * (4 * tid + 1)];
    du = &wrk[nx * (4 * tid + 2)];
    b  = &wrk[nx * (4 * tid + 3)];

    for(j = 0; j < ny; j++) { /* over y */

      if(j == ny-1) dl[ny-1] = 2.0;            /* Neumann outer boundary */
      else if(j > 0) dl[j] = 1.0; 
      
      if(j == 0) du[0] = 2.0;                  /* Neumann inner boundary */
      else if(j < ny-1) du[j] = 1.0;
      
      d[j] = ca;
    }

    /* create right-hand side vector */
    iny = i * ny;
    for(j = 1; j < ny - 1; j++) {
      ind = iny + j;
      /* (C - dr^2 laplace - (1/r) dr^2 grad) psi(t) */
      b[j] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + psi[ind - 1]);
    }
    
    /* neumann boundaries */
    j = 0; 
    ind = iny + j;
    b[j] = c * psi[ind] - (2.0 * psi[ind + 1] - 2.0 * psi[ind]);
 
    j = ny - 1; 
    ind = iny + j;
    b[j] = c * psi[ind] - (2.0 * psi[ind - 1] - 2.0 * psi[ind]);

    grid_solve_tridiagonal_system(ny, dl, d, du, b, dl);

    /* copy */
    for(j = 0; j < ny; j++)
      psi[iny + j] = dl[j];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy using periodic boundary condition (Crank-Nicolson).
 * 
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 * 
 */

EXPORT void grid2d_wf_propagate_kinetic_cn_pbc(wf2d *gwf, double complex time, cgrid2d *workspace) {

  /*
   * exp( -i (Tx + Ty) dt / hbar ) 
   *   = exp( -i Tx dt / hbar ) exp( -i Ty dt / hbar ) exp( -i Tz dt / hbar ) + O(dt^2)
   *   
   */

  grid2d_wf_propagate_kinetic_x_cn_pbc(gwf, time, workspace);
  grid2d_wf_propagate_kinetic_y_cn_pbc(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along x using periodic boundary condition (Crank-Nicolson).
 * 
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 * 
 */

EXPORT void grid2d_wf_propagate_kinetic_x_cn_pbc(wf2d *gwf, double complex time, cgrid2d *workspace) {

  double complex c, ca, boundary, *ad, *ld, *b, *w, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / hbar dt)
   */
  c = I * (4.0 * gwf->mass * step * step / HBAR) / time;
  
  /* create cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;
  ad = &wrk[omp_get_max_threads() * nx];
  
  /* create matrix A' */
  for(i = 0; i < nx; i++)
    ad[i] = ca;  
  /* periodic boundaries */
  ad[0]    -= 1.0;
  ad[nx-1] -= 1.0;
  
  /* create decomposition */
  grid_cholesky_decomposition(ad, nx, 0);
  /* matrix is replaced by decomposition */
  ld = ad;
  
  /* periodic boundaries */
  w = &wrk[omp_get_max_threads() * nx + nx];
  for(i = 0; i < nx; i++) 
    w[i] = 0.0;
  w[0] = w[nx-1] = 1.0;
  
  /* solve A' z = u */
  grid_cholesky_substitute(ld, w, nx, 0);
  
  /* solve problem for each y and z */
#pragma omp parallel for firstprivate(c,psi,nx,ny,ld,w,wrk) private(i,j,ind,tid,ca,b,boundary) default(none) schedule(runtime)
  for(j = 0; j < ny; j++) {
    tid = omp_get_thread_num();
      
    /* use diffent workspace for different threads */
    b = &wrk[tid * nx];
    
    /* create right-hand side vector */
    for(i = 1; i < nx-1; i++) {
      ind = i * ny + j;
      /* (C - dx^2 laplace) psi(t) */
      b[i] = c * psi[ind] - (psi[ind + ny] - 2.0 * psi[ind] + psi[ind - ny]);
    }
    
    /* periodic boundaries */
    i = nx-1; boundary = psi[i * ny + j];
    i = 0; ind = i * ny + j;
    b[i] = c * psi[ind] - (psi[ind + ny] - 2.0 * psi[ind] + boundary);
    
    i = 0; boundary = psi[i * ny + j];
    i = nx-1; ind = i * ny + j;
    b[i] = c * psi[ind] - (boundary - 2.0 * psi[ind] + psi[ind - ny]);
    
    /* solve A' x' = b */
    grid_cholesky_substitute(ld, b, nx, 0);  

    /* Shermann-Morrison form NR in C
     *             v . y'                y'_0 + y'_(n-1)
     * y = y' - ----------- w = y' - -------------------- z
     *           1 + v . w             1 + w_0 + w_(n-1)
     *
     * where v = (1, 0, 0, 0, ..., 0, 1) 
     */	  
    ca = (b[0] + b[nx-1]) / (1.0 + w[0] + w[nx-1]);
    for(i = 0; i < nx; i++)
      b[i] -= ca * w[i];
    
    /* copy back */
    for(i = 0; i < nx; i++)
      psi[i * ny + j] = b[i];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along y using periodic boundary condition (Crank-Nicolson).
 * 
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 * 
 */

EXPORT void grid2d_wf_propagate_kinetic_y_cn_pbc(wf2d *gwf, double complex time, cgrid2d *workspace) {

  double complex c, ca, boundary, *ad, *ld, *b, *w, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind, iny;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / hbar dt)
   */
  c = I * (4.0 * gwf->mass * step * step / HBAR) / time;
  
  /* create cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;
  ad = &wrk[omp_get_max_threads() * ny];
  
  /* create matrix A' */
  for(j = 0; j < ny; j++)
    ad[j] = ca;  
  /* periodic boundaries */
  ad[0] -= 1.0;
  ad[ny-1] -= 1.0;
  
  /* create decomposition */
  grid_cholesky_decomposition(ad, ny, 0);
  /* matrix is replaced by decomposition */
  ld = ad;  
  
  /* periodic boundaries */
  w = &wrk[omp_get_max_threads() * ny + ny];
  for(i = 0; i < ny; i++) 
    w[i] = 0.0;
  w[0] = w[ny-1] = 1.0;
  
  /* solve A' z = u */
  grid_cholesky_substitute(ld, w, ny, 0);
  
  /* solve problem for each x */
#pragma omp parallel for firstprivate(c,psi,nx,ny,ld,wrk,w) private(i,j,ind,iny,tid,ca,b,boundary) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {
    iny = i * ny;

    tid = omp_get_thread_num();
      
    /* use diffent workspace for different threads  */
    b = &wrk[tid * ny];
    
    /* create right-hand side vector */
    for(j = 1; j < ny-1; j++) {
      ind = iny + j;
      /* (C - dx^2 laplace) psi(t) */
      b[j] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + psi[ind - 1]);
    }
      
    /* periodic boundaries */
    j = ny-1; boundary = psi[i * ny + j];
    j = 0; ind = iny + j;
    b[j] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + boundary);
    
    j = 0; boundary = psi[i * ny + j];
    j = ny-1; ind = iny + j;
    b[j] = c * psi[ind] - (boundary - 2.0 * psi[ind] + psi[ind - 1]);
    
    /* solve A' y' = b */
    grid_cholesky_substitute(ld, b, ny, 0);

    /* Shermann-Morrison form NR in C
     *             v . y'                y'_0 + y'_(n-1)
     * y = y' - ----------- w = y' - -------------------- z
     *           1 + v . w             1 + w_0 + w_(n-1)
     *
     * where v = (1, 0, 0, 0, ..., 0, 1) 
     */	  
    ca = (b[0] + b[ny-1]) / (1.0 + w[0] + w[ny-1]);
    for(i = 0; i < ny; i++)
      b[i] -= ca * w[i];

    /* copy back */
    for(j = 0; j < ny; j++)
      psi[iny + j] = b[j];
  }
}

/*
 * Produce density grid from a given wavefunction.
 *
 * gwf     = wavefunction (wf2d *).
 * density = output density grid (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT inline void grid2d_wf_density(const wf2d *gwf, rgrid2d *density) {

  long ij, nxy = gwf->grid->nx * gwf->grid->ny;
  double complex *avalue = gwf->grid->value;
  double *cvalue = density->value;
  
#pragma omp parallel for firstprivate(nxy,avalue,cvalue) private(ij) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++)
    cvalue[ij] = conj(avalue[ij]) * avalue[ij];
}

/*
 * Zero wavefunction.
 *
 * gwf = wavefunction to be zeroed (wf2d *).
 *
 * No return value.
 *
 */

EXPORT inline void grid2d_wf_zero(wf2d *gwf) { 

  cgrid2d_zero(gwf->grid); 
}

/*
 * Set wavefunction to some constant value.
 *
 * gwf = wavefunction to be set (wf2d *).
 * c   = value (double complex).
 *
 * No return value.
 *
 */

EXPORT inline void grid2d_wf_constant(wf2d *gwf, double complex c) { 

  cgrid2d_constant(gwf->grid, c); 
}

/*
 * Map a given function on a wavefunction.
 *
 * gwf  = wavefunction where function will be mapped to (wf2d *).
 * func = function providing the mapping (double complex (*)(void *, double)).
 * farg = optional argument for passing parameters to func (void *).
 *
 * No return value.
 *
 */

EXPORT inline void grid2d_wf_map(wf2d *gwf, double complex (*func)(void *arg, double x, double y), void *farg) { 

  cgrid2d_map(gwf->grid, func, farg); 
}

/*
 * Calculate the norm of the given wavefunction.
 *
 * gwf = wavefunction for the calculation (wf2d *).
 *
 * Returns the norm (double).
 *
 */

EXPORT inline double grid2d_wf_norm(const wf2d *gwf) { 

  return cgrid2d_integral_of_square(gwf->grid); 
}

/*
 * Normalize wavefunction (to the value given in gwf->norm).
 *
 * gwf = wavefunction to be normalized (wf2d *).
 *
 * Returns the normalization constant (double).
 *
 */

EXPORT inline double grid2d_wf_normalize(wf2d *gwf) { 

  double norm = grid2d_wf_norm(gwf);

  cgrid2d_multiply(gwf->grid, sqrt(gwf->norm / norm));
  return norm; 
}

/*
 * Calculate overlap between two wavefunctions.
 *
 * gwfa = 1st wavefunction (wf2d *).
 * gwfb = 2nd wavefunction (wf2d *).
 *
 * Returns the overlap (double complex).
 *
 */

EXPORT inline double complex grid2d_wf_overlap(const wf2d *gwfa, const wf2d *gwfb) { 

  return cgrid2d_integral_of_conjugate_product(gwfa->grid, gwfb->grid); 
}

/*
 * Output wavefunction.
 *
 * gwf = wavefunction to be printed (wf2d *).
 * out = output file pointer (FILE *).
 *
 * No return value.
 *
 */

EXPORT inline void grid2d_wf_print(const wf2d *gwf, FILE *out) { 

  cgrid2d_print(gwf->grid, out); 
}
