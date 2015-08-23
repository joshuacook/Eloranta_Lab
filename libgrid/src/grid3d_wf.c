/*
 * Routines for handling 3D wavefunctions.
 *
 */

#include "grid.h"
#include "private.h"
#include "private3d.h"

/*
 * Allocate a 3D wavefunction.
 *
 * nx         = number of spatial grid points along x (int).
 * ny         = number of spatial grid points along y (int).
 * nz         = number of spatial grid points along z (int).
 * step       = spatial step size (double).
 * mass       = mass of the particle corresponding to this wavefunction (double).
 * boundary   = boundary condition (int):
 *              WF3D_DIRICHLET_BOUNDARY = Dirichlet boundary condition.
 *              WF3D_NEUMANN_BOUNDARY   = Neumann boundary condition.
 *              WF3D_PERIODIC_BOUNDARY  = Periodic boundary condition.
 * propagator = which time propagator to use for this wavefunction:
 *              WF3D_2ND_ORDER_PROPAGATOR = 2nd order in time.
 *              WF3D_4TH_ORDER_PROPAGATOR = 4th order in time.
 *
 * Return value is a pointer to the allocated wavefunction.
 * This routine returns NULL if allocation fails.
 *
 */

EXPORT wf3d *grid3d_wf_alloc(long nx, long ny, long nz, double step, double mass, int boundary, int propagator) {

  wf3d *wf;
  double complex (*value_outside)(const struct cgrid3d_struct *grid, long i, long j, long k);
  
  if(boundary != WF3D_DIRICHLET_BOUNDARY 
     && boundary != WF3D_NEUMANN_BOUNDARY 
     && boundary != WF3D_PERIODIC_BOUNDARY
     && boundary != WF3D_VORTEX_X_BOUNDARY
     && boundary != WF3D_VORTEX_Y_BOUNDARY
     && boundary != WF3D_VORTEX_Z_BOUNDARY) {
    fprintf(stderr, "libgrid: Error in grid3d_wf_alloc(). Unknown boundary condition (index = %d).\n", boundary);
    return 0;
  }
  
  if(propagator != WF3D_2ND_ORDER_PROPAGATOR 
       && propagator != WF3D_4TH_ORDER_PROPAGATOR) {
    fprintf(stderr, "libgrid: Error in grid3d_wf_alloc(). Unknown propagator (index = %d).\n", propagator);
    return 0;
  }
  
  if ((boundary == WF3D_DIRICHLET_BOUNDARY || boundary == WF3D_NEUMANN_BOUNDARY)
      && propagator == WF3D_4TH_ORDER_PROPAGATOR) {
    fprintf(stderr, "libgrid: Error in grid3d_wf_alloc(). Invalid boundary condition - propagator combination. 4th order propagator can be used only with periodic boundary conditions.\n");
    return 0;
  }
  
  wf = (wf3d *) malloc(sizeof(wf3d));
  if (!wf) {
    fprintf(stderr, "libgrid: Error in grid3d_wf_alloc(). Could not allocate memory for wf3d.\n");
    return 0;
  }
  
  value_outside = NULL;
  if (boundary == WF3D_DIRICHLET_BOUNDARY)
    value_outside = cgrid3d_value_outside_constantdirichlet;
  else if (boundary == WF3D_NEUMANN_BOUNDARY)
    value_outside = cgrid3d_value_outside_neumann;
  else if (boundary == WF3D_PERIODIC_BOUNDARY)
    value_outside = cgrid3d_value_outside_periodic;
  else if (boundary == WF3D_VORTEX_X_BOUNDARY)
    value_outside = cgrid3d_value_outside_vortex_x;
  else if (boundary == WF3D_VORTEX_Y_BOUNDARY)
    value_outside = cgrid3d_value_outside_vortex_y;
  else if (boundary == WF3D_VORTEX_Z_BOUNDARY)
    value_outside = cgrid3d_value_outside_vortex_z;
  
  wf->grid = cgrid3d_alloc(nx, ny, nz, step, value_outside, 0);
  
  if (!wf->grid) {
    fprintf(stderr, "libgrid: Error in grid3d_wf_alloc(). Could not allocate memory for wf3d->grid.\n");
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
 * Free 3D wavefunction.
 *
 * gwf = wavefunction to be freed (wf3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_free(wf3d *gwf) {

  if (gwf) {
    if (gwf->grid) cgrid3d_free(gwf->grid);
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
 * potential   = absorbing potential (cgrid3d *; output; overwritten).
 * density     = current density (rgrid3d *; input).
 * rho0        = Desired density at the boundary (double, input).
 * region      = function that will multiply the absorbing potential (double) (*); input).
 *               Return: 0.0 = no absorption, \approx 0.1 = full absorption (values this
 *               large often makes the time propagation unstable).
 *               Even a linear function starting from the absorbing region edge (value 0) 
 *               and then gradually increasing up to 0.1 will often work well.
 * workspace   = temporary space needed for the operation (rgrid3d *).
 * c           = Multiply the absorption term by this constant (double complex) (usually 1.0).
 *
 * The calculated potential should be added to the potential that will be
 * propagated. Note that it is important that this is also included in the
 * possible predict-correct cycle as that improves the numericla stability.
 *
 */

EXPORT void grid3d_wf_absorb(cgrid3d *potential, rgrid3d *density, double rho0, double (*region)(void *, double, double, double), rgrid3d *workspace, double complex c) {

  rgrid3d_copy(workspace, density);
  rgrid3d_add(workspace, -rho0);
  rgrid3d_product_func(workspace, region, NULL);
  rgrid3d_multiply(workspace, -c);  // was -1.0
  grid3d_add_real_to_complex_im(potential, workspace);
}

/*
 * Damp the excitations of a wavefunction.
 *
 * Similar to the absorbing boundary potential (grid3d_wf_absorb), this function damps 
 * excitations entering the region so that in that region the wavefunction is 
 * psi = rho0 + 0 * i. 
 * This will damp the density AND the phase, so waves will reduce their amplitude and slow down.
 * The function is based on a Newton-Raphson minimization of the cost function
 * 
 *  C[Psi] = 0.5 * region(r) * |Psi-Psi_0|^2 
 *
 * and performs the following operation:
 *
 * 	Psi(r) -> Psi(r) - damping*region(r)*( Psi - Psi_0)
 *
 * (using imaginary potential is equivalent to this procedure but with a cost function
 * depending only on |psi|^2 and not on psi).
 *
 * wf		= wavefunction to be damped (wf3d)
 * rho0		= desired density value of the wavefunction, psi_0 = sqrt(rho0). (double)
 * damping	= strenght of the damping. (double) 
 * 		Too small and the excitations will reach the end of the box and bounce back.
 * 		Too large and the damping will effectively be just a hard wall and
 * 		the excitations will also bounce back.
 * 		Typical value used in testing is 0.01 * time_step.
 * region	= Region where the damping takes place (double complex function, but really real numbers!).
 * 		Expected to return values from 0 to 1, 0 meaning no damping at all and 1 meaning
 * 		max damping.
 * 		To go from 0 to 1 one can use a linear function or 1+tanh(x) or any other smooth function.
 * cworkspace	= temporary space needed for the operation (cgrid3d *).
 * farg		=  extra arguments needed for region (void *) . Can be NULL.
 *
 */
EXPORT void grid3d_damp_wf(wf3d *wf , double rho0 , double damping, double complex (*region)(void *, double, double, double) , cgrid3d *cworkspace, void *farg ){
	        double complex psi0 = sqrt(rho0) ;
		cgrid3d_copy(cworkspace, wf->grid ) ;
		cgrid3d_add_and_multiply(cworkspace, -psi0 , damping ) ;
		cgrid3d_product_func(cworkspace , region , farg ) ;
		cgrid3d_difference(wf->grid, wf->grid , cworkspace ) ;
}


/*
 * Calculate the probability flux x component.
 * 
 * gwf       = wavefunction for the operation (wf3d *).
 * flux_x    = x output grid containing the flux (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_probability_flux_x(const wf3d *gwf, rgrid3d *flux_x) {

  cgrid3d *grid = gwf->grid;
  long i, j, k, ij, ijnz, ny = grid->ny, nz = grid->nz, nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step), mass = gwf->mass;
  double *lvalue = flux_x->value;
  
#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta,grid,mass) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = cimag(conj(cgrid3d_value_at_index(grid, i, j, k)) * inv_delta * (cgrid3d_value_at_index(grid, i+1, j, k) - cgrid3d_value_at_index(grid, i-1, j, k))) * (HBAR / mass);
  }
}

/*
 * Calculate the probability flux y component.
 * 
 * gwf       = wavefunction for the operation (wf3d *).
 * flux_y    = y output grid containing the flux (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_probability_flux_y(const wf3d *gwf, rgrid3d *flux_y) {

  cgrid3d *grid = gwf->grid;
  long i, j, k, ij, ijnz, ny = grid->ny, nz = grid->nz, nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step), mass = gwf->mass;
  double *lvalue = flux_y->value;
  
#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta,grid,mass) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = cimag(conj(cgrid3d_value_at_index(grid, i, j, k)) * inv_delta * (cgrid3d_value_at_index(grid, i, j+1, k) - cgrid3d_value_at_index(grid, i, j-1, k))) * (HBAR / mass);
  }
}

/*
 * Calculate the probability flux z component.
 * 
 * gwf       = wavefunction for the operation (wf3d *).
 * flux_z    = z output grid containing the flux (rgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_probability_flux_z(const wf3d *gwf, rgrid3d *flux_z) {

  cgrid3d *grid = gwf->grid;
  long i, j, k, ij, ijnz, ny = grid->ny, nz = grid->nz, nxy = grid->nx * grid->ny;
  double inv_delta = 1.0 / (2.0 * grid->step), mass = gwf->mass;
  double *lvalue = flux_z->value;
  
#pragma omp parallel for firstprivate(ny,nz,nxy,lvalue,inv_delta,grid,mass) private(ij,ijnz,i,j,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    i = ij / ny;
    j = ij % ny;
    for(k = 0; k < nz; k++)
      lvalue[ijnz + k] = cimag(conj(cgrid3d_value_at_index(grid, i, j, k)) * inv_delta * (cgrid3d_value_at_index(grid, i, j, k+1) - cgrid3d_value_at_index(grid, i, j, k-1))) * (HBAR / mass);
  }
}

/*
 * Calculate the probability flux.
 *
 * gwf        = wavefunction for the operation (wf3d *).
 * flux_x     = x output grid containing the flux (rgrid3d *).
 * flux_y     = y output grid containing the flux (rgrid3d *).
 * flux_z     = z output grid containing the flux (rgrid3d *).
 *
 * No return value.
 *
 * NOTE: This is not the liquid velocity! Divide by density (rho) to get v (velocity):
 *       v_i = flux_i / rho (i = x, y, z).
 *
 */

EXPORT void grid3d_wf_probability_flux(const wf3d *gwf, rgrid3d *flux_x, rgrid3d *flux_y, rgrid3d *flux_z) {
  
  /*
   * J(r) = -i (hbar/2m) ( psi^* grad psi - psi grad psi^* )
   *      = (hbar/m) Im[ psi^* grad psi ] 
   */
  grid3d_wf_probability_flux_x(gwf, flux_x);
  grid3d_wf_probability_flux_y(gwf, flux_y);
  grid3d_wf_probability_flux_z(gwf, flux_z);
}

/*
 * Calculate the momentum x component.
 * 
 * gwf        = wavefunction for the operation (wf3d *).
 * momentum_x = output grid containing the momentum (cgrid3d *).
 * workspace  = temporary storage needed for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_momentum_x(const wf3d *gwf, cgrid3d *momentum_x, cgrid3d *workspace) {

  cgrid3d_copy(workspace, gwf->grid);
  cgrid3d_fft(workspace);
  cgrid3d_fft_gradient_x(workspace, momentum_x);
  cgrid3d_inverse_fft(momentum_x);
  cgrid3d_multiply(momentum_x, -I*HBAR / (2.0 * gwf->mass));
}

/*
 * Calculate the momentum y component.
 * 
 * gwf        = wavefunction for the operation (wf3d *).
 * momentum_y = output grid containing the momentum (cgrid3d *).
 * workspace  = temporary storage needed for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_momentum_y(const wf3d *gwf, cgrid3d *momentum_y, cgrid3d *workspace) {

  cgrid3d_copy(workspace, gwf->grid);
  cgrid3d_fft(workspace);
  cgrid3d_fft_gradient_y(workspace, momentum_y);
  cgrid3d_inverse_fft(momentum_y);
  cgrid3d_multiply(momentum_y, -I*HBAR / (2.0 * gwf->mass));
}

/*
 * Calculate the probability momentum z component.
 * 
 * gwf        = wavefunction for the operation (wf3d *).
 * momentum_z = output grid containing the momentum (cgrid3d *).
 * workspace  = temporary storage needed for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_momentum_z(const wf3d *gwf, cgrid3d *momentum_z, cgrid3d *workspace) {

  cgrid3d_copy(workspace, gwf->grid);
  cgrid3d_fft(workspace);
  cgrid3d_fft_gradient_z(workspace, momentum_z);
  cgrid3d_inverse_fft(momentum_z);
  cgrid3d_multiply(momentum_z, -I*HBAR / (2.0 * gwf->mass));
}

/*
 * Calculate the momentum.
 *
 * gwf        = wavefunction for the operation (wf3d *).
 * momentum_x = x output grid containing the momentum (cgrid3d *).
 * momentum_y = y output grid containing the momentum (cgrid3d *).
 * momentum_z = z output grid containing the momentum (cgrid3d *).
 * workspace  = additional storage needed for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_momentum(const wf3d *gwf, cgrid3d *momentum_x, cgrid3d *momentum_y, cgrid3d *momentum_z, cgrid3d *workspace) {
  
  grid3d_wf_momentum_x(gwf, momentum_x, workspace);
  grid3d_wf_momentum_y(gwf, momentum_y, workspace);
  grid3d_wf_momentum_z(gwf, momentum_z, workspace);
}

/*
 * Calculate energy for the wavefunction. 
 * Includes -Ekin *n if the frame of reference has momentum != 0.
 *
 * gwf       = wavefunction for the energy calculation (wf3d *).
 * potential = grid containing the potential (cgrid3d *).
 * workspace = additional storage needed (cgrid3d *).
 *
 * Returns the energy (double).
 *
 */

EXPORT double grid3d_wf_energy(const wf3d *gwf, const cgrid3d *potential, cgrid3d *workspace) {

  double mass=gwf->mass , kx = gwf->grid->kx0 , ky = gwf->grid->ky0 , kz = gwf->grid->kz0 ;
  double ekin = - HBAR * HBAR * (kx * kx + ky * ky + kz * kz) / (2. * mass ) ;

  if(ekin!=0.0)
	  ekin *= creal(cgrid3d_integral_of_square(gwf->grid)) ; 
  if (gwf->boundary == WF3D_DIRICHLET_BOUNDARY )
    return grid3d_wf_energy_cn(gwf, gwf, potential, workspace) + ekin;
  else if (gwf->boundary == WF3D_PERIODIC_BOUNDARY
		  || gwf->boundary == WF3D_NEUMANN_BOUNDARY
		  || gwf->boundary == WF3D_VORTEX_X_BOUNDARY
		  || gwf->boundary == WF3D_VORTEX_Y_BOUNDARY
		  || gwf->boundary == WF3D_VORTEX_Z_BOUNDARY )
      	  return grid3d_wf_energy_fft(gwf, potential, workspace) + ekin;
  else
    abort();
}

/*
 * Auxiliary routine for calculating the energy (Crank-Nicolson).
 * Users should rather call grid3d_wf_energy().
 *
 * gwfa      = (left) wavefunction for the energy calculation (wf3d *).
 * gwfb      = (right) wavefunction for the energy calculation (wf3d *).
 *             Normally gwfa = gwfb.
 * potential = Potential grid (cgrid3d *).
 * workspace = Additional workspace needed for the operation (cgrid3d *).
 * 
 * Returns the energy.
 *
 */

EXPORT double grid3d_wf_energy_cn(const wf3d *gwfa, const wf3d *gwfb, const cgrid3d *potential, cgrid3d *workspace) {
  
  long i;

  /* (-2m/hbar^2) T psi */
  cgrid3d_fd_laplace(gwfb->grid, workspace);
  cgrid3d_multiply(workspace, -HBAR * HBAR / (2.0 * gwfb->mass));

  /* V psi */
  if(potential)
    for(i = 0; i < gwfb->grid->nx * gwfb->grid->ny * gwfb->grid->nz; i++)
      workspace->value[i] += potential->value[i] * gwfb->grid->value[i];
  
  /* int psi^* (T + V) psi d^3r */
  return creal(cgrid3d_integral_of_conjugate_product(gwfa->grid, workspace));
}

/*
 * Auxiliary routine for calculating the energy (FFT).
 * Users should rather call grid3d_wf_energy().
 *
 * gwf       = wavefunction for the energy calculation (wf3d *).
 * potential = Potential grid (cgrid3d *).
 * workspace = Additional workspace needed for the operation (cgrid3d *).
 * 
 * Returns the energy.
 *
 */

EXPORT double grid3d_wf_energy_fft(const wf3d *gwf, const cgrid3d *potential, cgrid3d *workspace) {

  double en;

  /* delta (- k^2) fft[f(x)] / N */
  cgrid3d_copy(workspace, gwf->grid);
  cgrid3d_fft(workspace);
  
  en = -HBAR * HBAR / (2.0 * gwf->mass) * cgrid3d_fft_laplace_expectation_value(workspace, workspace);
  if(potential)
    en += creal(cgrid3d_grid_expectation_value(gwf->grid, potential));
  return en;
}

/*
 * Auxiliary routine for calculating kinetic energy (FFT).
 * This is used by grid3d_wf_energy_fft().
 * 
 * gwf       = wavefunction for the kinetic energy calculation (wf3d *).
 * workspace = additional workspace required for the operation (cgrid3d *).
 *
 * Returns the kinetic energy.
 *
 */

EXPORT double grid3d_wf_kinetic_energy_fft(const wf3d *gwf, cgrid3d *workspace) {

  /* delta (- k^2) fft[f(x)] / N */
  cgrid3d_copy(workspace, gwf->grid);
  cgrid3d_fft(workspace);
  
  return -HBAR*HBAR / (2.0 * gwf->mass) * cgrid3d_fft_laplace_expectation_value(workspace, workspace);
}

/*
 * Auxiliary routine for calculating potential energy.
 * 
 * gwf       = wavefunction for the kinetic energy calculation (wf3d *).
 * workspace = additional workspace required for the operation (cgrid3d *).
 *
 * Returns the potential energy.
 *
 */

EXPORT double grid3d_wf_potential_energy(const wf3d *gwf, const cgrid3d *potential) {

  return creal(cgrid3d_grid_expectation_value(gwf->grid, potential));
}

/*
 * Calculate energy (E) and the error (dE).
 *
 * gwf       = wavefunction for the operation (wf3d *).
 * potential = grid containing the potential (cgrid3d *).
 * workspace = additional workspace required for the operation (cgrid3d *).
 * error     = error estimate for energy (double *).
 * 
 * Returns the energy.
 *
 */

EXPORT double grid3d_wf_energy_and_error(const wf3d *gwf, const cgrid3d *potential, cgrid3d *workspace, double *error) { 

  double energy;
  
/* 
 * energy and its error
 * dE^2 = int dE^2 |psi|^2 dtau = ... = int E_local^2 |psi|^2 dtau - E_avg^2
 * dE = E_local - E_avg
 * E_local psi = H psi
 */

  /* T psi */
  if (gwf->boundary == WF3D_DIRICHLET_BOUNDARY || gwf->boundary == WF3D_NEUMANN_BOUNDARY) {
    cgrid3d_fd_laplace(gwf->grid, workspace);
    cgrid3d_multiply(workspace, -HBAR*HBAR / (2.0 * gwf->mass));
  }
  else if (gwf->boundary == WF3D_PERIODIC_BOUNDARY) {
    cgrid3d_copy(workspace, gwf->grid);
    cgrid3d_fft(workspace);
    cgrid3d_fft_laplace(workspace, workspace);
    cgrid3d_scaled_inverse_fft(workspace, -HBAR*HBAR / (2.0 * gwf->mass));
  }
  else {
    fprintf(stderr, "libgrid: Error in grid3d_wf_energy_and_error(). Invalid boundary condition, index = %d\n", gwf->boundary);
    abort();
  }
  
  /* H psi */
  cgrid3d_add_scaled_product(workspace, 1.0, potential, gwf->grid);
  
  /* int E_local^2 |psi|^2 dtau */
  *error = cgrid3d_integral_of_square(workspace);
  
  /* int E_local |psi|^2 dtau */
  cgrid3d_conjugate_product(workspace, gwf->grid, workspace);
  energy = creal(cgrid3d_integral(workspace));
  
  /* sqrt( int E_local^2 |psi|^2 dtau - ( int E_local |psi|^2 dtau )^2 ) */
  *error = sqrt(*error - energy * energy);
  
  return energy;
}

/*
 * Propagate wavefunction in time subject to given potential.
 *
 * gwf         = wavefunction to be propagated (wf3d *).
 * potential   = grid containing the potential (cgrid3d *).
 * sq_grad_pot = grid containing square of potential gradient (cgrid3d *).
 * time        = time step (double complex). Note this may be either real or imaginary.
 * workspace   = additional workspace needed for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate(wf3d *gwf, const cgrid3d *potential, const cgrid3d *sq_grad_pot, double complex time, cgrid3d *workspace) {  
  
  double complex half_time = 0.5 * time;
  double complex one_sixth_time = time / 6.0;
  double complex two_thirds_time = 2.0 * time / 3.0;
  
  if (gwf->boundary == WF3D_DIRICHLET_BOUNDARY && gwf->propagator == WF3D_2ND_ORDER_PROPAGATOR) {    
    grid3d_wf_propagate_potential(gwf, potential, half_time);
    grid3d_wf_propagate_kinetic_cn(gwf, time, workspace);
    grid3d_wf_propagate_potential(gwf, potential, half_time);
  } else if ((gwf->boundary == WF3D_PERIODIC_BOUNDARY || gwf->boundary == WF3D_NEUMANN_BOUNDARY || gwf->boundary == WF3D_VORTEX_X_BOUNDARY || gwf->boundary == WF3D_VORTEX_Y_BOUNDARY || gwf->boundary == WF3D_VORTEX_Z_BOUNDARY)
	     && gwf->propagator == WF3D_2ND_ORDER_PROPAGATOR) {
    grid3d_wf_propagate_potential(gwf, potential, half_time);
    grid3d_wf_propagate_kinetic_fft(gwf, time);
    grid3d_wf_propagate_potential(gwf, potential, half_time);
  } else if (gwf->boundary == WF3D_PERIODIC_BOUNDARY 
            && gwf->propagator == WF3D_4TH_ORDER_PROPAGATOR) {    
    grid3d_wf_propagate_potential(gwf, potential, one_sixth_time);
    grid3d_wf_propagate_kinetic_fft(gwf, half_time);    
    cgrid3d_copy(workspace, potential);
    cgrid3d_add_scaled(workspace, (1/48.0 * HBAR * HBAR / gwf->mass) * sqnorm(time), sq_grad_pot);	
    grid3d_wf_propagate_potential(gwf, workspace, two_thirds_time);
    
    grid3d_wf_propagate_kinetic_fft(gwf, half_time);
    grid3d_wf_propagate_potential(gwf, potential, one_sixth_time);
  } else {
    fprintf(stderr, "libgrid: Error in grid3d_wf_propagate(). Unknown propagator - boundary value combination (propagator index = %d, boundary index = %d).\n", gwf->propagator, gwf->boundary);
    abort();
  }
}

/*
 * Auxiliary routine to propagate kinetic energy using FFT.
 *
 * gwf  = wavefunction to be propagated (wf3d *).
 * time = time step (double complex).
 *
 * TODO: Dirichlet BC could be implemented using sin transformation and Neumann by cos transformation. This
 * would eliminate the current CN based routines completely. The FFT methods have the best accuracy.
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_fft(wf3d *gwf, double complex time) {

  long i,j,k, ij, ijnz, nx,ny,nz, nxy;
  double kx, ky, kz, lz, step, norm, mass;
  double kx0 = gwf->grid->kx0 , ky0 = gwf->grid->ky0 , kz0 = gwf->grid->kz0 ;
  double complex *value = gwf->grid->value;
  
  cgrid3d_fft(gwf->grid);
  
  nx = gwf->grid->nx;
  ny = gwf->grid->ny;
  nz = gwf->grid->nz;
  nxy = nx * ny;
  step = gwf->grid->step;
  mass = gwf->mass;
  
  /* f(x) = ifft[fft[f(x)]] / N */
  norm = gwf->grid->fft_norm;
  if( gwf->boundary == WF3D_NEUMANN_BOUNDARY  ||
      gwf->boundary == WF3D_VORTEX_X_BOUNDARY ||
      gwf->boundary == WF3D_VORTEX_Y_BOUNDARY ||
      gwf->boundary == WF3D_VORTEX_Z_BOUNDARY)
    {
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,value,time,mass,kx0,ky0,kz0) private(i,j,ij,ijnz,k,kx,ky,kz,lz) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      i = ij / ny;
      j = ij % ny;
      ijnz = ij * nz;
      
      kx = M_PI * i / (nx * step) - kx0;
      ky = M_PI * j / (ny * step) - ky0;
      lz = nz * step;
      
      for(k = 0; k < nz; k++) {
        kz = M_PI * k / lz - kz0;
        
        /* psi(t+dt) = psi(t) exp( - i (hbar^2 * k^2 / 2m) dt / hbar ) */	  
        value[ijnz + k] *= norm * cexp(-I * time * HBAR * HBAR * (kx*kx + ky*ky + kz*kz) / (2.0 * mass * HBAR));
      }
    }
  }
  else{
#pragma omp parallel for firstprivate(norm,nx,ny,nz,nxy,step,value,time,mass,kx0,ky0,kz0) private(i,j,ij,ijnz,k,kx,ky,kz,lz) default(none) schedule(runtime)
    for(ij = 0; ij < nxy; ij++) {
      i = ij / ny;
      j = ij % ny;
      ijnz = ij * nz;
      
      /* 
       * k = 2 pi n / L 
       * if k < n/2, k = k
       * else k = -k
       */
      if (i < nx / 2)
        kx = 2.0 * M_PI * i / (nx * step) - kx0;
      else 
        kx = 2.0 * M_PI * (i - nx) / (nx * step) -kx0;
      
      if (j < ny / 2)
        ky = 2.0 * M_PI * j / (ny * step) - ky0;
      else 
        ky = 2.0 * M_PI * (j - ny) / (ny * step) - ky0;
      lz = nz * step;
      
      for(k = 0; k < nz; k++) {
        if (k < nz / 2)
          kz = 2.0 * M_PI * k / lz - kz0; 
        else 
          kz = 2.0 * M_PI * (k - nz) / lz - kz0;
        
        /* psi(t+dt) = psi(t) exp( - i (hbar^2 * k^2 / 2m) dt / hbar ) */	  
        value[ijnz + k] *= norm * cexp(-I * time * HBAR * HBAR * (kx*kx + ky*ky + kz*kz) / (2 * mass * HBAR));
      }
    } 
  }
  
  cgrid3d_inverse_fft(gwf->grid);
}

/*
 * Auxiliary routine to propagate kinetic energy (Crank-Nicolson).
 * Users should call grid3d_wf_propagate().
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional workspace required for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_cn(wf3d *gwf, double complex time, cgrid3d *workspace) {

  if (gwf->boundary == WF3D_DIRICHLET_BOUNDARY)
    grid3d_wf_propagate_kinetic_cn_dbc(gwf, time, workspace);
  else if (gwf->boundary == WF3D_NEUMANN_BOUNDARY)
    grid3d_wf_propagate_kinetic_cn_nbc(gwf, time, workspace);
  else if (gwf->boundary == WF3D_PERIODIC_BOUNDARY)
    grid3d_wf_propagate_kinetic_cn_pbc(gwf, time, workspace);
  else {
    fprintf(stderr, "libgrid: Error in grid3d_wf_propagate_cn(). Unknown boundary condition (index = %d).\n", gwf->boundary);
    abort();
  }
}

/*
 * Auxiliary routine to propagate potential energy.
 * Users should call grid3d_wf_propaget().
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * potential = grid containing the potential (cgrid3d *).
 * time      = time step (double complex).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_potential(wf3d *gwf, const cgrid3d *potential, double complex time) {

  long ij, ijnz, k, nxy = gwf->grid->nx * gwf->grid->ny, nz = gwf->grid->nz;
  double complex c, *psi = gwf->grid->value, *pot = potential->value;
  
  /* c = - i dt / hbar */
  c = -I * time / HBAR;
  
#pragma omp parallel for firstprivate(nz,nxy,psi,pot,c) private(k, ij, ijnz) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      /* psi(t+dt) = exp(- i V dt / hbar) psi(t) */
      psi[ijnz + k] *= cexp(c * pot[ijnz + k]);
  }
}

/*
 * Project "gwfb" out from "gwfa".
 *
 * gwfa = input wavefunction (wf3d *).
 * gwfb = this will be projected out from gwfa (wf3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_project_out(wf3d *gwfa, const wf3d *gwfb) {

  double complex overlap = grid3d_wf_overlap(gwfa, gwfb);

  cgrid3d_add_scaled(gwfa->grid, -overlap, gwfb->grid);
}

/*
 * "traditional" diagonalization of Hamiltonian.
 *
 * gwf    = an array of wavefunctions (wf3d **).
 * states = number of states (int).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_diagonalize(wf3d **gwf, int states) {

  int i,j;
  double *eigenvalue = (double *) malloc(states * sizeof(double));
  double complex *overlap = (double complex *) malloc(states * states * sizeof(double complex));
  wf3d *gwf_tmp;
  
  if (states == 1) {
    grid3d_wf_normalize(gwf[0]);
    return;
  }
  
  /* overlap matrix */
  for(i = 0; i < states; i++) {
    for(j = 0; j <= i; j++) {
      /* fortran (column major) matrix order, i is row (minor) index, j is column (major) index */
      overlap[i + j * states] = grid3d_wf_overlap(gwf[i], gwf[j]);
      overlap[j + i * states] = conj(overlap[i + j * states]);
    }
  }
  
  /* diagonalize */
  grid_hermitian_eigenvalue_problem(eigenvalue, overlap, states);
  
  /* phi_i = 1 / sqrt(m_i) C_ij psi_j, C (row major) matrix order ???is it??? */
  for(i = 0; i < states; i++)
    for(j = 0; j < states; j++)
      overlap[i * states + j] /= sqrt(eigenvalue[i]);
  
  grid3d_wf_linear_transform(gwf, overlap, states);
  
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
 * gwf       = an array of wavefunctions (wf3d **).
 * transform = transformation matrix (double complex *).
 * states    = number of states (int).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_linear_transform(wf3d **gwf, double complex *transform, int states) {

  int p,q, offset;
  long ij,k, nxy, nz, ijnz;
  double complex **value, *tmp;
  
  nxy = gwf[0]->grid->nx * gwf[0]->grid->ny;
  nz = gwf[0]->grid->nz;
  
  /* + 16 to prevent "write locks" */
  tmp = (double complex *) malloc(omp_get_max_threads() * (states + 16) * sizeof(double complex));
  
  value = (double complex **) malloc(states * sizeof(double complex *));

  for(p = 0; p < states; p++)
    value[p] = gwf[p]->grid->value;
  
#pragma omp parallel for firstprivate(nz,nxy,value,states,transform,tmp) private(ij,k,p,q,ijnz,offset) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij*nz;
    
    offset = (states + 16) * omp_get_thread_num();
    
    for(k = 0; k < nz; k++) {
      for(p = 0; p < states; p++)
        tmp[offset + p] = 0.0;
      
      for(p = 0; p < states; p++)
        for(q = 0; q < states; q++)
          tmp[offset + p] += transform[p * states + q] * value[q][ijnz + k];
      
      for(p = 0; p < states; p++)
        value[p][ijnz + k] = tmp[offset + p];
    }
  }
  
  free(value);
  free(tmp);
}

/*
 * Calcuate square of potential gradient.
 *
 * sq_grad_pot = output grid (cgrid3d *).
 * potential   = potental input grid (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_square_of_potential_gradient(cgrid3d *sq_grad_pot, const cgrid3d *potential, cgrid3d *workspace, cgrid3d *workspace2) {

  cgrid3d_copy(sq_grad_pot, potential);
  cgrid3d_fft(sq_grad_pot);
  cgrid3d_fft_gradient(sq_grad_pot, sq_grad_pot, workspace, workspace2);
  
  cgrid3d_inverse_fft(sq_grad_pot);
  cgrid3d_inverse_fft(workspace);
  cgrid3d_inverse_fft(workspace2);
  
  cgrid3d_conjugate_product(sq_grad_pot, sq_grad_pot, sq_grad_pot);
  cgrid3d_conjugate_product(workspace, workspace, workspace);
  cgrid3d_conjugate_product(workspace2, workspace2, workspace2);
  
  cgrid3d_sum(workspace, workspace, workspace2);
  cgrid3d_sum(sq_grad_pot, sq_grad_pot, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy using Dirichlet boundaries (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

/* TODO: REMOVE LAPACK LATER */ 
/* #define USE_LAPACK */

EXPORT void grid3d_wf_propagate_kinetic_cn_dbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  /*
   * exp( -i (Tx + Ty + Tz) dt / hbar ) 
   *   = exp( -i Tx dt / hbar ) exp( -i Ty dt / hbar ) exp( -i Tz dt / hbar ) + O(dt^2)
   *   
   */

  grid3d_wf_propagate_kinetic_x_cn_dbc(gwf, time, workspace);
  grid3d_wf_propagate_kinetic_y_cn_dbc(gwf, time, workspace);
  grid3d_wf_propagate_kinetic_z_cn_dbc(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along x using Dirichlet boundaries (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_x_cn_dbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, boundary, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind, nynz;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, nz = gwf->grid->nz;
#ifdef USE_LAPACK
  /* NOTE: parallel exec does not work yet */
  char TRANS = 'N';
  int N = nx, NRHS = 1, ldb = nx, info = 0;
  double complex *dl = &wrk[0];
  double complex *d = &wrk[nx];
  double complex *du = &wrk[2*nx];
  double complex *du2 = &wrk[3*nx];
  int *ipiv = &wrk[4*nx];
  double complex *b = &wrk[5*nx];
#else
  double complex *b, *ld, *ad;
#endif

  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* create cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;
  
#ifdef USE_LAPACK
  for(i = 0; i < nx; i++) {
    dl[i] = 1.0;
    d[i] = ca;
    du[i] = 1.0;
    ipiv[i] = i;
  }
  zgttrf_(&N, dl, d, du, du2, ipiv, &info);
  if(info != 0) {
    fprintf(stderr, "zgttrf info = %d\n", info);
    exit(1);
  }
#else
  ad = &wrk[omp_get_max_threads() * nx];
  /* create matrix */
  for(i = 0; i < nx; i++)
    ad[i] = ca;  
  
  /* create decomposition */
  grid_cholesky_decomposition(ad, nx, 0);
  /* matrix is replaced by decomposition */
  ld = ad;
#endif

  /* solve problem for each y and z*/
#pragma omp parallel for firstprivate(c,psi,nx,ny,nz,ld,wrk) private(i,j,k,ind,nynz,tid,b,boundary) default(none) schedule(runtime)
  for(j = 0; j < ny; j++) {
    for(k = 0; k < nz; k++) {
      nynz = ny*nz;
#ifndef USE_LAPACK
      tid = omp_get_thread_num();

      /* use diffent workspace for different threads  */
      b = &wrk[tid * nx];
#endif
      
      /* create right-hand side vector */
      for(i = 1; i < nx-1; i++) {
        ind = (i*ny + j)*nz + k;
        /* (C - dx^2 laplace) psi(t) */
	b[i] = c * psi[ind] - (psi[ind + nynz] - 2.0 * psi[ind] + psi[ind - nynz]);
      }
      
      /* dirichlet boundaries */
      boundary = 0.0;
      i = 0; 
      ind = (i*ny + j)*nz + k;
      b[i] = c * psi[ind] - (psi[ind + nynz] - 2.0 * psi[ind] + boundary);
      i = nx-1;
      ind = (i*ny + j)*nz + k;
      b[i] = c * psi[ind] - (boundary - 2.0 * psi[ind] + psi[ind - nynz]);
      
#ifdef USE_LAPACK
      zgttrs_(&TRANS, &N, &NRHS, dl, d, du, du2, ipiv, b, &ldb, &info);
      if(info != 0) {
	fprintf(stderr, "ZGTTRS info = %d\n", info);
	exit(1);
      }
#else
      /* solve */
      grid_cholesky_substitute(ld, b, nx, 0);
#endif
      /* copy */
      for(i = 0; i < nx; i++)
        psi[(i*ny + j)*nz + k] = b[i];
    }
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along y using Dirichlet boundaries (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_y_cn_dbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, boundary, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind, iny;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, nz = gwf->grid->nz;
#ifdef USE_LAPACK
  /* NOTE: parallel exec does not work yet */
  char TRANS = 'N';
  int N = ny, NRHS = 1, ldb = ny, info = 0;
  double complex *dl = &wrk[0];
  double complex *d = &wrk[ny];
  double complex *du = &wrk[2*ny];
  double complex *du2 = &wrk[3*ny];
  int *ipiv = &wrk[4*ny];
  double complex *b = &wrk[5*ny];
#else
  double complex *b, *ld, *ad;
#endif
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);
  
  /* create cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

#ifdef USE_LAPACK
  for(i = 0; i < ny; i++) {
    dl[i] = 1.0;
    d[i] = ca;
    du[i] = 1.0;
    ipiv[i] = i;
  }
  zgttrf_(&N, dl, d, du, du2, ipiv, &info);
  if(info != 0) {
    fprintf(stderr, "zgttrf info = %d\n", info);
    exit(1);
  }
#else
  ad = &wrk[omp_get_max_threads() * ny];
  
  /* create matrix */
  for(j = 0; j < ny; j++)
    ad[j] = ca;  
  
  /* create decomposition */
  grid_cholesky_decomposition(ad, ny, 0);
  /* matrix is replaced by decomposition */
  ld = ad;
#endif
  
  /* solve problem for each x and z*/
#pragma omp parallel for firstprivate(c,psi,nx,ny,nz,ld,wrk) private(i,j,k,ind,iny,tid,b,boundary) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {
    iny = i * ny;
    for(k = 0; k < nz; k++) {
#ifndef USE_LAPACK
      tid = omp_get_thread_num();
      
      /* use diffent workspace for different threads  */
      b = &wrk[tid * ny];
#endif
      
      /* create right-hand side vector */
      for(j = 1; j < ny-1; j++) {
        ind = (iny + j)*nz + k;
        /* (C - dx^2 laplace) psi(t) */
	b[j] = c * psi[ind] - (psi[ind + nz] - 2.0 * psi[ind] + psi[ind - nz]);
      }
      
      /* dirichlet boundaries */
      boundary = 0.0;
      j = 0;
      ind = (iny + j)*nz + k;
      b[j] = c * psi[ind] - (psi[ind + nz] - 2.0 * psi[ind] + boundary);
      j = ny-1; 
      ind = (iny + j)*nz + k;
      b[j] = c * psi[ind] - (boundary - 2.0 * psi[ind] + psi[ind - nz]);
      
#ifdef USE_LAPACK
      zgttrs_(&TRANS, &N, &NRHS, dl, d, du, du2, ipiv, b, &ldb, &info);
      if(info != 0) {
	fprintf(stderr, "ZGTTRS info = %d\n", info);
	exit(1);
      }
#else
      /* solve */
      grid_cholesky_substitute(ld, b, ny, 0);
#endif      
      /* copy */
      for(j = 0; j < ny; j++)
        psi[(iny + j)*nz + k] = b[j];
    }
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along z using Dirichlet boundaries (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_z_cn_dbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, boundary, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind, inyjnz;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, nz = gwf->grid->nz;
#ifdef USE_LAPACK
  /* NOTE: parallel exec does not work yet */
  char TRANS = 'N';
  int N = nz, NRHS = 1, ldb = nz, info = 0;
  double complex *dl = &wrk[0];
  double complex *d = &wrk[nz];
  double complex *du = &wrk[2*nz];
  double complex *du2 = &wrk[3*nz];
  double complex *ipiv = &wrk[4*nz];
  double complex *b = &wrk[5*nz];
#else
  double complex *b, *ld, *ad;
#endif
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);
  
  /* create cholesky decompostions for kinetic x */
  
  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;
#ifdef USE_LAPACK
  for(i = 0; i < nz; i++) {
    dl[i] = 1.0;
    d[i] = ca;
    du[i] = 1.0;
    ipiv[i] = i;
  }
  zgttrf_(&N, dl, d, du, du2, ipiv, &info);
  if(info != 0) {
    fprintf(stderr, "zgttrf info = %d\n", info);
    exit(1);
  }
#else
  ad = &wrk[omp_get_max_threads() * nz];
  
  /* create matrix */
  for(k = 0; k < nz; k++)
    ad[k] = ca;  
  
  /* create decomposition */
  grid_cholesky_decomposition(ad, nz, 0);
  /* matrix is replaced by decomposition */
  ld = ad;
#endif

  /* solve problem for each x and y */
#pragma omp parallel for firstprivate(c,psi,nx,ny,nz,ld,wrk) private(i,j,k,ind,inyjnz,tid,b,boundary) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) {
    for(j = 0; j < ny; j++) {
      inyjnz = (i*ny+j)*nz;
#ifndef USE_LAPACK
      tid = omp_get_thread_num();
      
      /* use diffent workspace for different threads  */
      b = &wrk[tid * nz];
#endif
      
      /* create right-hand side vector */
      for(k = 1; k < nz-1; k++) {
        ind = inyjnz + k;
        /* (C - dx^2 laplace) psi(t) */
	b[k] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + psi[ind-1]);
      }
      
      /* dirichlet boundaries */
      boundary = 0.0;
      k = 0;
      ind = inyjnz + k;
      b[k] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + boundary);
      k = nz-1;
      ind = inyjnz + k;
      b[k] = c * psi[ind] - (boundary - 2.0 * psi[ind] + psi[ind-1]);
      
#ifdef USE_LAPACK
      zgttrs_(&TRANS, &N, &NRHS, dl, d, du, du2, ipiv, b, &ldb, &info);
      if(info != 0) {
	fprintf(stderr, "ZGTTRS info = %d\n", info);
	exit(1);
      }
#else
      /* solve */
      grid_cholesky_substitute(ld, b, nz, 0);
#endif
      
      /* copy */
      for(k = 0; k < nz; k++)
        psi[inyjnz + k] = b[k];
    }
  }
}

/*
 * Auxiliary routine for propagating kinetic energy using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_cn_nbc(wf3d *gwf, double complex time, cgrid3d *workspace) {
  
  /*
   * exp( -i (Tx + Ty + Tz) dt / hbar ) 
   *   = exp( -i Tx dt / hbar ) exp( -i Ty dt / hbar ) exp( -i Tz dt / hbar ) + O(dt^2)
   *   
   */
  
  grid3d_wf_propagate_kinetic_x_cn_nbc(gwf, time, workspace);
  grid3d_wf_propagate_kinetic_y_cn_nbc(gwf, time, workspace);
  grid3d_wf_propagate_kinetic_z_cn_nbc(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along x using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_x_cn_nbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, jk;
  long nz = gwf->grid->nz, nyz = ny * nz;
  double complex *du, *dl, *d, *b;
  
  if((4 * (grid_threads()-1) + 3) >= nyz) {
    fprintf(stderr, "libgrid: insufficient workspace in grid3d_wf_propagate_kinetic_x_cn_nbc().\n");
    exit(1);
  }

  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nx,nz,nyz,psi,wrk,c) private(tid,i,j,k,jk,dl,d,du,b,ind) default(none) schedule(runtime)
  for (jk = 0; jk < nyz; jk++) {  /* for each (y,z) */
    j = jk / nz;
    k = jk % nz;
    tid = omp_get_thread_num();
    dl = &wrk[nx * (4 * tid + 0)];
    d  = &wrk[nx * (4 * tid + 1)];
    du = &wrk[nx * (4 * tid + 2)];
    b  = &wrk[nx * (4 * tid + 3)];

    for(i = 0; i < nx; i++) {
      if(i == nx-1) dl[nx-1] = 2.0;             /* Neumann outer boundary */
      else if(i > 0) dl[i] = 1.0;
      if(i == 0) du[0] = 2.0;                   /* Neumann inner boundary */
      else if(i < nx-1) du[i] = 1.0;
      d[i] = ca;
    }

    /* create right-hand side vector */
    for(i = 1; i < nx - 1; i++) {
      ind = i * nyz + j * nz + k;
      /* (C - dx^2 laplace) psi(t) */
      b[i] = c * psi[ind] - (psi[ind + nyz] - 2.0 * psi[ind] + psi[ind - nyz]);
    }
 
    /* neumann boundaries */
    i = 0; ind = i * nyz + j * nz + k;
    b[i] = c * psi[ind] - (2.0 * psi[ind + nyz] - 2.0 * psi[ind]);
    i = nx - 1; ind = i * nyz + j * nz + k;
    b[i] = c * psi[ind] - (2.0 * psi[ind - nyz] - 2.0 * psi[ind]);
    
    grid_solve_tridiagonal_system(nx, dl, d, du, b, dl);

    for(i = 0; i < nx; i++)
      psi[i * nyz + j * nz + k] = dl[i];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along y using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_y_cn_nbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, ik, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, nz = gwf->grid->nz, nxz = nx * nz, nyz = ny * nz;
  double complex *du, *dl, *d, *b;
  
  if((4 * (grid_threads()-1) + 3) >= nyz) {
    fprintf(stderr, "libgrid: insufficient workspace in grid3d_wf_propagate_kinetic_y_cn_nbc().\n");
    exit(1);
  }

  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nx,ny,nz,nxz,nyz,psi,wrk,c) private(tid,i,j,k,ik,dl,d,du,b,ind) default(none) schedule(runtime)
  for (ik = 0; ik < nxz; ik++) {  /* for each (x,z) */
    i = ik / nz;
    k = ik % nz;
    tid = omp_get_thread_num();
    dl = &wrk[ny * (4 * tid + 0)];
    d  = &wrk[ny * (4 * tid + 1)];
    du = &wrk[ny * (4 * tid + 2)];
    b  = &wrk[ny * (4 * tid + 3)];

    for(j = 0; j < ny; j++) {
      if(j == ny-1) dl[ny-1] = 2.0;              /* Neumann outer boundary */
      else if(j > 0) dl[j] = 1.0;
      if(j == 0) du[0] = 2.0;                   /* Neumann inner boundary */
      else if(j < ny-1) du[j] = 1.0;
      d[j] = ca;
    }

    /* create right-hand side vector */
    for(j = 1; j < ny - 1; j++) {
      ind = i * nyz + j * nz + k;
      /* (C - dx^2 laplace) psi(t) */
      b[j] = c * psi[ind] - (psi[ind + nz] - 2.0 * psi[ind] + psi[ind - nz]);
    }
 
    /* neumann boundaries */
    j = 0; 
    ind = i * nyz + j * nz + k;
    b[j] = c * psi[ind] - (2.0 * psi[ind + nz] - 2.0 * psi[ind]);

    j = ny - 1;
    ind = i * nyz + j * nz + k;
    b[j] = c * psi[ind] - (2.0 * psi[ind - nz] - 2.0 * psi[ind]);
    
    grid_solve_tridiagonal_system(ny, dl, d, du, b, dl);

    for(j = 0; j < ny; j++)
      psi[i * nyz + j * nz + k] = dl[j];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along z using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_z_cn_nbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind, inz;
  long nx = gwf->grid->nx;
  long ny = gwf->grid->ny, nz = gwf->grid->nz;
  long ij, k;
  double complex *dl, *d, *du, *b;
  
  if((4 * (grid_threads()-1) + 3) >= nx*ny) {
    fprintf(stderr, "libgrid: insufficient workspace in grid3d_wf_propagate_kinetic_z_cn_nbc().\n");
    exit(1);
  }

  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = 4 i m dx^2 / (hbar dt)
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);
  
  /* (C + dr^2 laplace + (1/r) dr^2 grad) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nx,ny,nz,psi,wrk,c,step) private(tid,ij,k,dl,d,du,b,ind,inz) default(none) schedule(runtime)
  for(ij = 0; ij < nx*ny; ij++) { /* for each (x,y) */
    tid = omp_get_thread_num();
    dl = &wrk[nz * (4 * tid + 0)];
    d  = &wrk[nz * (4 * tid + 1)];
    du = &wrk[nz * (4 * tid + 2)];
    b  = &wrk[nz * (4 * tid + 3)];

    for(k = 0; k < nz; k++) { /* over z */
      if(k == nz-1) dl[nz-1] = 2.0;            /* Neumann outer boundary */
      else if(k > 0) dl[k] = 1.0; 
      if(k == 0) du[0] = 2.0;                  /* Neumann inner boundary */
      else if(k < nz-1) du[k] = 1.0;
      d[k] = ca;
    }

    /* create right-hand side vector */
    inz = ij * nz;
    for(k = 1; k < nz - 1; k++) {
      ind = inz + k;
      /* (C - dr^2 laplace - (1/r) dr^2 grad) psi(t) */
      b[k] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + psi[ind - 1]);
    }
    
    /* neumann boundaries */
    k = 0; 
    ind = inz + k;
    b[k] = c * psi[ind] - (2.0 * psi[ind + 1] - 2.0 * psi[ind]);
 
    k = ny - 1; 
    ind = inz + k;
    b[k] = c * psi[ind] - (2.0 * psi[ind - 1] - 2.0 * psi[ind]);

    grid_solve_tridiagonal_system(nz, dl, d, du, b, dl);

    /* copy */
    for(k = 0; k < nz; k++)
      psi[inz + k] = dl[k];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy using Neumann boundary condition (Crank-Nicolson).
 * This is a special version of the CN-NBC routine that includes the necessary angular momentum terms,
 * which are needed to induce liquid rotation. The rotation is assumed to take place about the
 * z-axis (hardwired).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * omega     = angular velocity for rotation (double).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_cn_nbc_rot(wf3d *gwf, double complex time, double omega, cgrid3d *workspace) {
  
  /*
   * exp( -i (Tx + Ty + Tz) dt / hbar ) 
   *   = exp( -i Tx dt / hbar ) exp( -i Ty dt / hbar ) exp( -i Tz dt / hbar ) + O(dt^2)
   *   
   */
  
  grid3d_wf_propagate_kinetic_x_cn_nbc_rot(gwf, time, omega, workspace);
  grid3d_wf_propagate_kinetic_y_cn_nbc_rot(gwf, time, omega, workspace);
  grid3d_wf_propagate_kinetic_z_cn_nbc(gwf, time, workspace);    /* no changes needed for z */
}

/*
 * Auxiliary routine for propagating kinetic energy along x using Neumann boundary condition (Crank-Nicolson).
 * Includes rotation - see above.
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * omega     = angular velocity for rotation (double).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_x_cn_nbc_rot(wf3d *gwf, double complex time, double omega, cgrid3d *workspace) {

  double complex c, ca, c2, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step, y, y0 = gwf->grid->y0;
  long tid, ind;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, jk;
  long nz = gwf->grid->nz, nyz = ny * nz;
  double complex *du, *dl, *d, *b;

  if((4 * (grid_threads()-1) + 3) >= nyz) {
    fprintf(stderr, "libgrid: insufficient workspace in grid3d_wf_propagate_kinetic_x_cn_nbc_rot().\n");
    exit(1);
  }
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

  c2 = gwf->mass * omega * I * step / HBAR;

#pragma omp parallel for firstprivate(ca,c2,nx,ny,nz,nyz,psi,wrk,c,y0,step) private(tid,y,i,j,k,jk,dl,d,du,b,ind) default(none) schedule(runtime)
  for (jk = 0; jk < nyz; jk++) {  /* for each (y,z) */
    j = jk / nz;
    k = jk % nz;
    y = (j - ny/2) * step - y0;
    tid = omp_get_thread_num();
    dl = &wrk[nx * (4 * tid + 0)];
    d  = &wrk[nx * (4 * tid + 1)];
    du = &wrk[nx * (4 * tid + 2)];
    b  = &wrk[nx * (4 * tid + 3)];

    for(i = 0; i < nx; i++) {
      if(i == nx-1) dl[nx-1] = 2.0;             /* Neumann outer boundary */
      else if(i > 0) dl[i] = 1.0 - c2 * y;
      if(i == 0) du[0] = 2.0;                   /* Neumann inner boundary */
      else if(i < nx-1) du[i] = 1.0 + c2 * y;
      d[i] = ca;
    }

    /* create right-hand side vector */
    for(i = 1; i < nx - 1; i++) {
      ind = i * nyz + j * nz + k;
      /* (C - dx^2 laplace - C2 * grad) psi(t) */
      b[i] = c * psi[ind] - (psi[ind + nyz] - 2.0 * psi[ind] + psi[ind - nyz]) - c2 * y * (psi[ind + nyz] - psi[ind - nyz]);
    }
 
    /* neumann boundaries */
    i = 0; 
    ind = i * nyz + j * nz + k;
    b[i] = c * psi[ind] - (2.0 * psi[ind + nyz] - 2.0 * psi[ind]);

    i = nx - 1;
    ind = i * nyz + j * nz + k;
    b[i] = c * psi[ind] - (2.0 * psi[ind - nyz] - 2.0 * psi[ind]);
    
    grid_solve_tridiagonal_system(nx, dl, d, du, b, dl);

    for(i = 0; i < nx; i++)
      psi[i * nyz + j * nz + k] = dl[i];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along y using Neumann boundary condition (Crank-Nicolson).
 * Includes rotation - see above.
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * omega     = angular velocity for rotation (double).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_y_cn_nbc_rot(wf3d *gwf, double complex time, double omega, cgrid3d *workspace) {

  double complex c, ca, c2, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step, x, x0 = gwf->grid->x0;
  long tid, ind;
  long i, ik, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, nz = gwf->grid->nz, nxz = nx * nz, nyz = ny * nz;
  double complex *du, *dl, *d, *b;
  
  if((4 * (grid_threads()-1) + 3) >= nxz) {
    fprintf(stderr, "libgrid: insufficient workspace in grid3d_wf_propagate_kinetic_y_cn_nbc_rot().\n");
    exit(1);
  }

  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

  c2 = -gwf->mass * omega * I * step / HBAR;

#pragma omp parallel for firstprivate(ca,c2,nx,ny,nz,nxz,nyz,psi,wrk,c,x0,step) private(x,tid,i,j,k,ik,dl,d,du,b,ind) default(none) schedule(runtime)
  for (ik = 0; ik < nxz; ik++) {  /* for each (x,z) */
    i = ik / nz;
    k = ik % nz;
    x = (i - nx/2) * step - x0;
    tid = omp_get_thread_num();
    dl = &wrk[ny * (4 * tid + 0)];
    d  = &wrk[ny * (4 * tid + 1)];
    du = &wrk[ny * (4 * tid + 2)];
    b  = &wrk[ny * (4 * tid + 3)];

    for(j = 0; j < ny; j++) {
      if(j == ny-1) dl[ny-1] = 2.0;             /* Neumann outer boundary */
      else if(j > 0) dl[j] = 1.0 - c2 * x;
      if(j == 0) du[0] = 2.0;                   /* Neumann inner boundary */
      else if(j < ny-1) du[j] = 1.0 + c2 * x;
      d[j] = ca;
    }

    /* create right-hand side vector */
    for(j = 1; j < ny - 1; j++) {
      ind = i * nyz + j * nz + k;
      /* (C - dx^2 laplace) psi(t) */
      b[j] = c * psi[ind] - (psi[ind + nz] - 2.0 * psi[ind] + psi[ind - nz]) - c2 * x * (psi[ind + nz] - psi[ind - nz]);
    }
 
    /* neumann boundaries */
    j = 0; 
    ind = i * nyz + j * nz + k;
    b[j] = c * psi[ind] - (2.0 * psi[ind + nz] - 2.0 * psi[ind]);

    j = ny - 1;
    ind = i * nyz + j * nz + k;
    b[j] = c * psi[ind] - (2.0 * psi[ind - nz] - 2.0 * psi[ind]);
    
    grid_solve_tridiagonal_system(ny, dl, d, du, b, dl);

    for(j = 0; j < ny; j++)
      psi[i * nyz + j * nz + k] = dl[j];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy using periodic boundary condition (Crank-Nicolson).
 * 
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 * 
 */

EXPORT void grid3d_wf_propagate_kinetic_cn_pbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  /*
   * exp( -i (Tx + Ty + Tz) dt / hbar ) 
   *   = exp( -i Tx dt / hbar ) exp( -i Ty dt / hbar ) exp( -i Tz dt / hbar ) + O(dt^2)
   *   
   */

  grid3d_wf_propagate_kinetic_x_cn_pbc(gwf, time, workspace);
  grid3d_wf_propagate_kinetic_y_cn_pbc(gwf, time, workspace);
  grid3d_wf_propagate_kinetic_z_cn_pbc(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along x using periodic boundary condition (Crank-Nicolson).
 * 
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 * 
 */

EXPORT void grid3d_wf_propagate_kinetic_x_cn_pbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, jk;
  long nz = gwf->grid->nz, nyz = ny * nz;
  double complex *du, *dl, *d, *b, *x, *z;
  double complex gamma, alpha, beta, fact;

  /* Periodic */
  alpha = 1.0; beta = 1.0;   // change to -1 to get antisymm.
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nx,nz,nyz,psi,wrk,c,alpha,beta) private(tid,i,j,k,jk,dl,d,du,b,ind,x,z,gamma,fact) default(none) schedule(runtime)
  for (jk = 0; jk < nyz; jk++) {  /* for each (y,z) */
    j = jk / nz;
    k = jk % nz;
    tid = omp_get_thread_num();
    dl = &wrk[nx * (6 * tid + 0)];   /* 6 = # of arrays: dl, d, du, b, x, z */
    d  = &wrk[nx * (6 * tid + 1)];
    du = &wrk[nx * (6 * tid + 2)];
    b  = &wrk[nx * (6 * tid + 3)];
    x  = &wrk[nx * (6 * tid + 4)];
    z  = &wrk[nx * (6 * tid + 5)];

    /* create right-hand side vector */
    for(i = 1; i < nx - 1; i++) {
      ind = i * nyz + j * nz + k;
      /* (C - dx^2 laplace) psi(t) */
      b[i] = c * psi[ind] - (psi[ind + nyz] - 2.0 * psi[ind] + psi[ind - nyz]);
    }
    i = 0; ind = i * nyz + j * nz + k;     /* periodic boundaries */
    b[i] = c * psi[ind] - (psi[ind + nyz] - 2.0 * psi[ind] + beta * psi[ind + (nx-1) * nyz]);
    i = nx - 1; ind = i * nyz + j * nz + k;
    b[i] = c * psi[ind] - (psi[ind - nyz] - 2.0 * psi[ind] + alpha * psi[ind - (nx-1) * nyz]);
    
    /* Sherman-Morrison (See Num. Recp.) */
    for(i = 0; i < nx; i++) {
      dl[i] = 1.0;
      du[i] = 1.0;
      d[i] = ca;
    }
    gamma = -d[0];
    d[0] -= gamma;
    d[nx-1] -= alpha * beta / gamma;
    grid_solve_tridiagonal_system(nx, dl, d, du, b, x); // d & b overwritten   (debug last dl was x)
    for (i = 0; i < nx; i++) {
      b[i] = 0.0;
      dl[i] = 1.0;
      du[i] = 1.0;
      d[i] = ca;
    }
    b[0] = gamma;
    b[nx-1] = alpha;
    d[0] -= gamma;
    d[nx-1] -= alpha * beta / gamma;
    grid_solve_tridiagonal_system(nx, dl, d, du, b, z); // d and b overwritten
    fact = (x[0] + beta * x[nx-1] / gamma) / (1.0 + z[0] + beta * z[nx-1] / gamma);
    for(i = 0; i < nx; i++)
      psi[i * nyz + j * nz + k] = x[i] - fact * z[i];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along y using periodic boundary condition (Crank-Nicolson).
 * 
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 * 
 */

EXPORT void grid3d_wf_propagate_kinetic_y_cn_pbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, ik;
  long nz = gwf->grid->nz, nxz = nx * nz, nyz = ny * nz;
  double complex *du, *dl, *d, *b, *x, *z;
  double complex gamma, alpha, beta, fact;

  /* Periodic */
  alpha = 1.0; beta = 1.0;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nx,ny,nz,nxz,nyz,psi,wrk,c,alpha,beta) private(tid,i,j,k,ik,dl,d,du,b,ind,x,z,gamma,fact) default(none) schedule(runtime)
  for (ik = 0; ik < nxz; ik++) {  /* for each (x,z) */
    i = ik / nz;
    k = ik % nz;
    tid = omp_get_thread_num();
    dl = &wrk[ny * (6 * tid + 0)];
    d  = &wrk[ny * (6 * tid + 1)];
    du = &wrk[ny * (6 * tid + 2)];
    b  = &wrk[ny * (6 * tid + 3)];
    x  = &wrk[ny * (6 * tid + 4)];
    z  = &wrk[ny * (6 * tid + 5)];

    /* create right-hand side vector */
    for(j = 1; j < ny - 1; j++) {
      ind = i * nyz + j * nz + k;
      /* (C - dx^2 laplace) psi(t) */
      b[j] = c * psi[ind] - (psi[ind + nz] - 2.0 * psi[ind] + psi[ind - nz]);
    }
 
    /* periodic boundaries */
    j = 0; ind = i * nyz + j * nz + k;
    b[j] = c * psi[ind] - (psi[ind + nz] - 2.0 * psi[ind] + alpha * psi[ind + (ny-1) * nz]);
    j = ny - 1; ind = i * nyz + j * nz + k;
    b[j] = c * psi[ind] - (psi[ind - nz] - 2.0 * psi[ind] + beta * psi[ind - (ny-1) * nz]);
    
    /* Sherman-Morrison (See Num. Recp.) */
    for(j = 0; j < ny; j++) {
      dl[j] = 1.0;
      du[j] = 1.0;
      d[j] = ca;
    }
    gamma = -d[0];
    d[0] -= gamma;
    d[ny-1] -= alpha * beta / gamma;
    grid_solve_tridiagonal_system(ny, dl, d, du, b, x);
    for (j = 0; j < ny; j++) {
      b[j] = 0.0;
      dl[j] = 1.0;
      du[j] = 1.0;
      d[j] = ca;
    }      
    b[0] = gamma;
    b[ny-1] = alpha;
    d[0] -= gamma;
    d[ny-1] -= alpha * beta / gamma;
    grid_solve_tridiagonal_system(ny, dl, d, du, b, z);
    fact = (x[0] + beta * x[ny-1] / gamma) / (1.0 + z[0] + beta * z[ny-1] / gamma);
    for(j = 0; j < ny; j++)
      psi[i * nyz + j * nz + k] = x[j] - fact * z[j];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along z using periodic boundary condition (Crank-Nicolson).
 * 
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 * 
 */

EXPORT void grid3d_wf_propagate_kinetic_z_cn_pbc(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  long k, ij;
  long nz = gwf->grid->nz, nyz = ny * nz, nxy = nx * ny;
  double complex *du, *dl, *d, *b, *x, *z;
  double complex gamma, alpha, beta, fact;
  
  /* Periodic */
  alpha = 1.0; beta = 1.0;

  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nx,nz,nyz,nxy,psi,wrk,c,alpha,beta) private(tid,i,j,k,ij,dl,d,du,b,ind,x,z,gamma,fact) default(none) schedule(runtime)
  for (ij = 0; ij < nxy; ij++) {  /* for each (x,y) */
    i = ij / nz;
    j = ij % nz;
    tid = omp_get_thread_num();
    dl = &wrk[nz * (6 * tid + 0)];
    d  = &wrk[nz * (6 * tid + 1)];
    du = &wrk[nz * (6 * tid + 2)];
    b  = &wrk[nz * (6 * tid + 3)];
    x  = &wrk[nz * (6 * tid + 4)];
    z  = &wrk[nz * (6 * tid + 5)];

    /* create right-hand side vector */
    for(k = 1; k < nz - 1; k++) {
      ind = i * nyz + j * nz + k;
      /* (C - dx^2 laplace) psi(t) */
      b[k] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + psi[ind - 1]);
    }
 
    /* periodic boundaries */
    k = 0; 
    ind = i * nyz + j * nz + k;
    b[k] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + alpha * psi[ind + (nz-1)]);
    k = nz - 1;
    ind = i * nyz + j * nz + k;
    b[k] = c * psi[ind] - (psi[ind - 1] - 2.0 * psi[ind] + beta * psi[ind - (nz-1)]);
    
    /* Sherman-Morrison (See Num. Recp.) */
    for(k = 0; k < nz; k++) {
      dl[k] = 1.0;
      du[k] = 1.0;
      d[k] = ca;
    }
    gamma = -d[0];
    d[0] -= gamma;
    d[nz-1] -= alpha * beta / gamma;
    grid_solve_tridiagonal_system(nz, dl, d, du, b, x);
    for (k = 0; k < nz; k++) {
      b[k] = 0.0;
      dl[k] = 1.0;
      du[k] = 1.0;
      d[k] = ca;
    }
    b[0] = gamma;
    b[nz-1] = alpha;
    d[0] -= gamma;
    d[nz-1] -= alpha * beta / gamma;
    grid_solve_tridiagonal_system(nz, dl, d, du, b, z);
    fact = (x[0] + beta * x[nz-1] / gamma) / (1.0 + z[0] + beta * z[nz-1] / gamma);
    for(k = 0; k < nz; k++)
      psi[i * nyz + j * nz + k] = x[k] - fact * z[k];
  }
}

/*
 * Produce density grid from a given wavefunction.
 *
 * gwf     = wavefunction (wf3d *).
 * density = output density grid (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT inline void grid3d_wf_density(const wf3d *gwf, rgrid3d *density) {

  long ij, k, ijnz, nxy = gwf->grid->nx * gwf->grid->ny, nz = gwf->grid->nz;
  double complex *avalue = gwf->grid->value;
  double *cvalue = density->value;
  
#pragma omp parallel for firstprivate(nxy,nz,avalue,cvalue) private(ij,ijnz,k) default(none) schedule(runtime)
  for(ij = 0; ij < nxy; ij++) {
    ijnz = ij * nz;
    for(k = 0; k < nz; k++)
      cvalue[ijnz + k] = conj(avalue[ijnz + k]) * avalue[ijnz + k];
  }
}

/*
 * Zero wavefunction.
 *
 * gwf = wavefunction to be zeroed (wf3d *).
 *
 * No return value.
 *
 */

EXPORT inline void grid3d_wf_zero(wf3d *gwf) { 

  cgrid3d_zero(gwf->grid); 
}

/*
 * Set wavefunction to some constant value.
 *
 * gwf = wavefunction to be set (wf3d *).
 * c   = value (double complex).
 *
 * No return value.
 *
 */

EXPORT inline void grid3d_wf_constant(wf3d *gwf, double complex c) { 

  cgrid3d_constant(gwf->grid, c); 
}

/*
 * Map a given function on a wavefunction.
 *
 * gwf  = wavefunction where function will be mapped to (wf3d *).
 * func = function providing the mapping (double complex (*)(void *, double)).
 * farg = optional argument for passing parameters to func (void *).
 *
 * No return value.
 *
 */

EXPORT inline void grid3d_wf_map(wf3d *gwf, double complex (*func)(void *arg, double x, double y, double z), void *farg) { 

  cgrid3d_map(gwf->grid, func, farg); 
}

/*
 * Calculate the norm of the given wavefunction.
 *
 * gwf = wavefunction for the calculation (wf3d *).
 *
 * Returns the norm (double).
 *
 */

EXPORT inline double grid3d_wf_norm(const wf3d *gwf) { 

  return cgrid3d_integral_of_square(gwf->grid); 
}

/*
 * Normalize wavefunction (to the value given in gwf->norm).
 *
 * gwf = wavefunction to be normalized (wf3d *).
 *
 * Returns the normalization constant (double).
 *
 */

EXPORT inline double grid3d_wf_normalize(wf3d *gwf) { 

  double norm = grid3d_wf_norm(gwf);

  cgrid3d_multiply(gwf->grid, sqrt(gwf->norm / norm));
  return norm; 
}

/*
 * Calculate overlap between two wavefunctions.
 *
 * gwfa = 1st wavefunction (wf3d *).
 * gwfb = 2nd wavefunction (wf3d *).
 *
 * Returns the overlap (double complex).
 *
 */

EXPORT inline double complex grid3d_wf_overlap(const wf3d *gwfa, const wf3d *gwfb) { 

  return cgrid3d_integral_of_conjugate_product(gwfa->grid, gwfb->grid); 
}

/*
 * Output wavefunction.
 *
 * gwf = wavefunction to be printed (wf3d *).
 * out = output file pointer (FILE *).
 *
 * No return value.
 *
 */

EXPORT inline void grid3d_wf_print(const wf3d *gwf, FILE *out) { 

  cgrid3d_print(gwf->grid, out); 
}
