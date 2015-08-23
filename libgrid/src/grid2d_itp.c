/*
 * Routines to find eigenfunctions by the imaginary time propagation method in 2D.
 *
 * To debug define __DEBUG__.
 *
 */

#include "grid.h"
#include "private.h"

/*
 * Solve eigenfunctions of a Hamiltonian with a linear potential function.
 *
 * gwf            = an array of dimension "states" for storing the resulting eigenfunctions (wf2d **).
 * states         = number of states requested (int).
 * virtuals       = number of virtual states for orthogonalization (int).
 * potential      = potential grid (cgrid2d *).
 * tau            = initial imaginary time step length (double). This will be adjusted to smaller values dynamically.
 *                  If tau is given as negative number, its absolute value will be used and the time step will not
 *                  be adjusted dynamically.
 * threshold      = convergence threshold (dE / E < threshold) (double).
 * max_iterations = maximum number of iterations allowed (int).
 * rtau           = the final (adjusted) imaginary time step (double *).
 * riterations    = final number of iterations (int *).
 * 
 * Return value is the relative error (dE / E).
 *
 */

EXPORT double grid2d_itp_linear(wf2d **gwf, int states, int virtuals, const cgrid2d *potential, double tau, double threshold, int max_iterations, double *rtau, int *riterations) {

  int i, l;
  double cfactor = 0.5, max_cfactor = 0.9;
  double erel = 0.0, erms_long, erms_short, derms_long, derms_short;
  double complex time_step;
  double *energy_long = 0, *energy_short = 0;
  double *error_long = 0, *error_short = 0;
  double error; 

  wf2d **gwf_short;
  wf2d **gwf_long = gwf;
  
  cgrid2d *sq_grad_pot = 0;
  cgrid2d *workspace  = 0;
  cgrid2d *workspace2 = 0;
    
  /* if tau is negative, use constant time step */
  if (tau < 0.0) {
    cfactor = 1.0; tau = -tau;
  }
  
  /* allocate another set of wave functions for workspace */
  gwf_short = (wf2d **) malloc(states * sizeof(wf2d *));
  for(i = 0; i < states; i++) {
    gwf_short[i] = grid2d_wf_alloc(gwf[i]->grid->nx, gwf[i]->grid->ny, gwf[i]->grid->step,
				 gwf[i]->mass, gwf[i]->boundary, gwf[i]->propagator);
    
    if (!gwf_short[i]) {
      fprintf(stderr, "libgrid: Error in grid2d_itp_linear(). Could not allocate memory for workspaces.\n");
      abort();
    }
  }
  
  /* allocate grids for |grad pot|^2 and workspaces */
  sq_grad_pot = cgrid2d_alloc(potential->nx, potential->ny, potential->step,
			     potential->value_outside, potential->outside_params_ptr);
  workspace = cgrid2d_alloc(potential->nx, potential->ny, potential->step,
			   potential->value_outside, potential->outside_params_ptr);
  workspace2 = cgrid2d_alloc(potential->nx, potential->ny, potential->step,
			    potential->value_outside, potential->outside_params_ptr);
  
  if (!sq_grad_pot || !workspace || !workspace2) {
    fprintf(stderr, "libgrid: Error in grid2d_itp_linear(). Could not allocate memory for workspaces.\n");
    abort();
  }
  
  energy_long = (double *) malloc(states * sizeof(double));
  error_long = (double *) malloc(states * sizeof(double));
  energy_short = (double *) malloc(states * sizeof(double));
  error_short = (double *) malloc(states * sizeof(double));
  if (!energy_long || !energy_short || !error_short || !error_long) {
    fprintf(stderr, "libgrid: Error in grid2d_itp_linear(). Could not allocate memory for workspaces.\n");
    abort();
  }
  
  /* copy wave functions */
  for(i = 0; i < states; i++)
    cgrid2d_copy(gwf_short[i]->grid, gwf_long[i]->grid);
  
  /* |grad potential|^2 */
  grid2d_wf_square_of_potential_gradient(sq_grad_pot, potential, workspace);
  
  /* iteration loop */
  for(l = 0; l < max_iterations; l++) {
    
    /* propagate t = - i tau */
    time_step = -I * tau;
    for(i = 0; i < states; i++)
      grid2d_wf_propagate(gwf_long[i], potential, sq_grad_pot, time_step, workspace);
    
    grid2d_wf_diagonalize(gwf_long, states);
    
    for(i = 0; i < states; i++)
      energy_long[i] = grid2d_wf_energy_and_error(gwf_long[i], potential, workspace, &error_long[i]); 
    
    /* if constant time step, skip time step adjusting */
    if (cfactor > max_cfactor) {
#if __DEBUG__
      fprintf(stderr, "El %d ", l);
      for(i = 0; i < states; i++)
        fprintf(stderr, "%20.15le ", energy_long[i]);
      fprintf(stderr, "\n");
      fprintf(stderr, "el %d ", l);
      for(i = 0; i < states; i++)
        fprintf(stderr, "%20.15le ", error_long[i]);
      fprintf(stderr, "\n");
#endif      
      continue;
    }
    
    /* propagate t = - i c tau */
    time_step = -I * cfactor * tau;
    for(i = 0; i < states; i++)
      grid2d_wf_propagate(gwf_short[i], potential, sq_grad_pot, time_step, workspace);
    
    grid2d_wf_diagonalize(gwf_short, states);
    
    for( i = 0; i < states; i++ ) 
      energy_short[i] = grid2d_wf_energy_and_error(gwf_short[i], potential, workspace, &error_short[i]);
    
    /* relative error, dE / E */
    erel = 0;
    for(i = 0; i < states - virtuals; i++)
      erel += 2.0 * (error_long[i] / energy_long[i]) * (error_long[i] / energy_long[i]);
    erel = sqrt(erel);
    
    /* check convergence */
    if (erel < threshold) break;
    
    /* rms of absolute energy and error */
    erms_long = erms_short = 0.0;
    derms_long = derms_short = 0.0;
    for(i = 0; i < states - virtuals; i++) {
      erms_long += energy_long[i]  * energy_long[i];
      erms_short += energy_short[i] * energy_short[i];
      
      derms_long  += error_long[i]  * error_long[i];
      derms_short += error_short[i] * error_short[i];
    }
    erms_long = sqrt(erms_long / (states - virtuals));
    erms_short = sqrt(erms_short / (states - virtuals));
    derms_long = sqrt(derms_long / (states - virtuals));
    derms_short = sqrt(derms_short / (states - virtuals));
    
    /* if long time step gives better energy or error, use it and corresponding wave function for next iteration */
    if (erms_long < erms_short || derms_long < derms_short) {
      for(i = 0; i < states; i++)
        cgrid2d_copy(gwf_short[i]->grid, gwf_long[i]->grid);
      
      /* try smaller time step reduce */
      cfactor = sqrt(cfactor);
      if (cfactor > max_cfactor) cfactor = max_cfactor;
    }
    /* else use short time step and corresponing wave function */
    else {
      for(i = 0; i < states; i++)
        cgrid2d_copy(gwf_long[i]->grid, gwf_short[i]->grid);
      
      /* try shorter time step */
      tau = tau * cfactor;
      /* and larger time step reduce */
      cfactor *= cfactor;
    }
    
#if __DEBUG__
    fprintf(stderr, "T    %d  %lf\n", l, tau);
    fprintf(stderr, "dE   %d  %le\n", l, erel);
    
    fprintf(stderr, "El   %d ", l);
    for(i = 0; i < states; i++)
      fprintf(stderr, "%20.15le ", energy_long[i]);
    fprintf(stderr, "\n");
    
    fprintf(stderr, "Errl %d ", l);
    for(i = 0; i < states; i++)
      fprintf(stderr, "%20.15le ", error_long[i]);
    fprintf(stderr, "\n");
#endif
  }
#if __DEBUG__
  fprintf(stderr, "\n");
#endif
  
  /* free workspaces */
  for(i = 0; i < states; i++)
    grid2d_wf_free(gwf_short[i]);
  free(gwf_short);
  
  cgrid2d_free(sq_grad_pot);
  cgrid2d_free(workspace);
  cgrid2d_free(workspace2);
  
  free(energy_long);
  free(error_long);
  free(energy_short);
  free(error_short);

  error = erel;
  *rtau = tau;
  *riterations = l;
  
  return error;
}

/*
 * Solve eigenfunctions of a Hamiltonian with a nonlinear potential function.
 *
 * gwf                  = an array of dimension "states" for storing the resulting eigenfunctions (wf2d **).
 * states               = number of states requested (int).
 * virtuals             = number of virtual states for orthogonalization (int).
 * calculate_potentials = nonlinear potential, which takes the current grid etc. as argument (void (*)(cgrid2d **, void *, wf2d **, int)).
 *                        Note that this is a pointer to the potential function.
 * tau                  = initial imaginary time step length (double). This will be adjusted to smaller values dynamically.
 *                        If tau is given as negative number, its absolute value will be used and the time step will not
 *                        be adjusted dynamically.
 * threshold            = convergence threshold (dE / E < threshold) (double).
 * max_iterations       = maximum number of iterations allowed (int).
 * rtau                 = the final (adjusted) imaginary time step (double *).
 * riterations          = final number of iterations (int *).
 * 
 * Return value is the relative error (dE / E).
 *
 */

EXPORT double grid2d_itp_nonlinear(wf2d **gwf, int states, int virtuals, void (*calculate_potentials)(cgrid2d **potential, void *arg, wf2d **gwf, int states), void *arg, double tau, double threshold, int max_iterations, double *rtau, int *riterations) {

  int i, l;
  double erel = 0.0, erms_long, erms_short, derms_long, derms_short;
  double cfactor = 0.5, max_cfactor = 0.9;
  double complex time_step;
  double error;
  
  wf2d **gwf_short;
  wf2d **gwf_long = gwf;
  
  cgrid2d **potential   = 0;
  cgrid2d **sq_grad_pot = 0;
  cgrid2d *workspace  = 0;
  cgrid2d *workspace2 = 0;
  
  double *energy_long = 0, *energy_short = 0;
  double *error_long = 0, *error_short = 0;
  
  /* if tau is negative, use constant time step */
  if (tau < 0.0) {
    cfactor = 1.0; tau = -tau;
  }
  
  /* allocate grids for another set of wave functions, potential, |grad pot|^2 and workspaces */
  gwf_short = (wf2d **) malloc(states * sizeof(wf2d *));
  potential = (cgrid2d **) malloc(states * sizeof(cgrid2d *));
  sq_grad_pot = (cgrid2d **) malloc(states * sizeof(cgrid2d *));
  
  for(i = 0; i < states; i++) {
    gwf_short[i] = grid2d_wf_alloc(gwf[i]->grid->nx, gwf[i]->grid->ny, gwf[i]->grid->step,
				   gwf[i]->mass, gwf[i]->boundary, gwf[i]->propagator);
    potential[i] = cgrid2d_alloc(gwf[i]->grid->nx, gwf[i]->grid->ny, gwf[i]->grid->step,
				gwf[i]->grid->value_outside, gwf[i]->grid->outside_params_ptr);
    sq_grad_pot[i] = cgrid2d_alloc(potential[i]->nx, potential[i]->ny, potential[i]->step,
				  potential[i]->value_outside, potential[i]->outside_params_ptr);
  }
  workspace = cgrid2d_alloc(potential[0]->nx, potential[0]->ny, potential[0]->step,
			   potential[0]->value_outside, potential[0]->outside_params_ptr);
  workspace2 = cgrid2d_alloc(potential[0]->nx, potential[0]->ny, potential[0]->step,
			    potential[0]->value_outside, potential[0]->outside_params_ptr);
  
  for(i = 0; i < states; i++)
    if (!gwf_short[i] || !potential[i] || !sq_grad_pot[i]) {
      fprintf(stderr, "libgrid: Error in grid2d_itp_nonlinear(). Could not allocate memory for workspaces.\n");
      abort();
    }
  
  if (!workspace || !workspace2) {
    fprintf(stderr, "libgrid: Error in grid2d_itp_nonlinear(). Could not allocate memory for workspaces.\n");
    abort();
  }
  
  energy_long = (double *) malloc(states * sizeof(double));
  error_long = (double *) malloc(states * sizeof(double));
  energy_short = (double *) malloc(states * sizeof(double));
  error_short = (double *) malloc(states * sizeof(double));
  if (!energy_long || !energy_short || !error_short || !error_long) {
    fprintf(stderr, "libgrid: Error in grid2d_itp_nonlinear(). Could not allocate memory for workspaces.\n");
    abort();
  }
  
  /* copy wave functions */
  for(i = 0; i < states; i++)
    cgrid2d_copy(gwf_short[i]->grid, gwf_long[i]->grid);
  
  /* iteration loop */
  for(l = 0; l < max_iterations; l++) {
    
    /* calculate potentials */
    (*calculate_potentials)(potential, arg, gwf_long, states);
    
    /* |grad potential|^2 */
    for(i = 0; i < states; i++)
      grid2d_wf_square_of_potential_gradient(sq_grad_pot[i], potential[i], workspace);
    
    /* propagate t = - i tau */
    time_step = -I * tau;
    for(i = 0; i < states; i++)
      grid2d_wf_propagate(gwf_long[i], potential[i], sq_grad_pot[i], time_step, workspace);
    
    grid2d_wf_diagonalize(gwf_long, states);
    
    for(i = 0; i < states; i++)
      energy_long[i] = grid2d_wf_energy_and_error(gwf_long[i], potential[i], workspace, &error_long[i]);
    
    /* if constant time step, skip step adjusting */
    if (cfactor > max_cfactor) {
#if __DEBUG__
      fprintf(stderr,"E  %4d ", l);
      for(i = 0; i < states; i++)
        fprintf(stderr," %20.15lf %20.15lf ", energy_long[i], error_long[i]);
      fprintf(stderr,"\n");
#endif
      continue;
    }
    
    /* propagate t = - i c tau */
    time_step = -I * cfactor * tau;
    for(i = 0; i < states; i++)
      grid2d_wf_propagate(gwf_short[i], potential[i], sq_grad_pot[i], time_step, workspace);
    
    grid2d_wf_diagonalize(gwf_short, states);
    
    for(i = 0; i < states; i++)
      energy_short[i] = grid2d_wf_energy_and_error(gwf_short[i], potential[i], workspace, &error_short[i]);
    
    /* max relative error, dE^2 / E^2 */
    erel = 0;
    for(i = 0; i < states - virtuals; i++)
      erel += 2.0 * (error_long[i] / energy_long[i]) * (error_long[i] / energy_long[i]);
    
    /* check convergence */
    if ( sqrt(erel) < threshold ) break;
   
    /* rms of energy absolute error */
    erms_long = erms_short = 0;
    derms_long = derms_short = 0;
    for( i = 0; i < states - virtuals; i++ ) {
      erms_long += energy_long[i] * energy_long[i];
      erms_short += energy_short[i] * energy_short[i];
      derms_long += error_long[i] * error_long[i];
      derms_short += error_short[i] * error_short[i];
    }
    erms_long = sqrt(erms_long / (states - virtuals));
    erms_short = sqrt(erms_short / (states - virtuals));
    derms_long = sqrt(derms_long / (states - virtuals));
    derms_short = sqrt(derms_short / (states - virtuals));
    
    /* if short time step gives better energy AND error, use it and corresponding wave function for next iteration */
    if (erms_short < erms_long && derms_short < derms_long) {
      for(i = 0; i < states; i++)
        cgrid2d_copy(gwf_long[i]->grid, gwf_short[i]->grid);
      
      /* try shorter time step */
      tau = tau * cfactor;
      /* and larger time step reduce */
      cfactor *= cfactor;
    }
    /* else use long time step and corresponing wave function */
    else {
      for(i = 0; i < states; i++)
        cgrid2d_copy(gwf_short[i]->grid, gwf_long[i]->grid);
      
      /* try smaller time step reduce */
      cfactor = sqrt(cfactor);
      if (cfactor > max_cfactor) cfactor = max_cfactor;
    }
    
#if __DEBUG__
    fprintf(stderr,"T   %4d %12.4lf %12.4le\n", l, tau, derms_long);
    fprintf(stderr,"El  %4d ", l);
    for(i = 0; i < states; i++)
      fprintf(stderr," %20.15lf %20.15lf ", energy_long[i], error_long[i]);
    fprintf(stderr,"\n");
    fprintf(stderr,"Es  %4d ", l);
    for(i = 0; i < states; i++)
      fprintf(stderr," %20.15lf %20.15lf ", energy_short[i], error_short[i]);
    fprintf(stderr,"\n" );
#endif
  }
  
  /* free workspaces */
  for(i = 0; i < states; i++) {
    grid2d_wf_free(gwf_short[i]);
    cgrid2d_free(potential[i]);
    cgrid2d_free(sq_grad_pot[i]);
  }
  free(gwf_short);
  free(potential);
  free(sq_grad_pot);
  
  cgrid2d_free(workspace);
  cgrid2d_free(workspace2);
  
  free(energy_long);
  free(error_long);
  free(energy_short);
  free(error_short);

  error = erel;
  *rtau = tau;
  *riterations = l;
  
  return error;
}
