/*
 * Additional routines for propagating in cylindrical coordinates (2D).
 * (Neumann boundaries)
 *
 */

#include "grid.h"
#include "private.h"
#include "private2d.h"

/*
 * Routine for propagating kinetic energy in cylindrical coordinates (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_cn_cyl(wf2d *gwf, double complex time, cgrid2d *workspace) {

  /*
   * exp(-i (Tx + Ty) dt / hbar) 
   *   = exp(-i Tx dt / hbar) exp(-i Ty dt / hbar) + O(dt^2)
   *   
   */

  grid2d_wf_propagate_kinetic_cn_cyl_z(gwf, time, workspace);
  grid2d_wf_propagate_kinetic_cn_cyl_r(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along z in cylindrical coordinates (Crank-Nicolson). Neumann BC.
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 *
 * NOTE: In principle one could use the private.h cholesky with neumann
 * boundaries but it does not seem to work quite correctly. Thus the 
 * code below uses LAPACK (just line along r further down).
 */

EXPORT void grid2d_wf_propagate_kinetic_cn_cyl_z(wf2d *gwf, double complex time, cgrid2d *workspace) {

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
  for (j = 0; j < ny; j++) {  /* for each r */
    tid = omp_get_thread_num();
    dl = &wrk[nx * (4 * tid + 0)];
    d  = &wrk[nx * (4 * tid + 1)];
    du = &wrk[nx * (4 * tid + 2)];
    b  = &wrk[nx * (4 * tid + 3)];

    for(i = 0; i < nx; i++) {     /* over z */
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
 * Auxiliary routine for propagating kinetic energy along r in cylindrical coordinates (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf2d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid2d *).
 *
 * No return value.
 *
 */

EXPORT void grid2d_wf_propagate_kinetic_cn_cyl_r(wf2d *gwf, double complex time, cgrid2d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step, r;
  long tid, ind, iny;
  long i, nx = gwf->grid->nx;
  long j, ny = gwf->grid->ny;
  double complex *dl, *d, *du, *b;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt hbar) psi(t)  <=>  A x = b
   * (C + dr^2 laplace + dr^2 grad/r) psi(t+dt) = (C - dr^2 laplace - dr^2 grad/r) psi(t)
   * where C = 4 i m dx^2 / (hbar dt)
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);
  
  /* (C + dr^2 laplace + (1/r) dr^2 grad) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nx,ny,psi,wrk,c,step) private(tid,i,j,dl,d,du,b,ind,r,iny) default(none) schedule(runtime)
  for(i = 0; i < nx; i++) { /* for each z */
    tid = omp_get_thread_num();
    dl = &wrk[ny * (4 * tid + 0)];
    d  = &wrk[ny * (4 * tid + 1)];
    du = &wrk[ny * (4 * tid + 2)];
    b  = &wrk[ny * (4 * tid + 3)];

    for(j = 0; j < ny; j++) { /* over r */
      r = step * (double) j;
      if(j == ny-1) dl[ny-1] = 2.0;            /* Neumann outer boundary */
      else if(j > 0) dl[j] = 1.0 - 0.5 * step / r;       
      if(j == 0) du[0] = 2.0;                  /* Neumann inner boundary */
      else if(j < ny-1) du[j] = 1.0 + 0.5 * step / r;      
      d[j] = ca;
    }

    /* create right-hand side vector */
    iny = i * ny;
    for(j = 1; j < ny - 1; j++) {
      ind = iny + j;
      r = step * (double) j;
      /* (C - dr^2 laplace - (1/r) dr^2 grad) psi(t) */
      b[j] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + psi[ind - 1]) - (0.5 * step / r) * (psi[ind + 1] - psi[ind - 1]);
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
 * Map a given function on a wavefunction.
 *
 * gwf  = wavefunction where function will be mapped to (wf2d *).
 * func = function providing the mapping (double complex (*)(void *, double)).
 * farg = optional argument for passing parameters to func (void *).
 *
 * No return value.
 *
 */

EXPORT inline void grid2d_wf_map_cyl(wf2d *gwf, double complex (*func)(void *arg, double z, double r), void *farg) { 

  cgrid2d_map_cyl(gwf->grid, func, farg); 
}

/*
 * Calculate the norm of the given wavefunction (2-D cylindrical coordinates).
 *
 * gwf = wavefunction for the calculation (wf2d *).
 *
 * Returns the norm (double).
 *
 */

EXPORT double grid2d_wf_norm_cyl(const wf2d *gwf) { 

  return cgrid2d_integral_of_square_cyl(gwf->grid); 
}

/*
 * Normalize a given wavefunction (2-D cylindrical coordinates).
 *
 * gwf = wavefuction to be normalized (wf2d *).
 *
 * Returns the normalization factor used.
 *
 */

EXPORT inline double grid2d_wf_normalize_cyl(wf2d *gwf) { 

  double norm = grid2d_wf_norm_cyl(gwf);

  norm = sqrt(gwf->norm / norm);
  cgrid2d_multiply(gwf->grid, norm);
  return norm; 
}

/*
 * Calculate energy for a given wavefunction (2-D cylindrical coords)
 *
 * gwf        = wavefunction for kinetic energy evaluation (wf2d *).
 * potential  = potential energy to be included (if any) (cgrid2d *).
 * cworkspace = complex workspace required (cgrid2d *).
 *
 * Returns the kinetic energy.
 *
 */

EXPORT double grid2d_wf_energy_cyl(const wf2d *gwf, const cgrid2d *potential, cgrid2d *cworkspace) {

  long i;

  /* (-2m/hbar^2) T psi */
  cgrid2d_fd_laplace_cyl(gwf->grid, cworkspace);
  cgrid2d_multiply(cworkspace, -HBAR * HBAR / (2.0 * gwf->mass));

  /* V psi V */
  if(potential) 
    for(i = 0; i < gwf->grid->nx * gwf->grid->ny; i++)
      cworkspace->value[i] += potential->value[i] * gwf->grid->value[i];
  
  /* int psi^* (T + V) psi d^3r */
  return creal(cgrid2d_integral_of_conjugate_product_cyl(gwf->grid, cworkspace));
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

EXPORT void grid2d_wf_absorb_cyl(cgrid2d *potential, rgrid2d *density, double rho0, double (*region)(void *, double, double), rgrid2d *workspace) {

  rgrid2d_copy(workspace, density);
  rgrid2d_add(workspace, -rho0);
  rgrid2d_product_func_cyl(workspace, region, NULL);
  rgrid2d_multiply(workspace, -1.0);
  grid2d_add_real_to_complex_im(potential, workspace);
}
