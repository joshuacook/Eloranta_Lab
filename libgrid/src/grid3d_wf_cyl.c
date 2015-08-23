/*
 * Routines for handling 3D wavefunctions.
 *
 */

#include "grid.h"
#include "private.h"
#include "private3d.h"


/*
 * Calculate energy for the wavefunction (3D cyl).
 *
 * gwf       = wavefunction for the energy calculation (wf3d *).
 * potential = grid containing the potential (cgrid3d *).
 * workspace = additional storage needed (cgrid3d *).
 *
 * Returns the energy (double).
 *
 */

EXPORT double grid3d_wf_energy_cyl(const wf3d *gwf, const cgrid3d *potential, cgrid3d *workspace) {

  return grid3d_wf_energy_cn_cyl(gwf, gwf, potential, workspace);
}

/*
 * Auxiliary routine for calculating the energy (Crank-Nicolson). 3D Cylindrical coords.
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

EXPORT double grid3d_wf_energy_cn_cyl(const wf3d *gwfa, const wf3d *gwfb, const cgrid3d *potential, cgrid3d *workspace) {
  
  long i;

  /* (-2m/hbar^2) T psi */
  cgrid3d_fd_laplace_cyl(gwfb->grid, workspace);
  cgrid3d_multiply(workspace, -HBAR * HBAR / (2.0 * gwfb->mass));

  /* V psi */
  if(potential)
    for(i = 0; i < gwfb->grid->nx * gwfb->grid->ny * gwfb->grid->nz; i++)
      workspace->value[i] += potential->value[i] * gwfb->grid->value[i];
  
  /* int psi^* (T + V) psi d^3r */
  return creal(cgrid3d_integral_of_conjugate_product_cyl(gwfa->grid, workspace));
}

/*
 * Auxiliary routine for calculating potential energy (3D cyl).
 * 
 * gwf       = wavefunction for the kinetic energy calculation (wf3d *).
 * workspace = additional workspace required for the operation (cgrid3d *).
 *
 * Returns the potential energy.
 *
 */

EXPORT double grid3d_wf_potential_energy_cyl(const wf3d *gwf, const cgrid3d *potential) {

  return creal(cgrid3d_grid_expectation_value_cyl(gwf->grid, potential));
}

/*
 * Routine to propagate kinetic energy (Crank-Nicolson). 3D cyl.
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional workspace required for the operation (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_cn_cyl(wf3d *gwf, double complex time, cgrid3d *workspace) {

  /*
   * exp(-i (Tx + Ty) dt / hbar) 
   *   = exp(-i Tx dt / hbar) exp(-i Ty dt / hbar) + O(dt^2)
   *   
   */

  grid3d_wf_propagate_kinetic_cn_cyl_r(gwf, time, workspace);
  grid3d_wf_propagate_kinetic_cn_cyl_phi(gwf, time, workspace);
  // debug
  //grid3d_wf_propagate_kinetic_cn_cyl_z(gwf, time, workspace);
}

/*
 * Auxiliary routine for propagating kinetic energy along r using Neumann boundary condition (Crank-Nicolson).
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_cn_cyl_r(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, nr = gwf->grid->nx;
  long j, nphi = gwf->grid->ny;
  long k, nz = gwf->grid->nz;
  long nphiz = nphi * nz, jk;
  double complex *du, *dl, *d, *b, r;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* (C + dr^2 laplace + dr/(2r) x d/dr) psi(t+dt) */
  ca = c - 2.0;
#pragma omp parallel for firstprivate(ca,nz,nr,nphi,nphiz,psi,wrk,c,step) private(tid,i,j,k,jk,dl,d,du,b,ind,r) default(none) schedule(runtime)
  for (jk = 0; jk < nphiz; jk++) {  /* for each (phi,z) */
    j = jk / nz;
    k = jk % nz;
    tid = omp_get_thread_num();
    dl = &wrk[nr * (4 * tid + 0)];
    d  = &wrk[nr * (4 * tid + 1)];
    du = &wrk[nr * (4 * tid + 2)];
    b  = &wrk[nr * (4 * tid + 3)];

#define EPS 1.0E-3

    for(i = 0; i < nr; i++) {
      r = step * (double) i;
      if(i == nr-1) dl[nr-1] = 2.0;          /* Neumann outer boundary */
      else if(i > 0) dl[i] = 1.0 - step / (2.0 * r);
      if(i == 0) du[0] = step / (r + EPS);  /* Neumann inner boundary */
      else if(i < nr-1) du[i] = 1.0 + step / (2.0 * r);
      d[i] = ca;
    }

    /* create right-hand side vector */
    for(i = 1; i < nr - 1; i++) {
      ind = i * nphiz + j * nz + k;
      r = step * (double) i;
      /* (C - dr^2 laplace - dr/(2r) d/dr) psi(r,t) */
      b[i] = c * psi[ind] - (psi[ind + nphiz] - 2.0 * psi[ind] + psi[ind - nphiz]) - (step / (2.0 * r)) * (psi[ind + nphiz] - psi[ind - nphiz]);
    }
 
    /* neumann inner boundary */
    i = 0; ind = i * nphiz + j * nz + k;
    r = 0.0;
    b[i] = c * psi[ind] - (step * psi[ind + nphiz] / (r + EPS) - 2.0 * psi[ind]);

    /* neumann outer boundary */
    i = nr - 1; ind = i * nphiz + j * nz + k;
    b[i] = c * psi[ind] - (2.0 * psi[ind - nphiz] - 2.0 * psi[ind]);
    
    grid_solve_tridiagonal_system(nr, dl, d, du, b, dl);

    for(i = 0; i < nr; i++)
      psi[i * nphiz + j * nz + k] = dl[i];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along phi using Neumann boundary condition (Crank-Nicolson). 3D cylindrical.
 * Note: cylindrical BC for phi.
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_cn_cyl_phi(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step, r;
  long tid, ind;
  long nr = gwf->grid->nx;
  long nphi = gwf->grid->ny, nz = gwf->grid->nz, nphiz = nphi * nz;
  long ik, k, i, j;
  double complex *dl, *d, *du, *b, *x, *zwrk, alpha, beta;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt hbar) psi(t)  <=>  A x = b
   * (C + dphi^2 laplace) psi(t+dt) = (C - dphi^2 laplace) psi(t)
   * where C = 4 i m dx^2 / (hbar dt)
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  alpha = beta = 1.0; /* periodic */
  
#pragma omp parallel for firstprivate(nphiz,nz,nr,nphi,psi,wrk,c,step,alpha,beta) private(tid,ik,k,dl,d,du,b,ind,r,x,zwrk,ca,i,j) default(none) schedule(runtime)
  for(ik = 0; ik < nr * nz; ik++) { /* for each (r,z) */
    tid = omp_get_thread_num();
    dl   = &wrk[nphi * (6 * tid + 0)];
    d    = &wrk[nphi * (6 * tid + 1)];
    du   = &wrk[nphi * (6 * tid + 2)];
    b    = &wrk[nphi * (6 * tid + 3)];
    x    = &wrk[nphi * (6 * tid + 4)];
    zwrk = &wrk[nphi * (6 * tid + 5)];
    i = ik / nz;
    k = ik % nz;
    r = step * (double) i;
    ca = r * r * c;

    /* Construct the left hand side matrix */
    /* (r^2 * C + dphi^2 laplace) psi(t+dt) */
    for(j = 0; j < nphi; j++) {
      dl[j] = 1.0;
      du[j] = 1.0;
      d[j] = ca - 2.0;
    }

    /* create right-hand side vector */
    for(j = 1; j < nphi - 1; j++) {
      ind = i * nphiz + j * nz + k;
      /* (r^2 * C - dphi^2 laplace) psi(t) */
      b[j] = ca * psi[ind] - (psi[ind + nz] - 2.0 * psi[ind] + psi[ind - nz]);
    }
    /* periodic boundaries */
    j = 0;
    ind = i * nphiz + j * nz + k;
    b[j] = ca * psi[ind] - (psi[ind + nz] - 2.0 * psi[ind] + alpha * psi[ind + nz*(nphi-1)]);
    j = nphi - 1;
    ind = i * nphiz + j * nz + k;
    b[j] = ca * psi[ind] - (psi[ind - nz] - 2.0 * psi[ind] + beta * psi[ind - nz*(nphi-1)]);

    grid_solve_tridiagonal_system_cyclic(nphi, du, d, dl, b, alpha, beta, x, zwrk);

    for(j = 0; j < nphi; j++)
      psi[i * nphiz + j * nz + k] = x[j];
  }
}

/*
 * Auxiliary routine for propagating kinetic energy along z using Neumann boundary condition (Crank-Nicolson). 3D cyl.
 *
 * gwf       = wavefunction to be propagated (wf3d *).
 * time      = time step (double complex).
 * workspace = additional storage space needed (cgrid3d *).
 *
 * No return value.
 *
 */

EXPORT void grid3d_wf_propagate_kinetic_cn_cyl_z(wf3d *gwf, double complex time, cgrid3d *workspace) {

  double complex c, ca, *psi = gwf->grid->value, *wrk = workspace->value;
  double step = gwf->grid->step;
  long tid, ind;
  long i, ij, nr = gwf->grid->nx;
  long j, nphi = gwf->grid->ny;
  long k, nz = gwf->grid->nz, nphiz = nphi * nz;
  double complex *du, *dl, *d, *b;
  
  /*
   * (1 + .5 i T dt / hbar) psi(t+dt) = (1 - .5 i T dt / hbar) psi(t)  <=>  A x = b
   * (C + dx^2 laplace) psi(t+dt) = (C - dx^2 laplace) psi(t)
   * where C = i (4 m dx^2 / (hbar dt))
   */
  c = I * 4.0 * gwf->mass * step * step / (HBAR * time);

  /* (C + dx^2 laplace) psi(t+dt) */
  ca = c - 2.0;

#pragma omp parallel for firstprivate(ca,nz,nr,nphi,nphiz,psi,wrk,c,step) private(tid,i,j,k,ij,dl,d,du,b,ind) default(none) schedule(runtime)
  for (ij = 0; ij < nr * nphi; ij++) {  /* for each (r,phi) */
    i = ij / nphi;
    j = ij % nphi;
    tid = omp_get_thread_num();
    dl = &wrk[nz * (4 * tid + 0)];
    d  = &wrk[nz * (4 * tid + 1)];
    du = &wrk[nz * (4 * tid + 2)];
    b  = &wrk[nz * (4 * tid + 3)];

    for(k = 0; k < nz; k++) {
      if(k == nz-1) dl[nz-1] = 2.0;              /* Neumann outer boundary */
      else if(k > 0) dl[k] = 1.0;
      if(k == 0) du[0] = 2.0;                    /* Neumann inner boundary */
      else if(k < nz-1) du[k] = 1.0;
      d[k] = ca;
    }

    /* create right-hand side vector */
    for(k = 1; k < nz - 1; k++) {
      ind = i * nphiz + j * nz + k;
      /* (C - dx^2 laplace) psi(t) */
      b[k] = c * psi[ind] - (psi[ind + 1] - 2.0 * psi[ind] + psi[ind - 1]);
    }
 
    /* neumann boundaries */
    k = 0; ind = i * nphiz + j * nz + k;
    b[k] = c * psi[ind] - (2.0 * psi[ind + 1] - 2.0 * psi[ind]);
    k = nz - 1; ind = i * nphiz + j * nz + k;
    b[k] = c * psi[ind] - (2.0 * psi[ind - 1] - 2.0 * psi[ind]);
    
    grid_solve_tridiagonal_system(nz, dl, d, du, b, dl);

    for(k = 0; k < nz; k++)
      psi[i * nphiz + j * nz + k] = dl[k];
  }
}

/*
 * Map a given function on a wavefunction. 3D cylindrical.
 *
 * gwf  = wavefunction where function will be mapped to (wf3d *).
 * func = function providing the mapping (double complex (*)(void *, double)).
 * farg = optional argument for passing parameters to func (void *).
 *
 * No return value.
 *
 */

EXPORT inline void grid3d_wf_map_cyl(wf3d *gwf, double complex (*func)(void *arg, double z, double r, double phi), void *farg) { 

  cgrid3d_map_cyl(gwf->grid, func, farg); 
}

/*
 * Calculate the norm of the given wavefunction. 3D cyl.
 *
 * gwf = wavefunction for the calculation (wf3d *).
 *
 * Returns the norm (double).
 *
 */

EXPORT inline double grid3d_wf_norm_cyl(const wf3d *gwf) { 

  return cgrid3d_integral_of_square_cyl(gwf->grid); 
}

/*
 * Normalize wavefunction (to the value given in gwf->norm). 3D cyl.
 *
 * gwf = wavefunction to be normalized (wf3d *).
 *
 * Returns the normalization constant (double).
 *
 */

EXPORT inline double grid3d_wf_normalize_cyl(wf3d *gwf) { 

  double norm = grid3d_wf_norm_cyl(gwf);

  cgrid3d_multiply(gwf->grid, sqrt(gwf->norm / norm));
  return norm; 
}

/*
 * Calculate overlap between two wavefunctions. 3D cyl.
 *
 * gwfa = 1st wavefunction (wf3d *).
 * gwfb = 2nd wavefunction (wf3d *).
 *
 * Returns the overlap (double complex).
 *
 */

EXPORT inline double complex grid3d_wf_overlap_cyl(const wf3d *gwfa, const wf3d *gwfb) { 

  return cgrid3d_integral_of_conjugate_product_cyl(gwfa->grid, gwfb->grid); 
}

/* 
 * Calculate absorbing boundary potential. 3D cyl.
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
 * c           = Multiply the absorption term by this constant (double complex) (usually 1.0).
 *
 * The calculated potential should be added to the potential that will be
 * propagated. Note that it is important that this is also included in the
 * possible predict-correct cycle as that improves the numericla stability.
 *
 */

EXPORT void grid3d_wf_absorb_cyl(cgrid3d *potential, rgrid3d *density, double rho0, double (*region)(void *, double, double, double), rgrid3d *workspace, double complex c) {

  rgrid3d_copy(workspace, density);
  rgrid3d_add(workspace, -rho0);
  rgrid3d_product_func_cyl(workspace, region, NULL);
  rgrid3d_multiply(workspace, -c);
  grid3d_add_real_to_complex_im(potential, workspace);
}
