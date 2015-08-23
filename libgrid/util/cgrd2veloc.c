/*
 * Convert complex grid format to netcdf probability flux field (3D only).
 * Divide this by density to get liquid velocity.
 *
 * Note: To use this with 3D cylindical, write a
 * Cartesian 3D grid instead.
 *
 * Usage: cgrd2veloc gridfile cdffile
 *
 */

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <netcdf.h>
#include <grid/grid.h>

int main(int argc, char **argv) {

  char *cdffile = argv[2], *gridfile = argv[1];
  FILE *fp;
  long nx, ny, nz, i;
  int ncid, varidx, varidy, varidz, varidxyz, x_dimid, y_dimid, z_dimid, dimids[3], retval;
  double step;
  wf3d *wf;
  rgrid3d *flux_x, *flux_y, *flux_z, *flux_mag;
  
  if(argc != 3) {
    fprintf(stderr, "Usage: cgrd2veloc cgridfile cdffile\n");
    exit(1);
  }

  /* Read in the complex grid */
  if(!(fp = fopen(gridfile, "r"))) {
    fprintf(stderr, "Can't open grid file %s.\n", gridfile);
    exit(1);
  }
  grid3d_read_peek(fp, &nx, &ny, &nz, &step);
  fprintf(stderr, "nx = %ld, ny = %ld, nz = %ld, step = %le\n", nx, ny, nz, step);
  wf = grid3d_wf_alloc(nx, ny, nz, step, 1.0, WF3D_NEUMANN_BOUNDARY, WF3D_2ND_ORDER_PROPAGATOR);  
  cgrid3d_read(wf->grid, fp);
  fclose(fp);

  /* Calculate the velocity field */
  flux_x = rgrid3d_alloc(nx, ny, nz, step, RGRID3D_NEUMANN_BOUNDARY, NULL);
  flux_y = rgrid3d_alloc(nx, ny, nz, step, RGRID3D_NEUMANN_BOUNDARY, NULL);
  flux_z = rgrid3d_alloc(nx, ny, nz, step, RGRID3D_NEUMANN_BOUNDARY, NULL);
  flux_mag = rgrid3d_alloc(nx, ny, nz, step, RGRID3D_NEUMANN_BOUNDARY, NULL);
  grid3d_wf_probability_flux(wf, flux_x, flux_y, flux_z);
  for (i = 0; i < nx*ny*nz; i++)
    flux_mag->value[i] = flux_x->value[i] * flux_x->value[i] + flux_y->value[i] * flux_y->value[i] + flux_z->value[i] * flux_z->value[i];

  /* Write CDF file */
  if((retval = nc_create(cdffile, NC_CLOBBER, &ncid))) {
    puts(nc_strerror(retval));
    fprintf(stderr, "Error in nc_open().\n");
    exit(1);
  }
  nc_def_dim(ncid, "x", nx, &x_dimid);
  nc_def_dim(ncid, "y", ny, &y_dimid);
  nc_def_dim(ncid, "z", nz, &z_dimid);
  dimids[0] = x_dimid;
  dimids[1] = y_dimid;
  dimids[2] = z_dimid;
  nc_def_var(ncid, "fx", NC_DOUBLE, 3, dimids, &varidx);
  nc_def_var(ncid, "fy", NC_DOUBLE, 3, dimids, &varidy);
  nc_def_var(ncid, "fz", NC_DOUBLE, 3, dimids, &varidz);
  nc_def_var(ncid, "fmag", NC_DOUBLE, 3, dimids, &varidxyz);
  nc_enddef(ncid);
  nc_put_var_double(ncid, varidx, flux_x->value);
  nc_put_var_double(ncid, varidy, flux_y->value);
  nc_put_var_double(ncid, varidz, flux_z->value);
  nc_put_var_double(ncid, varidxyz, flux_mag->value);
  nc_close(ncid);
  return 0;
}

