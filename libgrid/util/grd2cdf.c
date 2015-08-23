/*
 * Convert real grid format to netcdf (3D only).
 *
 * Usage: grd2cdf gridfile cdffile
 *
 * NOTE: This reverses the X and Z axes in order to be compatible
 *       with the way VisIT reads netcdf files.
 *
 */

#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <netcdf.h>

double *array, *array_rev, *x, *y, *z;
long nx, ny, nz;

void reverse_xy(double *in, double *out) {

  long i, j;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      out[j * nx + i] = in[i * ny + j];
}

void reverse_xz(double *in, double *out) {

  long i, j, k;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++)
	out[k * nx * ny + j * nx + i] = in[i * nz * ny + j * nz + k];
}

int main(int argc, char **argv) {

  char *cdffile = argv[3], *gridfile = argv[2];
  FILE *fp;
  long i;
  int ncid, varid0, varid1, varid2, varid3, dimids[3], retval, dim;
  double step;

  if(argc != 4) {
    fprintf(stderr, "Usage: grd2cdf {-1,-2,-3} gridfile cdffile\n");
    exit(1);
  }
  dim = 0;
  if(argv[1][1] == '1') dim = 1;
  if(argv[1][1] == '2') dim = 2;
  if(argv[1][1] == '3') dim = 3;
  if(!dim) {
    fprintf(stderr, "Use -1, -2, -3 to define the dimension.\n");
    exit(1);
  }
  if(!(fp = fopen(gridfile, "r"))) {
    fprintf(stderr, "Can't open grid file %s.\n", gridfile);
    exit(1);
  }

  nx = ny = nz = 1;
  fread(&nx, sizeof(long), 1, fp);
  if(dim > 1) fread(&ny, sizeof(long), 1, fp);
  if(dim > 2) fread(&nz, sizeof(long), 1, fp);
  fread(&step, sizeof(double), 1, fp);

  fprintf(stderr, "nx = %ld, ny = %ld, nz = %ld, step = %le\n", nx, ny, nz, step);
  if(!(array = (double *) malloc(sizeof(double) * nx * ny * nz))) {
    fprintf(stderr, "Out of memory.\n");
    exit(1);
  }
  if(!(array_rev = (double *) malloc(sizeof(double) * nx * ny * nz))) {
    fprintf(stderr, "Out of memory.\n");
    exit(1);
  }
  if(!(x = (double *) malloc(sizeof(double) * nx))) {
    fprintf(stderr, "Out of memory.\n");
    exit(1);
  }
  for (i = 0; i < nx; i++) x[i] = (i - nx/2) * step;

  if (dim > 1) {
    if(!(y = (double *) malloc(sizeof(double) * ny))) {
      fprintf(stderr, "Out of memory.\n");
      exit(1);
    }
    for (i = 0; i < ny; i++) y[i] = (i - ny/2) * step;
  }

  if (dim > 2) {
    if(!(z = (double *) malloc(sizeof(double) * nz))) {
      fprintf(stderr, "Out of memory.\n");
      exit(1);
    }
    for (i = 0; i < nz; i++) z[i] = (i - nz/2) * step;
  }
  
  fread(array, sizeof(double), nx * ny * nz, fp);
  fclose(fp);

  if(dim == 2) reverse_xy(array, array_rev);
  if(dim == 3) reverse_xz(array, array_rev);

  if((retval = nc_create(cdffile, NC_CLOBBER, &ncid))) {
    puts(nc_strerror(retval));
    fprintf(stderr, "Error in nc_open().\n");
    exit(1);
  }
  
  switch(dim) {
  case 3:
    nc_def_dim(ncid, "z", nz, &dimids[0]);
    nc_def_var(ncid, "z", NC_DOUBLE, 1, &dimids[0], &varid0);
    nc_def_dim(ncid, "y", ny, &dimids[1]);
    nc_def_var(ncid, "y", NC_DOUBLE, 1, &dimids[1], &varid1);
    nc_def_dim(ncid, "x", nx, &dimids[2]);
    nc_def_var(ncid, "x", NC_DOUBLE, 1, &dimids[2], &varid2);
    break;
  case 2:
    nc_def_dim(ncid, "y", ny, &dimids[0]);
    nc_def_var(ncid, "y", NC_DOUBLE, 1, &dimids[0], &varid0);
    nc_def_dim(ncid, "x", nx, &dimids[1]);
    nc_def_var(ncid, "x", NC_DOUBLE, 1, &dimids[1], &varid1);
    break;
  case 1:
    nc_def_dim(ncid, "x", nx, &dimids[0]);
    nc_def_var(ncid, "x", NC_DOUBLE, 1, &dimids[0], &varid0);
    break;
  }
  nc_def_var(ncid, "density", NC_DOUBLE, dim, dimids, &varid3);
  nc_enddef(ncid);
  nc_put_var_double(ncid, varid3, array_rev);
  switch(dim) {
  case 3:
    nc_put_var_double(ncid, varid2, x);
    nc_put_var_double(ncid, varid1, y);
    nc_put_var_double(ncid, varid0, z);
    break;
  case 2:
    nc_put_var_double(ncid, varid1, x);
    nc_put_var_double(ncid, varid0, y);
    break;
  case 1:
    nc_put_var_double(ncid, varid0, x);
    break;
  }
  nc_close(ncid);
  return 0;
}

