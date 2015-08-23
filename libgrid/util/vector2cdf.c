/*
 * Convert flux vector field to netcdf (3D only).
 *
 * Usage: vector2cdf density flux_x flux_y flux_z cdffile
 *
 * In visit, convert the fx, fy and fz fields to a vector field by:
 * 1) Controls -> Expressions
 * 2) Enter name: vecfield
 * 3) Type: Vector Mesh Variable
 * 4) Standard Editor: {fx,fy,fz}
 * 5) Apply.
 * 
 * Use vecfield variable for plotting the vectori field.
 * Note that density is also included and can be overlaid
 * with the vector field in visit.
 *
 * NOTE: This reverses the X and Z axes in order to be compatible
 *       with the way VisIT reads netcdf files.
 *
 */

#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <complex.h>
#include <stdio.h>
#include <netcdf.h>

double *array[4], *array_rev = NULL, *x, *y, *z;
long nx, ny, nz;

void reverse_xz(double *in, double *out) {

  long i, j, k;

  for (i = 0; i < nx; i++)
    for (j = 0; j < ny; j++)
      for (k = 0; k < nz; k++)
	out[k * nx * ny + j * nx + i] = in[i * nz * ny + j * nz + k];
}

void read_grid(int what, char *file) {

  FILE *fp;
  double step;
  long i;
  static int been_here = 0;

  if(!(fp = fopen(file, "r"))) {
    fprintf(stderr, "Can't open grid file %s.\n", file);
    exit(1);
  }

  nx = ny = nz = 1;
  fread(&nx, sizeof(long), 1, fp);
  fread(&ny, sizeof(long), 1, fp);
  fread(&nz, sizeof(long), 1, fp);
  fread(&step, sizeof(double), 1, fp);

  if(!been_here) {
    if(!(x = (double *) malloc(sizeof(double) * nx))) {
      fprintf(stderr, "Out of memory.\n");
      exit(1);
    }
    if(!(y = (double *) malloc(sizeof(double) * ny))) {
      fprintf(stderr, "Out of memory.\n");
      exit(1);
    }
    if(!(z = (double *) malloc(sizeof(double) * nz))) {
      fprintf(stderr, "Out of memory.\n");
      exit(1);
    }
    // TODO: the grid files do not store the grid center x0,y0,z0...
    // origin assumed to be at the center
    for (i = 0; i < nx; i++) x[i] = (i - nx/2) * step;
    for (i = 0; i < ny; i++) y[i] = (i - ny/2) * step;
    for (i = 0; i < nz; i++) z[i] = (i - nz/2) * step;
    been_here = 1;
  }

  fprintf(stderr, "File %s: nx = %ld, ny = %ld, nz = %ld, step = %le\n", file, nx, ny, nz, step);
  if(!(array[what] = (double *) malloc(sizeof(double) * nx * ny * nz))) {
    fprintf(stderr, "Out of memory.\n");
    exit(1);
  }
  fread(array[what], sizeof(double), nx * ny * nz, fp);
  fclose(fp);

  if(array_rev == NULL) {
    if(!(array_rev = (double *) malloc(sizeof(double) * nx * ny * nz))) {
      fprintf(stderr, "Out of memory.\n");
      exit(1);
    }
  }  
  reverse_xz(array[what], array_rev);
  bcopy(array_rev, array[what], sizeof(double) * nx * ny * nz);
}

int main(int argc, char **argv) {

  int ncid, varid1, varid2, varid3, varid4, varid5, varid6, varid7;
  int dimids[3], retval;

  if(argc != 6) {
    fprintf(stderr, "Usage: vector2cdf density flux_x flux_y flux_z cdffile\n");
    exit(1);
  }
  read_grid(0, argv[2]);  // x
  read_grid(1, argv[3]);  // y
  read_grid(2, argv[4]);  // z
  read_grid(3, argv[1]);  // density
  
  if((retval = nc_create(argv[5], NC_CLOBBER, &ncid))) {
    puts(nc_strerror(retval));
    fprintf(stderr, "Error in nc_open().\n");
    exit(1);
  }
  nc_def_dim(ncid, "z", nz, &dimids[0]);
  nc_def_var(ncid, "z", NC_DOUBLE, 1, &dimids[0], &varid5);
  nc_def_dim(ncid, "y", ny, &dimids[1]);
  nc_def_var(ncid, "y", NC_DOUBLE, 1, &dimids[1], &varid6);
  nc_def_dim(ncid, "x", nx, &dimids[2]);
  nc_def_var(ncid, "x", NC_DOUBLE, 1, &dimids[2], &varid7);
  nc_def_var(ncid, "fx", NC_DOUBLE, 3, dimids, &varid1);
  nc_def_var(ncid, "fy", NC_DOUBLE, 3, dimids, &varid2);
  nc_def_var(ncid, "fz", NC_DOUBLE, 3, dimids, &varid3);
  nc_def_var(ncid, "rho", NC_DOUBLE, 3, dimids, &varid4);
  nc_enddef(ncid);
  nc_put_var_double(ncid, varid1, array[0]);
  nc_put_var_double(ncid, varid2, array[1]);
  nc_put_var_double(ncid, varid3, array[2]);
  nc_put_var_double(ncid, varid4, array[3]);
  nc_put_var_double(ncid, varid5, z);
  nc_put_var_double(ncid, varid6, y);
  nc_put_var_double(ncid, varid7, x);
  nc_close(ncid);
  return 0;
}

