/*
 * Interface to Discrete Hankel and Fourier transform.
 *
 * Uses FFTW along the z-axis and DHT along r.
 *
 * This is used in computing Fourier transformation in cylindrical
 * coordinates.
 *
 * Notes: 
 * - index as grid(z, r), where r runs closest in memory.
 *   So "x = z" and "y = r" and "NZ = NX" and "NR = NY".
 * - Along z-axis the data is stored in halfcomplex format.
 *
 * TODO: This seems not to work well with opencc and icc!?
 *
 */

#include "grid.h"

static double **temp;
static long temp_len = 0;
static double *temp2 = NULL;
static long temp2_len = 0;

/* 
 * Clean up the artificial tail caused by the Hankel transformation.
 *
 * grid       = grid to be operated on (rgrid2d *, input/output).
 * hankel_pad = number of points to zero (long, input).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_fft_cylindrical_cleanup(rgrid2d *grid, long hankel_pad) {

  long nx = grid->nx, ny = grid->ny;
  long i, j, l, l2;
  double *val = grid->value;             /* this is due to imperfect */
                                         /* boundaries for Hankel */

  for(i = 0; i < nx; i++) { /* z */
    l = i * ny + ny - hankel_pad;
    l2 = i * ny;
    for (j = ny - hankel_pad; j < ny; j++) /* r */
      val[l2 + j] = val[l];
  }
}

/* In X space, regular to Bessel zeros */
static void interpol1(double *src, double *dst, long n, double step, gsl_dht *gh) {

  long i, idxl, idxu;
  double x, xl, a;

  for (i = 0; i < n; i++) {
    x = gsl_dht_x_sample(gh, i);
    idxl = (long) (x / step);
    idxu = idxl + 1;
    if(idxu == n) {
      a = 0.0;
      idxu = idxl;
    } else {
      xl = step * (double) idxl;
      a = (x - xl) / step;
    }
    dst[i] = (1.0 - a) * src[idxl] + a * src[idxu];
  }
}

/* In X space, Bessel zeros to regular */
static void interpol2(double *src, double *dst, long n, double step, gsl_dht *gh) {

  long i, j, idxl;
  double x, xu, xl, a;

  for (i = 0; i < n; i++) {
    x = step * (double) i;
    for (j = 0; j < n; j++)
      if(x <= (xu = gsl_dht_x_sample(gh, j))) break;
    if(j >= n-1) { /* can't go far enough, default to last point available */
      a = 1.0;
      j = n-1;
      idxl = 0;
    } if(j == 0) { /* the first point, use that value */
      a = 1.0;
      j = 0;
      idxl = 0;
    } else {
      xl = gsl_dht_x_sample(gh, j-1);
      a = (x - xl) / (xu - xl);
      idxl = j-1;
    }
    dst[i] = a * src[j] + (1.0 - a) * src[idxl];
  }
}

/* allocate transform structure (only used internally) */
EXPORT void rgrid2d_fft_cylindrical_alloc(rgrid2d *grid) {

  int nx = grid->nx, ny = grid->ny;
  double ival = grid->step * ny;
  long i, s = nx * ny * sizeof(double), s2 = ny * sizeof(double);
  fftw_r2r_kind r2hc = FFTW_R2HC, hc2r = FFTW_HC2R;

  if(s > temp2_len) {
    if(temp2_len) free(temp2);
    if(!(temp2 = fftw_malloc(s))) {
      fprintf(stderr, "libgrid: Out of memory in rgrid2d_fft().\n");
      return;
    }
    temp2_len = s;
  }

  if(s2 > temp_len) {
    if(temp_len)
      for (i = 0; i < omp_get_max_threads(); i++)
	free(temp[i]);
    else if(!(temp = (double **) malloc(sizeof(double *) * omp_get_max_threads()))) {
      fprintf(stderr, "libgrid: Out of memory in rgrid2d_fft_cylindrical_alloc().\n");
      return;
    }

    for(i = 0; i < omp_get_max_threads(); i++)
      if(!(temp[i] = fftw_malloc(s2))) {
	fprintf(stderr, "libgrid: Out of memory in rgrid2d_fft_cylindrical_alloc()().\n");
	return;
      }
    temp_len = s2;
  }

  memcpy(temp2, grid->value, s);
 
  fftw_plan_with_nthreads(grid_threads());

  /* TODO: is in-place transformation OK speedwise? */
  grid->plan_cyl =
    fftw_plan_many_r2r(1, &nx, ny, grid->value, NULL, ny, 1, grid->value, NULL, ny, 1, &r2hc, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);

  grid->iplan_cyl =
    fftw_plan_many_r2r(1, &nx, ny, grid->value, NULL, ny, 1, grid->value, NULL, ny, 1, &hc2r, GRID_FFTW_PLAN | FFTW_DESTROY_INPUT);

  memcpy(grid->value, temp2, s);

  grid->gh = gsl_dht_new(ny, 0.0, grid->step * (double) ny);
  grid->gh_norm = 2.0 * pow(gsl_sf_bessel_zero_J0(ny+1), 3.0) / (ival * ival * ival * ival * ny);
}

EXPORT void rgrid2d_fft_cylindrical_free(rgrid2d *grid) {

  if(grid->plan_cyl) fftw_destroy_plan(grid->plan_cyl);
  if(grid->iplan_cyl) fftw_destroy_plan(grid->iplan_cyl);
  if(grid->gh) gsl_dht_free(grid->gh);
}

/*
 * FFT in cylindrical coordinates.
 *
 * grid = grid to be transformed (rgrid2d *; input/output).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_fft_cylindrical(rgrid2d *grid) {

  long i, nx = grid->nx, ny = grid->ny, idx;
  int thr;
  double *val = grid->value, step = grid->step;
  gsl_dht *gh;

  if(grid->plan_cyl == NULL) rgrid2d_fft_cylindrical_alloc(grid);
  gh = grid->gh;

  /* execute Hankel transformation for each z-row */
#pragma omp parallel for firstprivate(nx,ny,step,val,temp,gh) private(i,idx,thr) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {
    thr = omp_get_thread_num();
    idx = i * ny;
    interpol1(val + idx, temp[thr], ny, step, gh);/* regular to zeros of Bessel funct */
    gsl_dht_apply(gh, temp[thr], val + idx);
    /* note: the spacing in k space is irregular ! */
  }
#pragma omp barrier
  fftw_execute(grid->plan_cyl);
}

/*
 * Inverse FFT in cylindrical coordinates.
 *
 * grid = grid to be transformed (rgrid2d *; input/output).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_inverse_fft_cylindrical(rgrid2d *grid) {

  long i, nx = grid->nx, ny = grid->ny, idx;
  int thr;
  double *val = grid->value, step = grid->step;
  gsl_dht *gh;

  if(grid->iplan_cyl == NULL) rgrid2d_fft_cylindrical_alloc(grid);
  gh = grid->gh;

  /* execute Hankel transformation for each z-row */
#pragma omp parallel for firstprivate(nx,ny,step,val,temp,gh) private(i,idx,thr) default(none) schedule(runtime)
  for (i = 0; i < nx; i++) {
    thr = omp_get_thread_num();
    idx = i * ny;
    gsl_dht_apply(gh, val + idx, temp[thr]);
    interpol2(temp[thr], val + idx, ny, step, gh);   /* back to regular */
  }
#pragma omp barrier
  fftw_execute(grid->iplan_cyl);
}

/*
 * Convolution in cylindrical coordinates (in Fourier space).
 * 
 * dest  = destination grid (rgrid2d *; output).
 * grid1 = source grid1 (rgrid2d *; input/output).
 * grid2 = source grid2 (rgrid2d *; input/output).
 *
 * No return value.
 *
 */

EXPORT void rgrid2d_fft_cylindrical_convolute(rgrid2d *dest, rgrid2d *grid1, rgrid2d *grid2) {

  long i, j, nx = grid1->nx, ny = grid1->ny, ii, ir, nx2 = nx/2, nxy = nx * ny;
  int even = nx & 1;
  double tot_nrm;
  double *g1 = grid1->value, *g2 = grid2->value, *dst = dest->value;
  double *orig_dest = NULL;

  tot_nrm = grid1->gh_norm * grid1->step / (double) nx;

  if(dest == grid1 || dest == grid2) {
    orig_dest = dest->value;
    dst = temp2;
  }

  /* multiply */
#pragma omp parallel for firstprivate(dst,g1,g2,ny) private(j) default(none) schedule(runtime)
  for (j = 0; j < ny; j++) /* DC */
    dst[j] = g1[j] * g2[j];    

  if(nx > 1) {
    if(even) { /* only real at N/2 */
      double asd = nx2 * ny;
#pragma omp parallel for firstprivate(dst,g1,g2,ny,asd) private(j,ir) default(none) schedule(runtime)
      for (j = 0; j < ny; j++) {
	ir = asd + j;
	dst[ir] = g1[ir] * g2[ir];     
      }
    } else { /* odd, both real and imag. at N/2 */
      double asd = nx2 * ny;
#pragma omp parallel for firstprivate(dst,g1,g2,ny,asd) private(j,ii,ir) default(none) schedule(runtime)
      for (j = 0; j < ny; j++) {
	ir = asd + j;
	ii = ir + ny;
	dst[ir] = g1[ir] * g2[ir] - g1[ii] * g2[ii];
	dst[ii] = g1[ir] * g2[ii] + g1[ii] * g2[ir];
	dst[ir] *= -1.0; /* wrap */
	dst[ii] *= -1.0; /* wrap */
      }
    }
#pragma omp parallel for firstprivate(dst,g1,g2,nx2,ny,nx) private(i,j,ir,ii) default(none) schedule(runtime)
    for (i = 1; i < nx2; i++)  /* z */
      for (j = 0; j < ny; j++) { /* r */
	ir = i * ny + j;
	ii = (nx - i) * ny + j;
	dst[ir] = g1[ir] * g2[ir] - g1[ii] * g2[ii];
	dst[ii] = g1[ir] * g2[ii] + g1[ii] * g2[ir];
	if(i & 1) { /* odd */
	  dst[ir] *= -1.0; /* wrap */
	  dst[ii] *= -1.0; /* wrap */
	}
      }
  }

  /* normalize */
#pragma omp parallel for firstprivate(nxy,dst,tot_nrm) private(i) default(none) schedule(runtime)
  for (i = 0; i < nxy; i++)
    dst[i] *= tot_nrm;

  if(orig_dest)
    memcpy(orig_dest, dst, sizeof(double) * nxy);
}
