#ifndef __GRID_STRUCTS__
#define __GRID_STRUCTS__

typedef struct cgrid1d_struct {
  double complex *value;
  long nx;
  double step;
  double complex (*value_outside)(const struct cgrid1d_struct *grid, long i);
  void *outside_params_ptr;
  double complex default_outside_params;
  fftw_plan plan, iplan;
} cgrid1d;

typedef struct rgrid1d_struct {
  double *value;
  cgrid1d *cint;  /* complex structure interface after FFT */
  long nx;
  double step;
  double (*value_outside)(const struct rgrid1d_struct *grid, long i);
  void *outside_params_ptr;
  double default_outside_params;
  fftw_plan plan, iplan;
} rgrid1d;

typedef struct cgrid2d_struct {
  double complex *value;
  long nx, ny;
  double step;
  double complex (*value_outside)(const struct cgrid2d_struct *grid, long i, long j);
  void *outside_params_ptr;
  double complex default_outside_params;
  fftw_plan plan, iplan;
} cgrid2d;

typedef struct rgrid2d_struct {
  double *value;
  cgrid2d *cint;  /* complex structure interface after FFT */
  long nx, ny;
  double step;
  double (*value_outside)(const struct rgrid2d_struct *grid, long i, long j);
  void *outside_params_ptr;
  double default_outside_params;
  fftw_plan plan, iplan;
  fftw_plan plan_cyl, iplan_cyl;
  gsl_dht *gh;     /* for hankel (nu = 0) */
  double gh_norm;  /* for hankel (nu = 0) */
} rgrid2d;

typedef struct cgrid3d_struct {
  double complex *value;
  long nx, ny, nz;
  double step;
  double x0, y0, z0;
  double kx0, ky0, kz0;
  double complex (*value_outside)(const struct cgrid3d_struct *grid, long i, long j, long k);
  void *outside_params_ptr;
  double complex default_outside_params;
  fftw_plan plan, iplan, implan, iimplan;
  double fft_norm;                           /* TODO: should probably do this for 1d & 2d too */
} cgrid3d;

typedef struct rgrid3d_struct {
  double *value;
  cgrid3d *cint;
  long nx, ny, nz;
  double step;
  double x0, y0, z0;
  double kx0, ky0, kz0;
  double (*value_outside)(const struct rgrid3d_struct *grid, long i, long j, long k);
  void *outside_params_ptr;
  double default_outside_params;
  fftw_plan plan, iplan;
  double fft_norm;
} rgrid3d;

typedef struct wf1d_struct {
  cgrid1d *grid;
  double mass;
  double norm;
  int boundary;
  int propagator;
} wf1d;

typedef struct wf2d_struct {
  cgrid2d *grid;
  double mass;
  double norm;
  int boundary;
  int propagator;
} wf2d;

typedef struct wf3d_struct {
  cgrid3d *grid;
  double mass;
  double norm;
  int boundary;
  int propagator;
} wf3d;

typedef struct grid_fcmpl {
  float re;
  float im;
} grid_fcmpl;

typedef struct grid_timer_struct {
  struct timeval zero_time;
  clock_t zero_clock;
} grid_timer;

typedef struct rotation_struct{
	rgrid3d *rgrid;
	cgrid3d *cgrid;
	double sinth ;
	double costh ;
} rotation;


#endif /* __GRID_STRUCTS__ */
