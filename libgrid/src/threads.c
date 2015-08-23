/*
 * Common functions.
 *
 */

#include "grid.h"
#include "private.h"

static int grid_nthreads = 0;

#ifndef _OPENMP
int omp_get_max_threads() { return 1;}
int omp_get_thread_num() {return 1;}
void omp_set_num_threads(int threads) {}
int omp_get_num_threads() { return 1;}
void ompc_set_dynamic(int asd) {}
void omp_set_dynamic(int asd) {}
#endif

#ifdef __GNUC__
#include <math.h>
#include <fenv.h>
extern int fegetexcept();
extern int feenableexcept(int);
#endif

EXPORT void grid_threads_init(int threads) {

  if (threads <= 0) {
    fprintf(stderr, "libgrid: Error in grid_threads_init(). Number of threads <= 0.\n");
    abort();
  }

  if (grid_nthreads && grid_nthreads != threads) {
    fprintf(stderr, "libgrid: Error in grid_threads_init(). This function was called twice with different number of threads.\n");
    abort();	
  }

#ifdef __GNUC__
  /* core dump on floating point exceptions */
  //feenableexcept(fegetexcept()|FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW|FE_UNDERFLOW);
  feenableexcept(fegetexcept()|FE_DIVBYZERO|FE_INVALID|FE_OVERFLOW);   // Do not attempt to catch underflow - this is just zero to us...
#endif
  
  /* Without openmp, use just one thread for everything */
#ifdef _OPENMP
  grid_nthreads = threads;
#else
  grid_nthreads = 1;
#endif
  // TODO: This was left here - is it still relevant (disabled for now)
  //  omp_set_dynamic(0);
  omp_set_num_threads(threads);
  fftw_init_threads();
  fprintf(stderr, "libgrid: Initialized with %d threads.\n", threads);
  fprintf(stderr, "libgrid: omp_num_threads = %d, omp_max_num_threads = %d.\n",
	  omp_get_num_threads(), omp_get_max_threads());
}

EXPORT int grid_threads() {

  return grid_nthreads;
}

