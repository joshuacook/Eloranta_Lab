
/* FFTW plan decision */
/* FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE */
#define GRID_FFTW_PLAN FFTW_MEASURE

#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <complex.h>

#include <time.h>
#include <sys/time.h>

#include <fftw3.h>
#include <gsl/gsl_dht.h>
#include <gsl/gsl_sf_bessel.h>

/* TODO: there should be grid/... */
#include "structs.h"
#include "proto.h"
#include "defs.h"

/* identifier for extracting prototypes */
#define EXPORT
