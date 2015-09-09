#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

double fourth_order_hilbert[] =
{
  1.0  , 1/2.0, 1/3.0, 1/4.0,
  1/2.0, 1/3.0, 1/4.0, 1/5.0,
  1/3.0, 1/4.0, 1/5.0, 1/6.0,
  1/4.0, 1/5.0, 1/6.0, 1/7.0
}

gsl_matrix_view fourth_order_view = gsl_matrix_view_array(
    fourth_order_hilbert, 4, 4);

