#include <stdio.h>
#include <gsl/gsl_blas.h>
#include "minunit.h"

int compare_vector();

int test_dgemm()

int
main (void)
{
  double
}

int compare_vector(gsl_vector * x, gsl_vector * x)
{

  for (int i = 0; i < DIMENSION_OF(expected_eigenvalues); i++)
  {
     float expected = expected_eigenvalues[i];
     float actual = round_to_5_places(gsl_vector_get(actual_eigenvalues, i));

     printf("\nexpected: %f actual: %f\n",
       expected,actual);
     mu_assert("error, eigenvalues do not match",
         expected == actual);
  }
}
