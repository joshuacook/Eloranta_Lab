#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "minunit.h"
#include "fixtures/test_gsl_eigenproblems"

#define DIMENSION_OF(a) (sizeof(a)/sizeof(a[0]))

gsl_vector * symmetric_eigenvalues(gsl_matrix_view * matrix, int m, int n);
gsl_matrix * symmetric_eigenvectors(gsl_matrix_view * matrix, int m, int n);

int tests_run = 0;

static char * test_solve_eigenproblems()
{
  gsl_vector * eigenvalues  = symmetric_eigenvalues(fourth_order_view,4,4);
  gsl_matrix * eigenvectors = symmetric_eigenvectors(fourth_order_view,4,4)

  for (int i = 0; i < DIMENSION_OF(expected_eigenvalues); i++)
  {
     float expected = expected_eigenvalues[i];
     float actual = round_to_5_places(gsl_vector_get(actual_eigenvalues, i));

     printf("\nexpected: %f actual: %f\n",
       expected,actual);
     mu_assert("error, eigenvalues do not match",
         expected == actual);
  }
  return 0;
}

static char * all_tests()
{
  mu_run_test(test_solve_eigenproblems);
  return 0;
}

int main(int argc, char **argv)
{
  char * result = all_tests();
  if (result != 0)
  {
    printf("%s\n", result);
  }
  else
  {
    printf("All Tests Passed\n");
  }
  printf("Tests run: %d\n", tests_run);

  return result != 0;
}

gsl_vector * symmetric_eigenvalues(gsl_matrix_view * matrix, int m, int n)
{
  gsl_vector * my_evals = gsl_vector_alloc (m);
  gsl_matrix * my_evecs = gsl_matrix_alloc (m, n);

  gsl_eigen_symmv_workspace * my_workspace = gsl_eigen_symmv_alloc (m);

  gsl_eigen_symmv (&matrix.matrix, my_evals, my_evecs, my_workspace);

  gsl_eigen_symmv_free (my_workspace);

  gsl_eigen_symmv_sort (my_evals, my_evecs, GSL_EIGEN_SORT_ABS_DESC);

  return my_evals;
}

gsl_matrix * symmetric_eigenvectors(gsl_matrix_view * matrix, int m, int n)
{
  gsl_vector * my_evals = gsl_vector_alloc (m);
  gsl_matrix * my_evecs = gsl_matrix_alloc (m, n);

  gsl_eigen_symmv_workspace * my_workspace = gsl_eigen_symmv_alloc (m);

  gsl_eigen_symmv (&matrix.matrix, my_evals, my_evecs, my_workspace);

  gsl_eigen_symmv_free (my_workspace);

  gsl_eigen_symmv_sort (my_evals, my_evecs, GSL_EIGEN_SORT_ABS_DESC);

  return my_evecs;
}

// Method for displaying eigenvalues and eigenvectors
//     int i;
//
//     for (i = 0; i < 4; i++)
//     {
//       double eval_i = gsl_vector_get (my_evals, i);
//       gsl_vector_view evec_i = gsl_matrix_column (my_evecs, i);
//
//       printf ("eigenvalue = %g\n", eval_i);
//       printf ("eigenvector = \n");
//       gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
//     }
//   }



