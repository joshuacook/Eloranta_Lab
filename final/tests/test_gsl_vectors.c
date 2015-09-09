#include <stdio.h>
#include <gsl/gsl_vector.h>

  int i;
  gsl_vector * v = gsl_vector_alloc (3);

static char * test_vector_allocation()
{
   w
}

static char * all_tests()
{
  mu_run_test(test_vector_allocation);
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
