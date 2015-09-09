#include "unity.h"
#include "../src/test_gsl_eigs.c"

#define DIMENSION_OF(a) (sizeof(a)/sizeof(a[0]))

void test_pass(void)
{
  TEST_ASSERT_EQUAL(3,3);
}

void test_eigenvalues_of_4th_order_Hilbert_Matrix(void)
{
  double expected_eigenvalues[] = {1.50021,1.69141e-01,6.73827e-03,9.67023e-05};
  int i;
  
  for (i=0;i < DIMENSION_OF(expected_eigenvalues); i++)
  {
    printf("\nexpected: %f actual: %f\n",expected_eigenvalues[i],gsl_eval_4th_order_hilbert_matrix()[i]);
    TEST_ASSERT_EQUAL(expected_eigenvalues[i],gsl_eval_4th_order_hilbert_matrix()[i]);
  }
}

int main(void) 
{
  UNITY_BEGIN();
  RUN_TEST(test_pass);
  RUN_TEST(test_eigenvalues_of_4th_order_Hilbert_Matrix);
  return UNITY_END();
}

