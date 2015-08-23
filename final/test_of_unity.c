#include "unity.h"

void test_one(void)
{
  TEST_ASSERT_EQUAL(1,1);
}

int main(void) 
{
  UNITY_BEGIN();
  RUN_TEST(test_one);
  return UNITY_END();
}

