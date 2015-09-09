#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_block.h>
#include <math.h>
#include <my_minunit.h>

int tests_run = 0;
float round_to_6_places(float num);

static 
char * test_gsl_0_order_bessel_function_of_the_first_kind ()
{
  double x = 5.0;
  double y = round_to_6_places(gsl_sf_bessel_J0 (x));
  double expected_y = round_to_6_places(-1.775967713143382920e-01);
  printf("y: %f \t expected_y: %f\n",y,expected_y);
  mu_assert
  (
    "value of zero order Bessel function of the first kind",
    y == expected_y,
    "y: not equal to expected_y" 
  );
  return 0;
}

static
char * test_gsl_rectangular_complex_number_struct()
{
  double x = 2.43728;
  double y = 3.23412;

  gsl_complex test_rect_complex_number = gsl_complex_rect ( x, y ); 

  mu_assert
  (
    "real part of a rectangular complex number",
    GSL_REAL(test_rect_complex_number) == x,
    "real part of rectangular complex number does not match expected" 
  );

  mu_assert
  (
    "imaginary part of a rectangular complex number",
    GSL_IMAG(test_rect_complex_number) == y,
    "imaginary part of rectangular complex number does not match expected" 
  );

  GSL_SET_REAL(&test_rect_complex_number,y);
  GSL_SET_IMAG(&test_rect_complex_number,x);
  
  mu_assert
  (
    "real part of a rectangular complex number after redefinition",
    GSL_REAL(test_rect_complex_number) == y,
    "redefined real part of rectangular complex number does not match expected" 
  );

  mu_assert
  (
    "imaginary part of a rectangular complex number after redefinition",
    GSL_IMAG(test_rect_complex_number) == x,
    "redefined imaginary part of rectangular complex number does not match expected" 
  );

  return 0;
}  

static 
char * test_gsl_block_struct()
{
  gsl_block * block_p = gsl_block_alloc(4);

  mu_assert("Confirm that the block was created", block_p != NULL, "The block was not created."); 
 
  FILE * my_block_file = fopen ("my_block", "w+");


  gsl_block_fscanf (stdin, block_p);
  
  gsl_block_fprintf (stdout, block_p, "%g");

  gsl_block_fwrite(my_block_file, block_p);

  // gsl_block_free(block_p);

  // gsl_block_fprintf (stdout, block_p, "%g");

  // printf("block pointer after being freed: %p %zu\n", block_p, block_p->size);

  // mu_assert("Confirm that the block was destroyed", block_p == NULL, "The block was not destroyed."); 

  return 0;
}

static
char * test_gsl_vector_struct()
{

}
static 
char * all_tests ()
{
  mu_run_test("Test of GSL", test_gsl_0_order_bessel_function_of_the_first_kind); 
  mu_run_test("Test of Rectangular Complex Number Struct", 
              test_gsl_rectangular_complex_number_struct);
  mu_run_test("Test of Block Struct", 
              test_gsl_block_struct);
  return 0;
}

int 
main(int argc, char **argv)
{
  char * result = all_tests();
  if (result != 0)
  {
    printf("%s\n", result);
  }
  else
  {
    printf("\nAll Tests Passed\n");
  }
  printf("Tests run: %d\n", tests_run);

  return result != 0;
}

float 
round_to_6_places(float num)
{
  float nearest = roundf(num * 10000000) / 10000000;
  return nearest;
}
