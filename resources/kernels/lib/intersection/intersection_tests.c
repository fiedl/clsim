#include <stdio.h>
#include <stdbool.h>
#include "math.h"
#include "intersection_tests.h"
#include "intersection.c"

inline floating_t my_sqrt(floating_t a) {return sqrt(a);}
inline floating_t sqr(floating_t a) {return a * a;}
inline floating_t my_nan() { return NAN; }
inline bool my_is_nan(floating_t a) { return isnan(a); }


IntersectionProblemParameters_t intersection_parameters()
{
  IntersectionProblemParameters_t parameters = {
    0.0, 0.5, // A
    5.0, 0.5, // B
    1.0, 1.0, // M
    0.5       // r
  };
  return parameters;
}

int test_intersection()
{
  IntersectionProblemParameters_t parameters = intersection_parameters();
  printf("Number of intersections: %i\n", number_of_intersections(parameters));
  printf("Intersection point 1: (%e,%e)\n", intersection_x1(parameters), intersection_y1(parameters));
  printf("Intersection point 2: (%e,%e)\n", intersection_x2(parameters), intersection_y2(parameters));
  return 0;
}

int test_intersection_speed()
{
  IntersectionProblemParameters_t parameters = intersection_parameters();
  for (int i = 0; i < 1e6; i++) {
    number_of_intersections(parameters);
    intersection_x1(parameters);
    intersection_y1(parameters);
    intersection_x2(parameters);
    intersection_y2(parameters);
  }
  return 0;
}

int main()
{
  test_intersection();
  test_intersection_speed();
  return 0;
}

