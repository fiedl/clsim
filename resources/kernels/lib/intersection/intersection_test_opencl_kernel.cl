typedef float floating_t;

inline floating_t my_sqrt(floating_t a) {
  return sqrt(a);
}

inline floating_t sqr(floating_t a) {
  return a * a;
}

inline floating_t my_nan() {
  return NAN;
}

inline bool my_is_nan(floating_t a) {
  return isnan(a);
}

#include "intersection.c"

__kernel void test_intersection_with_gpu(__global floating_t* results)
{
  printf("Hello from the kernel!\n");

  IntersectionProblemParameters_t parameters = {
    0.0, 0.5, // A
    5.0, 0.5, // B
    1.0, 1.0, // M
    0.5       // r
  };

  results[0] = number_of_intersections(parameters);
  results[1] = intersection_x1(parameters);
  results[2] = intersection_y1(parameters);
  results[3] = intersection_x2(parameters);
  results[4] = intersection_y2(parameters);
  results[5] = intersection_s1(parameters);
  results[6] = intersection_s2(parameters);

}
