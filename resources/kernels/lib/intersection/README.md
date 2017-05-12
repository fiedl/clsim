# Intersection

This code is used in gpu kernels to determine intersection points of photon trajectories and hole-ice cylinders in the x-y plane.

## Usage

Consider a trajectory from point A to point B and a circle around point M with radius r. The following snippet shows how to calculate the number of intersections and the intersection points X1 and X2. For a working example, please refer to the [intersection_tests.c](intersection_tests.c) file.

```c
#include "intersection.c"

int main()
{

  IntersectionProblemParameters_t parameters = {
    0.0, 0.5, // A
    5.0, 0.5, // B
    1.0, 1.0, // M
    0.5       // r
  };

  printf("Number of intersections: %i\n", number_of_intersections(parameters));
  printf("Intersection point 1: (%e,%e)\n", intersection_x1(parameters), intersection_y1(parameters));
  printf("Intersection point 2: (%e,%e)\n", intersection_x2(parameters), intersection_y2(parameters));
  return 0;
}
```

To get the intersection points as ratio of the trajectory length AB, use the `intersection_s1` and `intersection_s2` functions, which is useful, for example, when applying this to a projection of a 3-dimensional problem to the x-y plane.

```c
printf("AX1 / AB = %e\n", intersection_s1(parameters));
printf("AX2 / AB = %e\n", intersection_s2(parameters));
```

For an example how to use this with OpenCL on a GPU, have a look at [intersection_opencl_tests_kernel.cl](intersection_opencl_tests_kernel.cl) and the code that launches the kernel in [intersection_opencl_tests.c](intersection_opencl_tests.c).

You may also make use of these helper functions:

```c
floating_t squared_distance_from_center(floating_t X, floating_t Y, floating_t MX, floating_t MY)
bool intersecting_trajectory_starts_inside(IntersectionProblemParameters_t p)
inline bool intersecting_trajectory_starts_outside(IntersectionProblemParameters_t p)
floating_t intersection_ratio_inside(IntersectionProblemParameters_t p)
```

The most helpful one might be the latter, which calculates the fraction of the trajectory that is inside the circle.

## Requirements

In order to be able to adjust the accuracy, the type of the variables and return values is `floating_t`. You need to define this type either to be `float` or `double`:

```c
typedef double floating_t;
```

The calculation makes use of a `my_sqrt` function that needs to be defined by the main program to call the apropriate square-root function. This mechanism can be used to switch between math libraries or native gpu math functions. Example:

```c
inline floating_t my_sqrt(floating_t a) {return sqrt(a);}
```

Also, a `sqr` function is needed to calculate squared numbers.

```c
inline floating_t sqr(floating_t a) {return a * a;}
```

A function `my_nan` is needed to define the not-a-number result. Also, a function `my_is_nan` is needed to check whether a variable is not a number. For example, math.h provides a `NAN`.

```c
inline floating_t my_nan() { return NAN; }
inline bool my_is_nan(floating_t a) { return isnan(a); }
```

## Installation

To install this script on your development machine and run the automated tests, you may follow the these steps:

```bash
git clone git@github.com:fiedl/clsim.git
cd clsim/resources/kernels/lib/intersection
make test
```

## Author

Author: Sebastian Fiedlschuster, 2014
Repo: https://github.com/fiedl/clsim

