# Intersection

This code is used in gpu kernels to determine intersection points of photon trajectories and hole-ice cylinders in the x-y plane.

![image](https://user-images.githubusercontent.com/1679688/36549728-c613aa2c-17f3-11e8-856e-f00f4f32288f.png)

![math](https://user-images.githubusercontent.com/1679688/36590656-a6b4b940-188f-11e8-8f14-f0955b154204.png)

## Usage

Consider a trajectory from point A to point B and a circle around point M with radius r. The following snippet shows how to calculate the number of intersections and the intersection points X1 and X2. For a working example, please refer to the [intersection_tests.c](intersection_tests.c) file.

```c
#include "intersection.c"

int main()
{

  IntersectionProblemParameters_t p = {
    0.0, 0.5,              // A
    1.0, 1.0,              // M
    0.5,                   // r
    {1.0, 0.0, 0.0, 0.0},  // normalized direction from A to B
    5.0                    // distance from A to B
  };

  // The number of intersections is determined by the
  // intersection discriminant.
  if (intersection_discriminant(p) < 0) {
    printf("no intersections\n");
  } else if (intersection_discriminant(p) == 0) {
    printf("tangent point at (%e,%e)\n",
        intersection_x1(p), intersection_y1(p));
  } else {
    printf("intersection point 1: (%e,%e)\n",
        intersection_x1(p), intersection_y1(p));
    printf("intersection point 2: (%e,%e)\n",
        intersection_x2(p), intersection_y2(p));
  }

  return 0;
}
```

For technical reasons, `distance` is expected as 4-vector (`floating4_t`). As long as it is normalized, you may use just the first two components (x,y) or the first three components (x,z,z). If using three components, make sure that `distance` is the distance in three dimensions rather than the xy-projected distance.

To get the intersection points as ratio of the trajectory length AB, use the `intersection_s1` and `intersection_s2` functions:

```c
printf("AX1 / AB = %e\n", intersection_s1(p));
printf("AX2 / AB = %e\n", intersection_s2(p));
```

For an example how to use this with OpenCL on a GPU, have a look at [intersection_opencl_tests_kernel.cl](intersection_opencl_tests_kernel.cl) and the code that launches the kernel in [intersection_opencl_tests.c](intersection_opencl_tests.c).

You may also make use of these helper functions:

```c
inline bool intersecting_trajectory_starts_inside(IntersectionProblemParameters_t p)
inline bool intersecting_trajectory_starts_outside(IntersectionProblemParameters_t p)
inline bool intersecting_trajectory_ends_inside(IntersectionProblemParameters_t p)
```

## Requirements

### Data types

In order to be able to adjust the accuracy, the type of the variables and return values is `floating_t`. You need to define this type either to be `float` or `double`:

```c
typedef double floating_t;
```

Also, a 4-vector type is required. In OpenCL, you may just define:

```c
typedef float4 floating4_t;
```

When not using OpenCL, you need to define `floating4_t` as struct:

```c
struct floating4_t {
  floating_t x;
  floating_t y;
  floating_t z;
  floating_t w;
};
```

### Squares and square roots

The calculation makes use of a `my_sqrt` function that needs to be defined by the main program to call the apropriate square-root function. This mechanism can be used to switch between math libraries or native gpu math functions. Example:

```c
inline floating_t my_sqrt(floating_t a) {return sqrt(a);}
```

Also, a `sqr` function is needed to calculate squared numbers.

```c
inline floating_t sqr(floating_t a) {return a * a;}
```

### Not-a-number

A function `my_nan` is needed to define the not-a-number result. Also, a function `my_is_nan` is needed to check whether a variable is not a number. For example, math.h provides a `NAN`.

```c
inline floating_t my_nan() { return NAN; }
inline bool my_is_nan(floating_t a) { return (a != a); }
inline bool my_is_nan(floating_t a) { return isnan(a); }
```

Make sure not to use `isnan(a)` as this is buggy with some gpu drivers. See https://github.com/fiedl/hole-ice-study/issues/16.

### Minimize

A function `min` is required that returns the lesser value of two arguments.

```c
#include "math.h"
inline floating_t min(floating_t a, floating_t b) { return fmin(a, b); }
```

### Dot product

A `dot` function is required that calculates the scalar product of two vectors.

In OpenCL, `dot` is automatically defined. When not using OpenCL, you may define the scalar product directly:

```c
inline floating_t dot(floating4_t a, floating4_t b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
```

## Installation and Tests

To install this script on your development machine and run the automated tests, you may follow the these steps:

```bash
git clone git@github.com:fiedl/clsim.git
cd clsim/resources/kernels/lib/intersection
make test
```

## Author

License: MIT

Author: Sebastian Fiedlschuster, 2014, 2017-2018

Repo: https://github.com/fiedl/clsim

See also: https://github.com/fiedl/hole-ice-study
