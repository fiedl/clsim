#ifndef INTERSECTION_H
#define INTERSECTION_H

typedef struct IntersectionProblemParameters {

  // Input values
  //
  floating_t ax;
  floating_t ay;
  floating_t mx;
  floating_t my;
  floating_t r;
  floating4_t direction;
  floating_t distance;

  // Output values, which will be calculated in
  // `calculate_intersections()`.
  //
  floating_t discriminant;
  floating_t s1;
  floating_t s2;

} IntersectionProblemParameters_t;

#endif