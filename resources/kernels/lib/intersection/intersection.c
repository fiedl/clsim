#include "intersection.h"

inline void calculate_intersections(IntersectionProblemParameters_t *p)
{
  // Step 1
  const floating4_t vector_AM = {p->mx - p->ax, p->my - p->ay, 0.0, 0.0};
  const floating_t xy_projection_factor = my_sqrt(1 - sqr(p->direction.z));
  const floating_t length_AMprime = dot(vector_AM, p->direction) / xy_projection_factor;

  // Step 2
  p->discriminant = sqr(p->r) - dot(vector_AM, vector_AM) + sqr(length_AMprime);

  // Step 3
  const floating_t length_XMprime = my_sqrt(p->discriminant);

  // Step 4
  const floating_t length_AX1 = length_AMprime - length_XMprime;
  const floating_t length_AX2 = length_AMprime + length_XMprime;
  p->s1 = length_AX1 / p->distance / xy_projection_factor;
  p->s2 = length_AX2 / p->distance / xy_projection_factor;

  // printf("  calculate_intersections()\n");
  // printf("    distance = %f\n", p->distance);
  // printf("    projection factor = %f\n", xy_projection_factor);
  // printf("    r = %f\n", p->r);
  // printf("    AM*AM = %f\n", dot(vector_AM, vector_AM));
  // printf("    vector_AM = (%f,%f)\n", vector_AM.x, vector_AM.y);
  // printf("    length_AMprime = %f\n", length_AMprime);
  // printf("    length_XMprime = %f\n", length_XMprime);
  // printf("    length_AX1 = %f\n", length_AX1);
  // printf("    length_AX2 = %f\n", length_AX2);

}

inline floating_t intersection_s1(IntersectionProblemParameters_t p)
{
  return p.s1;
}

inline floating_t intersection_s2(IntersectionProblemParameters_t p)
{
  return p.s2;
}

inline floating_t intersection_discriminant(IntersectionProblemParameters_t p)
{
  return p.discriminant;
}

inline floating_t intersection_x1(IntersectionProblemParameters_t p)
{
  if ((p.s1 > 0) && (p.s1 < 1))
    return p.ax + p.direction.x * p.distance * p.s1;
  else
    return my_nan();
}

inline floating_t intersection_x2(IntersectionProblemParameters_t p)
{
  if ((p.s2 > 0) && (p.s2 < 1))
    return p.ax + p.direction.x * p.distance * p.s2;
  else
    return my_nan();
}

inline floating_t intersection_y1(IntersectionProblemParameters_t p)
{
  if ((p.s1 > 0) && (p.s1 < 1))
    return p.ay + p.direction.y * p.distance * p.s1;
  else
    return my_nan();
}

inline floating_t intersection_y2(IntersectionProblemParameters_t p)
{
  if ((p.s2 > 0) && (p.s2 < 1))
    return p.ay + p.direction.y * p.distance * p.s2;
  else
    return my_nan();
}

inline bool intersecting_trajectory_starts_inside(IntersectionProblemParameters_t p)
{
  return (intersection_s1(p) <= 0) &&
      (intersection_s2(p) > 0) &&
      (intersection_discriminant(p) > 0);
}

inline bool intersecting_trajectory_starts_outside(IntersectionProblemParameters_t p)
{
  return ( ! intersecting_trajectory_starts_inside(p));
}

inline bool intersecting_trajectory_ends_inside(IntersectionProblemParameters_t p)
{
  return (intersection_s1(p) < 1) &&
      (intersection_s2(p) >= 1) &&
      (intersection_discriminant(p) > 0);
}

