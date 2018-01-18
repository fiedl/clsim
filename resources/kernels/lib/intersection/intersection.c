#include "intersection.h"

inline floating_t intersection_alpha(IntersectionProblemParameters_t p)
{
  return sqr(p.by - p.ay) + sqr(p.bx - p.ax);
}

inline floating_t intersection_beta(IntersectionProblemParameters_t p)
{
  return 2 * p.ay * (p.by - p.ay)
    + 2 * p.ax * (p.bx - p.ax)
    - 2 * p.my * (p.by - p.ay)
    - 2 * p.mx * (p.bx - p.ax);
}

inline floating_t intersection_gamma(IntersectionProblemParameters_t p)
{
  return p.ay * p.ay - 2 * p.ay * p.my + p.my * p.my
    - p.r * p.r + p.ax * p.ax - 2 * p.ax * p.mx + p.mx * p.mx;
}

inline floating_t intersection_discriminant(IntersectionProblemParameters_t p)
{
  return sqr(intersection_beta(p)) - 4 * intersection_alpha(p) * intersection_gamma(p);
}

inline floating_t intersection_s_for_lines(IntersectionProblemParameters_t p, int sign)
{
  return (-intersection_beta(p) +
      sign * my_sqrt(intersection_discriminant(p))) / 2 / intersection_alpha(p);
}

inline floating_t intersection_s1_for_lines(IntersectionProblemParameters_t p)
{
  return intersection_s_for_lines(p, -1);
}

inline floating_t intersection_s2_for_lines(IntersectionProblemParameters_t p)
{
  return intersection_s_for_lines(p, +1);
}

inline floating_t intersection_s(IntersectionProblemParameters_t p, int sign)
{
  floating_t scale_parameter = intersection_s_for_lines(p, sign);

  // If the intersection point is outside, i.e. before or after the trajectory,
  // return 'not a number'.
  if (( scale_parameter <= 0.0 ) || ( scale_parameter >= 1.0 )) scale_parameter = my_nan();

  return scale_parameter;
}

inline floating_t intersection_s1(IntersectionProblemParameters_t p)
{
  return intersection_s(p, -1);
}

inline floating_t intersection_s2(IntersectionProblemParameters_t p)
{
  return intersection_s(p, +1);
}

inline floating_t intersection_x1(IntersectionProblemParameters_t p)
{
  return p.ax + (p.bx - p.ax) * intersection_s1(p);
}

inline floating_t intersection_x2(IntersectionProblemParameters_t p)
{
  return p.ax + (p.bx - p.ax) * intersection_s2(p);
}

inline floating_t intersection_y1(IntersectionProblemParameters_t p)
{
  return p.ay + (p.by - p.ay) * intersection_s1(p);
}

inline floating_t intersection_y2(IntersectionProblemParameters_t p)
{
  return p.ay + (p.by - p.ay) * intersection_s2(p);
}

inline bool intersecting_trajectory_starts_inside(IntersectionProblemParameters_t p)
{
  return (intersection_s1_for_lines(p) <= 0) &&
      (intersection_s2_for_lines(p) > 0) &&
      (intersection_discriminant(p) > 0);
}

inline bool intersecting_trajectory_starts_outside(IntersectionProblemParameters_t p)
{
  return ( ! intersecting_trajectory_starts_inside(p));
}

inline bool intersecting_trajectory_ends_inside(IntersectionProblemParameters_t p)
{
  return (intersection_s1_for_lines(p) < 1) &&
      (intersection_s2_for_lines(p) >= 1) &&
      (intersection_discriminant(p) > 0);
}

