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

inline floating_t squared_distance_from_center(floating_t X, floating_t Y, floating_t MX, floating_t MY)
{
    return (sqr(MX - X) + sqr(MY - Y));
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

inline bool is_tangent(IntersectionProblemParameters_t p)
{
  return intersection_s2(p) == intersection_s1(p);
}

inline int number_of_intersections(IntersectionProblemParameters_t p)
{
  const floating_t d = intersection_discriminant(p);
  const floating_t s1 = intersection_s1(p);
  const floating_t s2 = intersection_s2(p);

  if (d < 0) return 0;
  else if (d == 0) return 1;
  else if (d > 0) {
    // Both intersection points behind the trajectory starting point A:
    if (my_is_nan(s1) && my_is_nan(s2)) return 0;

    // One intersection point behind the trajectory starting point A:
    if (my_is_nan(s1) || my_is_nan(s2)) {
      // Here, a numerical issue may arise: See 2017-05-27 notes.
      // If the photon starts outside and ends outside, then there can only be
      // 0 or 2 intersection points, not 1. This can happen when the
      // trajectory starts near the circle radius. The intersection
      // point does not exist but results from a numerical issue.
      const bool start_inside = intersecting_trajectory_starts_inside(p);
      const bool end_inside = intersecting_trajectory_ends_inside(p);

      if (start_inside && end_inside) {
        return 0;
      } else if (( ! start_inside) && ( ! end_inside)) {
        return 0;
      } else {
        return 1;
      }
    }
    // Both intersection points on the positive trajectory:
    return 2;
  }

#ifdef PRINTF_ENABLED
  printf("ERROR: THIS POINT SHOULD NOT BE REACHED. in number_of_intersections().\n");
#endif
  return my_nan();
}

inline floating_t intersection_ratio_inside(IntersectionProblemParameters_t p)
{
    bool starts_inside = intersecting_trajectory_starts_inside(p);
    int num_of_intersections = number_of_intersections(p);

    // printf("HOLE ICE - INTERSECTION\n");
    //printf(" -> num of intersections = %i\n", num_of_intersections);
    //if (starts_inside) printf(" -> starts inside.\n");

    //printf("intersection_alpha = %f\n", intersection_alpha(p));
    //printf("intersection_beta = %f\n", intersection_beta(p));
    //printf("intersection_gamma = %f\n", intersection_gamma(p));
    //printf("intersection_discriminant = %f\n", intersection_discriminant(p));
    //printf("intersection_s1 = %f\n", intersection_s1(p));

    if (( ! starts_inside ) && ( num_of_intersections == 0 ))
        return 0.0;
    if (( ! starts_inside ) && ( num_of_intersections == 1 )) {
        if (is_tangent(p)) {
            return 0.0;
        } else {
            return 1.0 - intersection_s1(p);
        }
    }
    if (( ! starts_inside ) && ( num_of_intersections == 2 ))
        return intersection_s2(p) - intersection_s1(p);
    if (( starts_inside ) && ( num_of_intersections == 0 ))
        return 1.0;
    if (( starts_inside ) && ( num_of_intersections == 1 ))
        return intersection_s2(p);

#ifdef PRINTF_ENABLED
    printf("ERROR. This point should not be reached! in intersection_ratio_inside().\n");
    printf("starts_inside = %d\n", starts_inside);
    printf("num_of_intersections = %i\n", num_of_intersections);
#endif
    return my_nan();
}

inline floating_t intersection_trajectory_length(IntersectionProblemParameters_t p)
{
  return my_sqrt(sqr(p.ax - p.bx) + sqr(p.ay - p.by));
}

inline floating_t intersection_trajectory_length_inside(IntersectionProblemParameters_t p)
{
  return intersection_trajectory_length(p) * intersection_ratio_inside(p);
}