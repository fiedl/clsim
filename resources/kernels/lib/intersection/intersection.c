typedef struct IntersectionProblemParameters {
  floating_t ax;
  floating_t ay;
  floating_t bx;
  floating_t by;
  floating_t mx;
  floating_t my;
  floating_t r;
} IntersectionProblemParameters_t;


floating_t intersection_alpha(IntersectionProblemParameters_t p)
{
  return sqr(p.by - p.ay) + sqr(p.bx - p.ax);
}

floating_t intersection_beta(IntersectionProblemParameters_t p)
{
  return 2 * p.ay * (p.by - p.ay)
    + 2 * p.ax * (p.bx - p.ax)
    - 2 * p.my * (p.by - p.ay)
    - 2 * p.mx * (p.bx - p.ax);
}

floating_t intersection_gamma(IntersectionProblemParameters_t p)
{
  return p.ay * p.ay - 2 * p.ay * p.my + p.my * p.my
    - p.r * p.r + p.ax * p.ax - 2 * p.ax * p.mx + p.mx * p.mx;
}

floating_t intersection_discriminant(IntersectionProblemParameters_t p)
{
  return sqr(intersection_beta(p)) - 4 * intersection_alpha(p) * intersection_gamma(p);
}

floating_t intersection_s(IntersectionProblemParameters_t p, int sign)
{
  floating_t scale_parameter = (-intersection_beta(p) +
    sign * my_sqrt(intersection_discriminant(p))) / 2 / intersection_alpha(p);

  // If the intersection point is outside, i.e. before or after the trajectory,
  // return 'not a number'.
  if (( scale_parameter < 0.0 ) || ( scale_parameter > 1.0 )) scale_parameter = my_nan();

  return scale_parameter;
}

floating_t intersection_s1(IntersectionProblemParameters_t p)
{
  return intersection_s(p, -1);
}

floating_t intersection_s2(IntersectionProblemParameters_t p)
{
  return intersection_s(p, +1);
}

floating_t intersection_x1(IntersectionProblemParameters_t p)
{
  return p.ax + (p.bx - p.ax) * intersection_s1(p);
}

floating_t intersection_x2(IntersectionProblemParameters_t p)
{
  return p.ax + (p.bx - p.ax) * intersection_s2(p);
}

floating_t intersection_y1(IntersectionProblemParameters_t p)
{
  return p.ay + (p.by - p.ay) * intersection_s1(p);
}

floating_t intersection_y2(IntersectionProblemParameters_t p)
{
  return p.ay + (p.by - p.ay) * intersection_s2(p);
}

int number_of_intersections(IntersectionProblemParameters_t p)
{
  floating_t d = intersection_discriminant(p);

  if (d < 0) return 0;
  else if (d == 0) return 1;
  else if (d > 0) {
    // Both intersection points behind the trajectory starting point A:
    if (my_is_nan(intersection_s2(p))) return 0;
    // One intersection point behind the trajectory starting point A:
    if (my_is_nan(intersection_s1(p))) return 1;
    // Both intersection points on the positive trajectory:
    return 2;
  }

  printf("ERROR: THIS POINT SHOULD NOT BE REACHED. in number_of_intersections().\n");
  return my_nan();
}
