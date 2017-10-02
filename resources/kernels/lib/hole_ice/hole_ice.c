#ifndef HOLE_ICE_C
#define HOLE_ICE_C

#include "../intersection/intersection.c"

inline floating_t hole_ice_distance_correction(floating_t distance, floating_t interaction_length_factor, IntersectionProblemParameters_t p)
{
  // Depending on the fraction of the distance the photon is traveling
  // within the hole ice, there are six cases to consider.
  //
  // N  denotes the number of intersections.
  // H  means that the trajectory starts within hole ice.
  // !H means that the trajectory starts outside of the hole ice.
  //
  // Case 1: The trajectory is completely outside of the hole ice (!H, N=0).
  // Case 2: The trajectory is completely within the hole ice (H, N=0).
  // Case 3: The trajectory begins outside, but ends inside the hole ice (!H, N=1).
  // Case 4: The trajectory begins inside, but ends outside the hole ice (H, N=1).
  // Case 5: The trajectory starts and ends outside, but passes through the hole ice (!H, N=2).
  // Case 6: The trajectory begins within one hole-ice cylinder, passes through
  //           normal ice and ends in another hole-ice cylinder (H, N=2).
  //
  // For further information, please have look at:
  // https://github.com/fiedl/clsim/tree/sf/master/resources/kernels/lib/hole_ice

  bool starts_in_hole_ice = intersecting_trajectory_starts_inside(p);
  unsigned int num_of_intersections = number_of_intersections(p);

  // TODO: Histogram to find out which case is most probable.
  // The probable ones need to come first for best performane.

#ifdef PRINTF_ENABLED
  printf("PARAMETERS:\n");
  printf("intersecting_trajectory_starts_inside = %d\n", intersecting_trajectory_starts_inside(p));
  printf("intersection_ratio_inside = %f\n", intersection_ratio_inside(p));
  printf("intersection_alpha = %f\n", intersection_alpha(p));
  printf("intersection_beta = %f\n", intersection_beta(p));
  printf("intersection_gamma = %f\n", intersection_gamma(p));
  printf("intersection_discriminant = %f\n", intersection_discriminant(p));
  printf("intersection_s1 = %f\n", intersection_s1(p));
  printf("intersection_s2 = %f\n", intersection_s2(p));
  printf("intersecting_trajectory_starts_inside = %d\n", intersecting_trajectory_starts_inside(p));
  printf("number_of_intersections = %i\n", number_of_intersections(p));
  printf("squared radius = %f\n", sqr(p.r));
  printf("squared_distance_from_center = %f\n", squared_distance_from_center(p.ax, p.ay, p.mx, p.my));
  printf("A = (%f, %f)\n", p.ax, p.ay);
  printf("B = (%f, %f)\n", p.bx, p.by);
  printf("M = (%f, %f)\n", p.mx, p.my);
  printf("r = %f\n", p.r);
  printf("distance = %f\n", distance);
#endif

  // Case 1: The trajectory is completely outside of the hole ice.
  // Thus, needs no correction.
  if ((num_of_intersections == 0) && !starts_in_hole_ice) {
    return 0;
  }

  // Case 2: The trajectory is completely within the hole ice.
  if ((num_of_intersections == 0) && starts_in_hole_ice) {
    return (interaction_length_factor - 1.0) * distance;
  }

  const floating_t distance_within_hole_ice = distance * intersection_ratio_inside(p);

#ifdef PRINTF_ENABLED
  printf("distance_within_hole_ice = %f\n", distance_within_hole_ice);
#endif

  // Case 3: The trajectory begins outside, but ends inside the hole ice.
  if ((num_of_intersections == 1) && !starts_in_hole_ice) {
    return (interaction_length_factor - 1.0) * distance_within_hole_ice;
  }

  const floating_t ab = my_sqrt(sqr(p.bx - p.ax) + sqr(p.by - p.ay));
  const floating_t scalingFactorFor3D = distance / ab;

#ifdef PRINTF_ENABLED
  printf("ab = %f\n", ab);
  printf("scalingFactorFor3D = %f\n", scalingFactorFor3D);
#endif

  // Case 4: The trajectory begins inside, but ends outside the hole ice.
  if ((num_of_intersections == 1) && starts_in_hole_ice) {
    const floating_t ax = ab * intersection_s2(p);

#ifdef PRINTF_ENABLED
    printf("ax = %f\n", ax);
#endif

    if (interaction_length_factor * ab > ax) {
      return (1.0 - 1.0 / interaction_length_factor) * ax * scalingFactorFor3D;
    } else {
      // Scaled trajectory is too short for this case. Fall back to case 2.
      return (interaction_length_factor - 1.0) * distance;
    }
  }

  // Case 5: The trajectory starts and ends outside, but passes through the hole ice.
  if ((num_of_intersections == 2) && !starts_in_hole_ice) {
    const floating_t yb = ab * (1.0 - intersection_s1(p));
    const floating_t yx = ab * (intersection_s2(p) - intersection_s1(p));

#ifdef PRINTF_ENABLED
    printf("yb = %f\n", yb);
    printf("yx = %f\n", yx);
#endif

    if (interaction_length_factor * yb > yx) {
      return (1.0 - 1.0 / interaction_length_factor) * yx * scalingFactorFor3D;
    } else {
      // Scaled trajectory is too short for this case. Fall back to case 3.
      return (interaction_length_factor - 1.0) * yb * scalingFactorFor3D;
    }
  }

  // Case 6: The trajectory starts in one cylinder, ends in another one.
  // We ignore this case as the cylinders are too far away from each other.
  // If this case is reached, this should raise an error.
  if ((num_of_intersections == 2) && starts_in_hole_ice) {
#ifdef PRINTF_ENABLED
    printf("WARNING: INTERSECTION CASE 6 REACHED. This is not implemented, yet.");
#endif
    // TODO: Raise error.
  }

#ifdef PRINTF_ENABLED
  printf("WARNING: UNHANDLED INTERSECTION CASE. This point should not be reached.");
#endif
  // TODO: Raise error if `num_of_intersections` is anything else.

  return my_nan();
}

#endif