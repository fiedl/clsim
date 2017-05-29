#ifndef HOLE_ICE_C
#define HOLE_ICE_C

#include "../intersection/intersection.c"

inline floating_t hole_ice_corrected_distance(floating_t trajectory_length, floating_t interaction_length_factor, IntersectionProblemParameters_t p)
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

  // Case 1: The trajectory is completely outside of the hole ice.
  // Thus, needs no correction.
  if ((num_of_intersections == 0) && !starts_in_hole_ice) {
    return trajectory_length;
  }

  // Case 2: The trajectory is completely within the hole ice.
  if ((num_of_intersections == 0) && starts_in_hole_ice) {
    return interaction_length_factor * trajectory_length;
  }

  const floating_t trajectory_length_within_hole_ice = trajectory_length * intersection_ratio_inside(p);

  // Case 3: The trajectory begins outside, but ends inside the hole ice.
  if ((num_of_intersections == 1) && !starts_in_hole_ice) {
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
    printf("trajectory_length = %f\n", trajectory_length);


    printf("trajectory_length_within_hole_ice = %f\n", trajectory_length_within_hole_ice);
    return trajectory_length + (interaction_length_factor - 1.0) * trajectory_length_within_hole_ice;
  }

  return trajectory_length;

  // Case 4: The trajectory begins inside, but ends outside the hole ice.
  if ((num_of_intersections == 1) && starts_in_hole_ice) {
    if (interaction_length_factor * trajectory_length > trajectory_length_within_hole_ice) {
      return trajectory_length + (1.0 - 1.0 / interaction_length_factor) * trajectory_length_within_hole_ice;
    } else {
      // Scaled trajectory is too short for this case. Fall back to case 2.
      return interaction_length_factor * trajectory_length;
    }
  }

  // Case 5: The trajectory starts and ends outside, but passes through the hole ice.
  if ((num_of_intersections == 2) && !starts_in_hole_ice) {
    // YB = ((1 - intersection_s1(p)) * trajectory_length)
    if (interaction_length_factor * ((1 - intersection_s1(p)) * trajectory_length) > trajectory_length_within_hole_ice) {
      return trajectory_length + (1 - 1 / interaction_length_factor) * trajectory_length_within_hole_ice;
    } else {
      // Scaled trajectory is too short for this case. Fall back to case 3.
      return trajectory_length + (interaction_length_factor - 1.0) * trajectory_length_within_hole_ice;
    }
  }

  // Case 6: The trajectory starts in one cylinder, ends in another one.
  // We ignore this case as the cylinders are too far away from each other.
  // If this case is reached, this should raise an error.
  if ((num_of_intersections == 2) && starts_in_hole_ice) {
    printf("WARNING: INTERSECTION CASE 6 REACHED. This is not implemented, yet.");
    // TODO: Raise error.
  }

  printf("WARNING: UNHANDLED INTERSECTION CASE. This point should not be reached.");
  // TODO: Raise error if `num_of_intersections` is anything else.

  return my_nan();
}

#endif