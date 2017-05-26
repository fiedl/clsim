#include "../intersection/intersection.c"

void perform_hole_ice_correction_for(floating_t *distancePropagated, floating_t *distanceToAbsorption, floating_t *sca_step_left, floating_t *abs_lens_left, floating_t currentScaLen, floating_t currentAbsLen, floating_t hole_ice_scattering_length_factor, floating_t hole_ice_absorption_length_factor, IntersectionProblemParameters_t p)
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
  // https://github.com/fiedl/hole_ice_study
  // https://github.com/fiedl/clsim/tree/sf/master/resources/kernels/lib/hole_ice

  bool starts_in_hole_ice = intersecting_trajectory_starts_inside(p);
  unsigned int num_of_intersections = number_of_intersections(p);

  floating_t distance_correction_from_scattering = 0.0;
  floating_t distance_correction_from_absorption = 0.0;

  // Case 1: The trajectory is completely outside of the hole ice.
  // Thus, needs no correction.
  // // ((num_of_intersections == 0) && !starts_in_hole_ice)

  // FIXME: The following quantities are not correct for the absorption
  // problem. There really need to be two separate intersection problems,
  // one with `distancePropagated`, one with `distanceToAbsorption`.
  floating_t trajectory_length = intersection_trajectory_length(p);
  floating_t trajectory_length_within_hole_ice = intersection_trajectory_length_inside(p);

  // Case 2: The trajectory is completely within the hole ice.
  if ((num_of_intersections == 0) && starts_in_hole_ice) {
    distance_correction_from_scattering = (hole_ice_scattering_length_factor - 1) * trajectory_length;
    distance_correction_from_absorption = (hole_ice_absorption_length_factor - 1) * trajectory_length;
  }

  // Case 3: The trajectory begins outside, but ends inside the hole ice.
  if ((num_of_intersections == 1) && !starts_in_hole_ice) {
    distance_correction_from_scattering = (hole_ice_scattering_length_factor - 1) * trajectory_length_within_hole_ice;
    distance_correction_from_absorption = (hole_ice_absorption_length_factor - 1) * trajectory_length_within_hole_ice;
  }

  // Case 4: The trajectory begins inside, but ends outside the hole ice.
  if ((num_of_intersections == 1) && starts_in_hole_ice) {
    distance_correction_from_scattering = (1 - 1 / hole_ice_scattering_length_factor) * trajectory_length_within_hole_ice;
    distance_correction_from_absorption = (1 - 1 / hole_ice_absorption_length_factor) * trajectory_length_within_hole_ice;

    // TODO: Sub case if trajectory is too short.
  }

  // Case 5: The trajectory starts and ends outside, but passes through the hole ice.
  if ((num_of_intersections == 2) && !starts_in_hole_ice) {
    distance_correction_from_scattering = (1 - 1 / hole_ice_scattering_length_factor) * trajectory_length_within_hole_ice;
    distance_correction_from_absorption = (1 - 1 / hole_ice_absorption_length_factor) * trajectory_length_within_hole_ice;

    // TODO: Sub case if trajectory is too short.
  }

  // Case 6: The trajectory starts in one cylinder, ends in another one.
  // We ignore this case as the cylinders are too far away from each other.
  // If this case is reached, this should raise an error.
  if ((num_of_intersections == 2) && starts_in_hole_ice) {
    // TODO: Raise error.
  }

  // TODO: Raise error if `num_of_intersections` is anything else.

  // Now, apply corrections to the variables given by the pointers,
  // which will pass the corrections to the original variables outside
  // the scope of this function.
  *distancePropagated = *distancePropagated + distance_correction_from_scattering;
  *distanceToAbsorption = *distanceToAbsorption + distance_correction_from_absorption;
  *sca_step_left = *sca_step_left + (distance_correction_from_scattering / currentScaLen);
  *abs_lens_left = *abs_lens_left + (distance_correction_from_absorption / currentAbsLen);

}