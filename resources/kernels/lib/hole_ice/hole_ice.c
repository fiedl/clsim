#ifndef HOLE_ICE_C
#define HOLE_ICE_C

#include "hole_ice.h"
#include "../intersection/intersection.c"

inline bool is_between_zero_and_one(floating_t a) {
  return ((!my_is_nan(a)) && (a > 0.0) && (a < 1.0));
}

inline bool not_between_zero_and_one(floating_t a) {
  return (! is_between_zero_and_one(a));
}

inline unsigned int number_of_medium_changes(HoleIceProblemParameters_t p)
{
  // These are ordered by their frequency of occurrance in order
  // to optimize for performance.
  if (not_between_zero_and_one(p.entry_point_ratio) && not_between_zero_and_one(p.termination_point_ratio)) return 0;
  if (not_between_zero_and_one(p.entry_point_ratio) || not_between_zero_and_one(p.termination_point_ratio)) return 1;
  if (p.entry_point_ratio == p.termination_point_ratio) return 0; // tangent
  return 2;
}

inline floating_t hole_ice_distance_correction(HoleIceProblemParameters_t p)
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

  p.number_of_medium_changes = number_of_medium_changes(p);

  // TODO: Histogram to find out which case is most probable.
  // The probable ones need to come first for best performane.

  // #ifdef PRINTF_ENABLED
  //   printf("HOLE ICE DISTANCE CORRECTION DEBUG:\n");
  //   printf("  p.number_of_medium_changes = %i\n", p.number_of_medium_changes);
  //   printf("  p.distance = %f\n", p.distance);
  // #endif

  //printf("p.number_of_medium_changes = %i\n", p.number_of_medium_changes);

  // Case 1: The trajectory is completely outside of the hole ice.
  // Thus, needs no correction.
  if ((p.number_of_medium_changes == 0) && !p.starts_within_hole_ice) {
    printf("FALL 1\n");
    return 0;
  }

  const floating_t ac = p.distance * p.termination_point_ratio;

  if (p.starts_within_hole_ice) {

    // Case 4: The trajectory begins inside, but ends outside the hole ice.
    if (p.interaction_length_factor * p.distance > ac) {
      printf("FALL 4\n");

      return (1.0 - 1.0 / p.interaction_length_factor) * ac;

    // Case 2: The trajectory is completely within the hole ice.
    } else {
      printf("FALL 2\n");

      return (p.interaction_length_factor - 1.0) * p.distance;
    }

  } else {

    const floating_t yb = p.distance * (1.0 - p.entry_point_ratio);
    const floating_t yc = p.distance * (p.termination_point_ratio - p.entry_point_ratio);

    // Case 5: The trajectory starts and ends outside, but passes through the hole ice.
    if (p.interaction_length_factor * yb > yc) {
      printf("FALL 5\n");

      return (1.0 - 1.0 / p.interaction_length_factor) * yc;
    } else {
      printf("FALL 3\n");

      return (p.interaction_length_factor - 1.0) * p.distance * (1.0 - p.entry_point_ratio);
    }
  }

#ifdef PRINTF_ENABLED
  printf("WARNING: UNHANDLED INTERSECTION CASE. This point should not be reached.");
#endif
  // TODO: Raise error if `p.number_of_medium_changes` is anything else.

  printf("FALL 0\n");

  return my_nan();
}

inline floating_t hole_ice_distance_correction_for_intersection_problem(floating_t distance, floating_t interaction_length_factor, IntersectionProblemParameters_t p)
{
  // #ifdef PRINTF_ENABLED
  //   printf("HOLE ICE DISTANCE CORRECTION INTERSECTION PROBLEM DEBUG:\n");
  //   printf("  intersecting_trajectory_starts_inside(p) = %d\n", intersecting_trajectory_starts_inside(p));
  //   printf("  intersection_alpha(p) = %f\n", intersection_alpha(p));
  //   printf("  intersection_beta(p) = %f\n", intersection_beta(p));
  //   printf("  intersection_gamma(p) = %f\n", intersection_gamma(p));
  //   printf("  intersection_discriminant(p) = %f\n", intersection_discriminant(p));
  //   printf("  intersection_s1(p) = %f\n", intersection_s1(p));
  //   printf("  intersection_s2(p) = %f\n", intersection_s2(p));
  //   printf("  intersection_s1_for_lines(p) = %f\n", intersection_s1_for_lines(p));
  //   printf("  intersection_s2_for_lines(p) = %f\n", intersection_s2_for_lines(p));
  //   printf("  A = (%f, %f)\n", p.ax, p.ay);
  //   printf("  B = (%f, %f)\n", p.bx, p.by);
  //   printf("  M = (%f, %f)\n", p.mx, p.my);
  //   printf("  r = %f\n", p.r);
  // #endif

  HoleIceProblemParameters_t hip = {
    distance,
    interaction_length_factor,
    intersection_s1_for_lines(p), // entry_point_ratio
    intersection_s2_for_lines(p), // termination_point_ratio
    intersecting_trajectory_starts_inside(p), // starts_within_hole_ice
    0 // number_of_medium_changes (will be calculated)
  };
  return hole_ice_distance_correction(hip);
}

#ifdef HOLE_ICE_TEST_H
inline floating_t apply_hole_ice_correction(floating4_t photonPosAndTime, floating4_t photonDirAndWlen, unsigned int numberOfCylinders, floating4_t *cylinderPositionsAndRadii, floating_t holeIceScatteringLengthFactor, floating_t holeIceAbsorptionLengthFactor, floating_t *distancePropagated, floating_t *distanceToAbsorption)
#endif
#ifndef HOLE_ICE_TEST_H
inline floating_t apply_hole_ice_correction(floating4_t photonPosAndTime, floating4_t photonDirAndWlen, unsigned int numberOfCylinders, __constant floating4_t *cylinderPositionsAndRadii, floating_t holeIceScatteringLengthFactor, floating_t holeIceAbsorptionLengthFactor, floating_t *distancePropagated, floating_t *distanceToAbsorption)
#endif
{
  // The algorithm for the hole ice corrections is as follows:
  //
  // 1. intersection problem p = (vec A, vec B, vec M, r)
  //
  //    vec A is the position of the photon at the beginning of this
  //    simulation step.
  //
  //    vec B is the position of the photon at the end of this simulation
  //    step:
  //
  //        vec B = vec photonPosition + vec photonDirection * distancePropagated
  //
  // 2. distancePropagated += hole_ice_distance_correction(distancePropagated,
  //        holeIceScatteringLengthFactor, p)
  //
  // 3. Gamma = min(distancePropagated, distanceToAbsorption)
  //
  //    distancePropagated is the already corrected distance until the next
  //    scattering point.
  //
  //    Gamma is the distance between vec A and vec B that is used to calculate
  //    the hole ice correction for the distance to absorption. We need to take
  //    the min here, because, if the photon is scattered away before reaching
  //    the point of absorption, only a fraction of the path, i.e. Gamma rather
  //    then the whole distance to absorption, is affected by the hole ice.
  //
  // 4. p.vec B = vec photonPosition + vec photonDirection * Gamma
  //
  // 5. distanceToAbsorption += hole_ice_distance_correction(Gamma,
  //        holeIceAbsorptionLengthFactor, p)
  //
  // After these steps, both `distancePropagated` and `distanceToAbsorption`
  // have been properly corrected for this simulation step.

  // For some reason, there are photons with photonPosAndTime coordinates
  // nan. I will have to ignore them.
  // TODO: Why?
  if (!(my_is_nan(photonPosAndTime.x) || my_is_nan(*distancePropagated))) {

    // TODO: Remove those after debugging:
    const floating_t distancePropagatedBeforeCorrection = *distancePropagated;
    const floating_t distanceToAbsorptionBeforeCorrection = *distanceToAbsorption;

    for (unsigned int i = 0; i < numberOfCylinders; i++) {

      // Is the cylinder in range?
      if (sqr(photonPosAndTime.x - cylinderPositionsAndRadii[i].x) +
          sqr(photonPosAndTime.y - cylinderPositionsAndRadii[i].y) <=
            sqr(*distancePropagated + cylinderPositionsAndRadii[i].w /* radius */))
      {

        // TODO: Update algorithm description above.

        printf("HOLE ICE DEBUG:\n");
        printf("  *distancePropagated = %f\n", *distancePropagated);
        printf("  *distanceToAbsorption = %f\n", *distanceToAbsorption);

        IntersectionProblemParameters_t p = {
          photonPosAndTime.x,
          photonPosAndTime.y,
          photonPosAndTime.x + photonDirAndWlen.x * *distancePropagated,
          photonPosAndTime.y + photonDirAndWlen.y * *distancePropagated,
          cylinderPositionsAndRadii[i].x,
          cylinderPositionsAndRadii[i].y,
          cylinderPositionsAndRadii[i].w // radius
        };

        // Are intersection points possible?
        if (intersection_discriminant(p) > 0) {

          const floating_t scatteringEntryPointRatio = intersection_s1_for_lines(p);
          const floating_t scatterintTerminationPointRatio = intersection_s2_for_lines(p);

          HoleIceProblemParameters_t scatteringCorrectionParameters = {
            *distancePropagated,
            holeIceScatteringLengthFactor,
            scatteringEntryPointRatio, // entry_point_ratio
            scatterintTerminationPointRatio, // termination_point_ratio
            intersecting_trajectory_starts_inside(p), // starts_within_hole_ice
            0 // number_of_medium_changes (will be calculated)
          };

          const floating_t scaCorrection = hole_ice_distance_correction(scatteringCorrectionParameters);
          *distancePropagated += scaCorrection;

          printf("  SCATTERING CORRECTION:\n");
          printf("    scaCorrection = %f\n", scaCorrection);
          printf("    *distancePropagated = %f\n", *distancePropagated);
          printf("    intersection_s1(p) = %f\n", intersection_s1(p));
          printf("    intersection_s2(p) = %f\n", intersection_s2(p));
          printf("    intersection_s1_for_lines(p) = %f\n", intersection_s1_for_lines(p));
          printf("    intersection_s2_for_lines(p) = %f\n", intersection_s2_for_lines(p));
          printf("    intersection_discriminant(p) = %f\n", intersection_discriminant(p));
          printf("    entry_point_ratio = %f\n", scatteringCorrectionParameters.entry_point_ratio);
          printf("    termination_point_ratio = %f\n", scatteringCorrectionParameters.termination_point_ratio);
          printf("    number_of_medium_changes = %i\n", number_of_medium_changes(scatteringCorrectionParameters));
          if (scatteringCorrectionParameters.starts_within_hole_ice) {
            printf("    starts_within_hole_ice = true\n");
          } else {
            printf("    starts_within_hole_ice = false\n");
          }


          // For the absorption, there are special cases where the photon is scattered before
          // reaching either the first or the second absorption intersection point.
          floating_t absorptionEntryPointRatio;
          floating_t absorptionTerminationPointRatio;
          floating_t absCorrection = 0.0;
          if (!(not_between_zero_and_one(scatteringCorrectionParameters.entry_point_ratio) && !scatteringCorrectionParameters.starts_within_hole_ice)) {
            // The photon reaches the hole ice, i.e. the absorption correction
            // needs to be calculated.
            p.bx = photonPosAndTime.x + photonDirAndWlen.x * *distanceToAbsorption;
            p.by = photonPosAndTime.y + photonDirAndWlen.y * *distanceToAbsorption;

            // If the photon is scattered away before reaching the far and of
            // the hole ice, the affected trajectory is limited by the
            // point where the photon is scattered away.
            absorptionEntryPointRatio = intersection_s1_for_lines(p);
            absorptionTerminationPointRatio = min(
              *distancePropagated / *distanceToAbsorption,
              intersection_s2_for_lines(p)
            );

            HoleIceProblemParameters_t absorptionCorrectionParameters = {
              *distanceToAbsorption,
              holeIceAbsorptionLengthFactor,
              absorptionEntryPointRatio, // entry_point_ratio
              absorptionTerminationPointRatio, // termination_point_ratio
              intersecting_trajectory_starts_inside(p), // starts_within_hole_ice
              0 // number_of_medium_changes (will be calculated)
            };

            absCorrection = hole_ice_distance_correction(absorptionCorrectionParameters);

            printf("  ABSORPTION CORRECTION:\n");
            printf("    absCorrection = %f\n", absCorrection);
            printf("    *distanceToAbsorption = %f\n", *distanceToAbsorption);
            printf("    intersection_s1(p) = %f\n", intersection_s1(p));
            printf("    intersection_s2(p) = %f\n", intersection_s2(p));
            printf("    intersection_s1_for_lines(p) = %f\n", intersection_s1_for_lines(p));
            printf("    intersection_s2_for_lines(p) = %f\n", intersection_s2_for_lines(p));
            printf("    intersection_discriminant(p) = %f\n", intersection_discriminant(p));
            printf("    number_of_medium_changes = %i\n", number_of_medium_changes(absorptionCorrectionParameters));
            printf("    entry_point_ratio = %f\n", absorptionCorrectionParameters.entry_point_ratio);
            printf("    termination_point_ratio = %f\n", absorptionCorrectionParameters.termination_point_ratio);
            if (absorptionCorrectionParameters.starts_within_hole_ice) {
              printf("    starts_within_hole_ice = true\n");
            } else {
              printf("    starts_within_hole_ice = false\n");
            }

          }
          *distanceToAbsorption += absCorrection;

          printf(
            "NAN DEBUG: "
            "scaCorrection=%f, "
            "absCorrection=%f, "
            "photonPosAndTime=(%f,%f,%f,.), "
            "photonDirAndWlen=(%f,%f,%f,.), "
            "cylinderPositionsAndRadii={{%f,%f,%f,%f}}, "
            "holeIceScatteringLengthFactor=%f, "
            "holeIceAbsorptionLengthFactor=%f, "
            "distancePropagatedBeforeCorrection=%f, "
            "distanceToAbsorptionBeforeCorrection=%f"
            "\n",
            scaCorrection,
            absCorrection,
            photonPosAndTime.x, photonPosAndTime.y, photonPosAndTime.z,
            photonDirAndWlen.x, photonDirAndWlen.y, photonDirAndWlen.z,
            cylinderPositionsAndRadii[0].x, cylinderPositionsAndRadii[0].y, cylinderPositionsAndRadii[0].z, cylinderPositionsAndRadii[0].w,
            holeIceScatteringLengthFactor,
            holeIceAbsorptionLengthFactor,
            distancePropagatedBeforeCorrection,
            distanceToAbsorptionBeforeCorrection
          );
        }


                  // We don't need to calculate a correction for `abs_lens_left` and `sca_step_left`, because
                  // `abs_lens_left` is recalculated after the hole-ice code, and `sca_step_left` is not used
                  // for this loop anymore.

                  // TODO: Test the code from 2014 with the hole_ice_tests.c
                  // before deleting the following code from 2014.
                  // Is the old code also correct but more efficient?

            //floating_t trajectory_ratio_inside_of_the_cylinder =
            //    intersection_ratio_inside(p);
            //
            //if ((!trajectory_ratio_inside_of_the_cylinder == ZERO) &
            //    (!my_is_nan(trajectory_ratio_inside_of_the_cylinder))) {
            //
            //  // printf("HOLE ICE -> trajectory inside: %f\n",
            //  //    trajectory_ratio_inside_of_the_cylinder);
            //
            //  // The propagated distance and the absorpotion lengths left have
            //  // to be corrected for the modified ice-properties within the hole
            //  // ice along the part of the trajectory that is within the
            //  // hole-ice cylinder.
            //
            //  // printf(" -> distancePropagated before: %f\n",
            //  // distancePropagated);
            //
            //  // Correct for the modified scattering length.
            //  floating_t distanceInsideTheCylinder =
            //      distancePropagated * trajectory_ratio_inside_of_the_cylinder;
            //  distancePropagated -= distanceInsideTheCylinder *
            //                        (ONE / holeIceScatteringLengthFactor - ONE);
            //  if (distancePropagated < ZERO)
            //    distancePropagated = ZERO;
            //  sca_step_left -= distanceInsideTheCylinder *
            //                   (ONE / holeIceScatteringLengthFactor - ONE) /
            //                   (currentScaLen * holeIceScatteringLengthFactor);
            //  if (sca_step_left < ZERO)
            //    sca_step_left = ZERO;
            //  abs_lens_left += distanceInsideTheCylinder *
            //                   (ONE / holeIceScatteringLengthFactor - ONE) /
            //                   (currentAbsLen * holeIceAbsorptionLengthFactor);
            //
            //  // printf(" -> distancePropagated AFTER: %f\n",
            //  // distancePropagated);
            //
            //  // Correct for the modified absorption length.
            //  abs_lens_left -= distanceInsideTheCylinder *
            //                   (ONE / holeIceAbsorptionLengthFactor - ONE) /
            //                   (currentAbsLen * holeIceAbsorptionLengthFactor);
            //  if (abs_lens_left < ZERO)
            //    abs_lens_left = ZERO;
            //
//#ifdef PRINTF_ENABLED
//              if (my_is_nan(abs_lens_left)) {
//                printf("WARNING: THIS SHOULD NOT BE REACHED. abs_lens_left == "
//                       "nan!\n");
//                printf("distance inside = %f\n", distanceInsideTheCylinder);
//                printf("absorption factor = %f\n",
//                       holeIceAbsorptionLengthFactor);
//                printf("currentAbsLen = %f\n", currentAbsLen);
//                printf("holeIceScatteringLengthFactor = %f\n",
//                       holeIceScatteringLengthFactor);
//                printf("distancePropagated = %f\n", distancePropagated);
//                printf("trajectory_ratio_inside_of_the_cylinder = %f\n",
//                       trajectory_ratio_inside_of_the_cylinder);
//              }
//#endif
            //}
                }
              }
	          }

}

#endif
