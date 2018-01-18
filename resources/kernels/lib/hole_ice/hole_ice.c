#ifndef HOLE_ICE_C
#define HOLE_ICE_C

#include "hole_ice.h"
#include "../intersection/intersection.c"

inline unsigned int number_of_medium_changes(HoleIceProblemParameters_t p)
{
  // These are ordered by their frequency of occurrance in order
  // to optimize for performance.
  if (my_is_nan(p.entry_point_ratio) && my_is_nan(p.termination_point_ratio)) return 0;
  if (my_is_nan(p.entry_point_ratio) || my_is_nan(p.termination_point_ratio)) return 1;
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
    return 0;
  }

  // Case 2: The trajectory is completely within the hole ice.
  if ((p.number_of_medium_changes == 0) && p.starts_within_hole_ice) {
    return (p.interaction_length_factor - 1.0) * p.distance;
  }

  // p.distance_ratio_inside_hole_ice = distance_ratio_inside_hole_ice(p);
  // const floating_t distance_within_hole_ice = p.distance * p.distance_ratio_inside_hole_ice;

  // #ifdef PRINTF_ENABLED
  //   printf("HOLE ICE DISTANCE CORRECTION DEBUG:\n");
  //   printf("  p.distance_ratio_inside_hole_ice = %f\n", p.distance_ratio_inside_hole_ice);
  //   printf("  p.distance_within_hole_ice = %f\n", distance_within_hole_ice);
  // #endif

  // Case 3: The trajectory begins outside, but ends inside the hole ice.
  if ((p.number_of_medium_changes == 1) && !p.starts_within_hole_ice) {
    return (p.interaction_length_factor - 1.0) * p.distance * (1.0 - p.entry_point_ratio);
  }

  const floating_t ab = p.distance;

  // Case 4: The trajectory begins inside, but ends outside the hole ice.
  if ((p.number_of_medium_changes == 1) && p.starts_within_hole_ice) {
    const floating_t ac = ab * p.termination_point_ratio;
    if (p.interaction_length_factor * ab > ac) {
      return (1.0 - 1.0 / p.interaction_length_factor) * ac;
    } else {
      // Scaled trajectory is too short for this case. Fall back to case 2.
      return (p.interaction_length_factor - 1.0) * p.distance;
    }
  }

  // Case 5: The trajectory starts and ends outside, but passes through the hole ice.
  if ((p.number_of_medium_changes == 2) && !p.starts_within_hole_ice) {
    const floating_t yb = ab * (1.0 - p.entry_point_ratio);
    const floating_t yc = ab * (p.termination_point_ratio - p.entry_point_ratio);

    if (p.interaction_length_factor * yb > yc) {
      return (1.0 - 1.0 / p.interaction_length_factor) * yc;
    } else {
      // Scaled trajectory is too short for this case. Fall back to case 3.
      return (p.interaction_length_factor - 1.0) * yb;
    }
  }

  // Case 6: The trajectory starts in one cylinder, ends in another one.
  // We ignore this case as the cylinders are too far away from each other.
  // If this case is reached, this should raise an error.
  if ((p.number_of_medium_changes == 2) && p.starts_within_hole_ice) {
#ifdef PRINTF_ENABLED
    printf("WARNING: INTERSECTION CASE 6 REACHED. This is not implemented, yet.");
#endif
    // TODO: Raise error.
  }

#ifdef PRINTF_ENABLED
  printf("WARNING: UNHANDLED INTERSECTION CASE. This point should not be reached.");
#endif
  // TODO: Raise error if `p.number_of_medium_changes` is anything else.

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
    intersection_s1(p), // entry_point_ratio
    intersection_s2(p), // termination_point_ratio
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

    for (unsigned int i = 0; i < numberOfCylinders; i++) {

      // Is the cylinder in range?
      if (sqr(photonPosAndTime.x - cylinderPositionsAndRadii[i].x) +
          sqr(photonPosAndTime.y - cylinderPositionsAndRadii[i].y) <=
            sqr(*distancePropagated + cylinderPositionsAndRadii[i].w /* radius */))
      {

        const floating_t xyProjectionFactor = my_sqrt(1 - sqr(photonDirAndWlen.z));
        const floating_t projectedDistancePropagated = *distancePropagated * xyProjectionFactor;

        // TODO: Update algorithm description above.

        // printf("HOLE ICE DEBUG:\n");
        // printf("  *distancePropagated = %f\n", *distancePropagated);
        // printf("  *distanceToAbsorption = %f\n", *distanceToAbsorption);

        IntersectionProblemParameters_t p = {
          photonPosAndTime.x,
          photonPosAndTime.y,
          photonPosAndTime.x + photonDirAndWlen.x * projectedDistancePropagated,
          photonPosAndTime.y + photonDirAndWlen.y * projectedDistancePropagated,
          cylinderPositionsAndRadii[i].x,
          cylinderPositionsAndRadii[i].y,
          cylinderPositionsAndRadii[i].w // radius
        };

        HoleIceProblemParameters_t scatteringCorrectionParameters = {
          *distancePropagated,
          holeIceScatteringLengthFactor,
          intersection_s1(p), // entry_point_ratio
          intersection_s2(p), // termination_point_ratio
          intersecting_trajectory_starts_inside(p), // starts_within_hole_ice
          0 // number_of_medium_changes (will be calculated)
        };

        const floating_t scaCorrection = hole_ice_distance_correction(scatteringCorrectionParameters);
        *distancePropagated += scaCorrection;

        // printf("  SCATTERING CORRECTION:\n");
        // printf("    scaCorrection = %f\n", scaCorrection);
        // printf("    *distancePropagated = %f\n", *distancePropagated);
        // printf("    intersection_s1(p) = %f\n", intersection_s1(p));
        // printf("    intersection_s2(p) = %f\n", intersection_s2(p));
        // printf("    intersection_discriminant(p) = %f\n", intersection_discriminant(p));
        // printf("    number_of_intersections(p) = %i\n", number_of_intersections(p));
        // printf("    number_of_medium_changes = %i\n", scatteringCorrectionParameters.number_of_medium_changes);
        // printf("    distance_ratio_inside_hole_ice = %f\n", scatteringCorrectionParameters.distance_ratio_inside_hole_ice);
        // printf("    entry_point_ratio = %f\n", scatteringCorrectionParameters.entry_point_ratio);
        // printf("    termination_point_ratio = %f\n", scatteringCorrectionParameters.termination_point_ratio);
        // if (intersecting_trajectory_starts_inside(p)) { printf("    intersecting_trajectory_starts_inside = true\n"); }

        const floating_t projectedDistanceToAbsorption = *distanceToAbsorption * xyProjectionFactor;
        p.bx = photonDirAndWlen.x * projectedDistanceToAbsorption;
        p.by = photonDirAndWlen.y * projectedDistanceToAbsorption;
        const floating_t absorptionTerminationPointRatio = my_is_nan(intersection_s2(p)) ? my_nan() : min(
          *distancePropagated / *distanceToAbsorption,
          intersection_s2(p)
        );

        HoleIceProblemParameters_t absorptionCorrectionParameters = {
          *distanceToAbsorption,
          holeIceAbsorptionLengthFactor,
          intersection_s1(p), // entry_point_ratio
          absorptionTerminationPointRatio, // termination_point_ratio
          intersecting_trajectory_starts_inside(p), // starts_within_hole_ice
          0 // number_of_medium_changes (will be calculated)
        };

        const floating_t absCorrection = hole_ice_distance_correction(absorptionCorrectionParameters);
        *distanceToAbsorption += absCorrection;

        printf("  ABSORPTION CORRECTION:\n");
        printf("    absCorrection = %f\n", absCorrection);
        printf("    *distanceToAbsorption = %f\n", *distanceToAbsorption);
        // printf("    intersection_s1(p) = %f\n", intersection_s1(p));
        // printf("    intersection_s2(p) = %f\n", intersection_s2(p));
        // printf("    intersection_discriminant(p) = %f\n", intersection_discriminant(p));
        // printf("    number_of_intersections(p) = %i\n", number_of_intersections(p));
        // printf("    number_of_medium_changes = %i\n", absorptionCorrectionParameters.number_of_medium_changes);
        // printf("    distance_ratio_inside_hole_ice = %f\n", absorptionCorrectionParameters.distance_ratio_inside_hole_ice);
        // printf("    entry_point_ratio = %f\n", absorptionCorrectionParameters.entry_point_ratio);
        // printf("    termination_point_ratio = %f\n", absorptionCorrectionParameters.termination_point_ratio);
        //
        // if (intersecting_trajectory_starts_inside(p)) {
        //   printf("    intersecting_trajectory_starts_inside = true\n");
        // } else {
        //   printf("    intersecting_trajectory_starts_inside = false\n");
        // }
        // printf("    squared_distance_from_center(p.ax, p.ay, p.mx, p.my) = %f\n", squared_distance_from_center(p.ax, p.ay, p.mx, p.my));
        // printf("    sqr(p.r) = %f\n", sqr(p.r));

//
//
//
//
// 		              //printf("distancePropagated = %f\n", distancePropagated);
//
//                   const floating_t scaCorrection = hole_ice_distance_correction(
//                     *distancePropagated,
//                     holeIceScatteringLengthFactor,
//                     p
//                   );
//
// 		              //printf("scaCorrection = %f\n", scaCorrection);
//
//                   // floating_t sum = *distancePropagated;
//                   //
//                   // floating_t d;
//                   // d = -0.1; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -0.2; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -0.3; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -0.4; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -0.5; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -0.6; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -0.7; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -0.8; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -0.9; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.0; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.1; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.2; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.3; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.4; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.5; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.6; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.7; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.8; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -1.9; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//                   // d = -2.0; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
//
//                   *distancePropagated += scaCorrection; // Seltsam. Ohne diese Zeile geht es.
//
// 		              //const floating_t sum = distancePropagated + scaCorrection;
// 		              //printf("sum = %f\n", sum);
//                   //*distancePropagated = sum;
//
// 		              //const floating_t foo = 0.9 * distancePropagated;  // geht
// 		              //distancePropagated = distancePropagated + foo;
//
// 		              //distancePropagated = distancePropagated; // geht.
// 		              //distancePropagated = distancePropagated + 0.0 * scaCorrection; // geht nicht.
//
// 	      	        // Wenn auch nur eine Zeile mit "geht nicht" drin ist, wird keiner der prints ausgeführt.
//                   //if (scaCorrection < -1.0) distancePropagated -= 1.0; // geht.
//
//                   // distancePropagated += round(scaCorrection); // geht nicht.
//
//                   // geht: d.h. es liegt nicht am Schreiben in die Variable, sondern
//                   // irgendwo bei der Weiterverarbeitung.
//                   //distancePropagated += scaCorrection;
//                   //distancePropagated = 1.0;
//
//                   // // geht:
//                   // distancePropagated *= holeIceScatteringLengthFactor;
//                   // printf("holeIceScatteringLengthFactor = %f\n", holeIceScatteringLengthFactor);
//
//                   //distancePropagated += 0 * scaCorrection;
//                   //distancePropagated *= holeIceScatteringLengthFactor;
//                   //printf("distancePropagated = %f NEW\n", distancePropagated);
//
//
//
//                   //printf("distancePropagated = %f before resetting to 1.0\n", distancePropagated);
//                   //distancePropagated = 1.0;
//
//
//                   // const floating_t distanceToConsiderForAbsorption = min( // Gamma
//                   //   *distancePropagated, // after hole-ice correction
//                   //   *distanceToAbsorption // before hole-ice correction
//                   // );
//                   //
//                   // p.bx = photonPosAndTime.x + photonDirAndWlen.x * distanceToConsiderForAbsorption;
//                   // p.by = photonPosAndTime.y + photonDirAndWlen.y * distanceToConsiderForAbsorption;
//                   //
//                   // const floating_t absCorrection = hole_ice_distance_correction(
//                   //   distanceToConsiderForAbsorption,
//                   //   holeIceAbsorptionLengthFactor,
//                   //   p
//                   // );
//
//                   const floating_t terminus = min(
//                     *distanceToAbsorption,
//                     *distancePropagated
//                   );
//
//                   p.bx = photonPosAndTime.x + photonDirAndWlen.x * *distanceToAbsorption;
//                   p.by = photonPosAndTime.y + photonDirAndWlen.y * *distanceToAbsorption;
//
//                   const floating_t absCorrection = hole_ice_distance_correction(
//                     *distanceToAbsorption,
//                     holeIceAbsorptionLengthFactor,
//                     p
//                   );
//
//                   printf("absCorrection = %f\n", absCorrection);
//                   printf("(*distanceToAbsorption - *distancePropagated) / holeIceAbsorptionLengthFactor = %f\n", (*distanceToAbsorption - *distancePropagated) / holeIceAbsorptionLengthFactor);
//
//                   *distanceToAbsorption += absCorrection;
//
//                   if (*distanceToAbsorption > *distancePropagated) {
//                     *distanceToAbsorption += (*distanceToAbsorption - *distancePropagated) / holeIceAbsorptionLengthFactor;
//                   }
//
//                   // sum = 0;
//                   // for (floating_t d = -0.0; d > -5.0; d -= 0.1) {
//                   //   if ((absCorrection < d) && (absCorrection > (d - 0.1))) sum = *distanceToAbsorption + d;
//                   // }
//                   // //printf("  absCorrection = %f\n", absCorrection);
//                   // printf("  sum = %f\n", sum);
//                   // *distanceToAbsorption = sum;
//
//                   //*distanceToAbsorption += absCorrection;
//
//                   // const floating_t distanceToInteraction = min(
//                   //   *distancePropagated + scaCorrection,
//                   //   *distanceToAbsorption + absCorrection
//                   // );
//                   //
//                   // printf("distanceToInteraction = %f\n", distanceToInteraction);
//                   //
//                   // if (distanceToInteraction > *distancePropagated + scaCorrection) {
//                   //   *distancePropagated += hole_ice_distance_correction(
//                   //     distanceToInteraction,
//                   //     holeIceScatteringLengthFactor,
//                   //     p
//                   //   );
//                   // } else {
//                   //   *distancePropagated += scaCorrection;
//                   // }
//                   //
//                   // if (distanceToInteraction > *distanceToAbsorption + absCorrection) {
//                   //   p.bx = photonPosAndTime.x + photonDirAndWlen.x * distanceToInteraction;
//                   //   p.by = photonPosAndTime.y + photonDirAndWlen.y * distanceToInteraction;
//                   //
//                   //   *distanceToAbsorption += hole_ice_distance_correction(
//                   //     distanceToInteraction,
//                   //     holeIceAbsorptionLengthFactor,
//                   //     p
//                   //   );
//                   // } else {
//                   //   *distanceToAbsorption += absCorrection;
//                   // }




// 		  printf("HOLE ICE DEBUG:\n");
// 		  printf("  distancePropagated = %f\n", distancePropagated);
// 		  printf("  correction         = %f\n", hole_ice_distance_correction(distancePropagated, holeIceScatteringLengthFactor, p));


// #ifdef PRINTF_ENABLED
//             printf("distanceToAbsorption = %f\n", distanceToAbsorption);
// #endif

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
