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

inline floating_t apply_hole_ice_correction(floating4_t photonPosAndTime, floating4_t photonDirAndWlen, unsigned int numberOfCylinders, __constant floating4_t *cylinderPositionsAndRadii, floating_t holeIceScatteringLengthFactor, floating_t holeIceAbsorptionLengthFactor, floating_t *distancePropagated, floating_t *distanceToAbsorption)
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

		              IntersectionProblemParameters_t p = {
                    photonPosAndTime.x,
                    photonPosAndTime.y,
                    photonPosAndTime.x + photonDirAndWlen.x * *distancePropagated,
                    photonPosAndTime.y + photonDirAndWlen.y * *distancePropagated,
                    cylinderPositionsAndRadii[i].x,
                    cylinderPositionsAndRadii[i].y,
                    cylinderPositionsAndRadii[i].w // radius
                  };
                  
		              //printf("distancePropagated = %f\n", distancePropagated);
		              
                  const floating_t scaCorrection = hole_ice_distance_correction(
                    *distancePropagated,
                    holeIceScatteringLengthFactor,
                    p
                  );
                  
		              //printf("scaCorrection = %f\n", scaCorrection);
                    
                  floating_t sum = *distancePropagated;
                    
                  floating_t d;
                  d = -0.1; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -0.2; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -0.3; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -0.4; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -0.5; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -0.6; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -0.7; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -0.8; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -0.9; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.0; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.1; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.2; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.3; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.4; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.5; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.6; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.7; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.8; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -1.9; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;
                  d = -2.0; if ((scaCorrection < d) && (scaCorrection > (d - 0.1))) sum = *distancePropagated + d;

                  //distancePropagated += scaCorrection; // Seltsam. Ohne diese Zeile geht es.
                  
		              //const floating_t sum = distancePropagated + scaCorrection;
		              //printf("sum = %f\n", sum);
                  *distancePropagated = sum;
                  
		              //const floating_t foo = 0.9 * distancePropagated;  // geht
		              //distancePropagated = distancePropagated + foo;
                  
		              //distancePropagated = distancePropagated; // geht.
		              //distancePropagated = distancePropagated + 0.0 * scaCorrection; // geht nicht.

	      	        // Wenn auch nur eine Zeile mit "geht nicht" drin ist, wird keiner der prints ausgeführt.
                  //if (scaCorrection < -1.0) distancePropagated -= 1.0; // geht.

                  // distancePropagated += round(scaCorrection); // geht nicht.

                  // geht: d.h. es liegt nicht am Schreiben in die Variable, sondern
                  // irgendwo bei der Weiterverarbeitung.
                  //distancePropagated += scaCorrection; 
                  //distancePropagated = 1.0;
                                    
                  // // geht:
                  // distancePropagated *= holeIceScatteringLengthFactor;
                  // printf("holeIceScatteringLengthFactor = %f\n", holeIceScatteringLengthFactor);
                  
                  //distancePropagated += 0 * scaCorrection;
                  //distancePropagated *= holeIceScatteringLengthFactor;
                  //printf("distancePropagated = %f NEW\n", distancePropagated);

                  

                  //printf("distancePropagated = %f before resetting to 1.0\n", distancePropagated);
                  //distancePropagated = 1.0;
      		  

                  const floating_t distanceToConsiderForAbsorption = min( // Gamma
                    *distancePropagated, // after hole-ice correction
                    *distanceToAbsorption // before hole-ice correction
                  );

                  p.bx = photonPosAndTime.x + photonDirAndWlen.x * distanceToConsiderForAbsorption;
                  p.by = photonPosAndTime.y + photonDirAndWlen.y * distanceToConsiderForAbsorption;

                  const floating_t absCorrection = hole_ice_distance_correction(
                    distanceToConsiderForAbsorption,
                    holeIceAbsorptionLengthFactor,
                    p
                  );
                  
                  sum = 0;  
                  for (floating_t d = -0.0; d > -5.0; d -= 0.1) {
                    if ((absCorrection < d) && (absCorrection > (d - 0.1))) sum = *distanceToAbsorption + d;
                  }
                  //printf("  absCorrection = %f\n", absCorrection);
                  printf("  sum = %f\n", sum);
                  *distanceToAbsorption = sum;
                    
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