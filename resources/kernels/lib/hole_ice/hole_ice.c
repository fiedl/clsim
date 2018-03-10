#ifndef HOLE_ICE_C
#define HOLE_ICE_C

#include "hole_ice.h"
#include "../intersection/intersection.c"

inline void add_hole_ice_cylinders_on_photon_path_to_medium_changes(floating4_t photonPosAndTime, floating4_t photonDirAndWlen, floating_t photonRange, unsigned int numberOfCylinders, __constant floating4_t *cylinderPositionsAndRadii, int *number_of_medium_changes, floating_t *distances_to_medium_changes, floating_t *local_scattering_lengths, floating_t *local_absorption_lengths)
{
  // Find out which cylinders are in range in a separate loop
  // in order to improve parallelism and thereby performance.
  //
  // See: https://github.com/fiedl/hole-ice-study/issues/30
  //
  #ifdef NUMBER_OF_CYLINDERS
    // When running this on OpenCL, defining arrays using a constant
    // as array size is not possible. Therefore, we need to use a
    // pre-processor makro here.
    //
    // See: https://github.com/fiedl/hole-ice-study/issues/38
    //
    int indices_of_cylinders_in_range[NUMBER_OF_CYLINDERS];
  #else
    int indices_of_cylinders_in_range[numberOfCylinders];
  #endif
  {
    unsigned int j = 0;
    for (unsigned int i = 0; i < numberOfCylinders; i++) {
      indices_of_cylinders_in_range[i] = -1;
    }
    for (unsigned int i = 0; i < numberOfCylinders; i++) {
      if (sqr(photonPosAndTime.x - cylinderPositionsAndRadii[i].x) +
          sqr(photonPosAndTime.y - cylinderPositionsAndRadii[i].y) <=
          sqr(photonRange + cylinderPositionsAndRadii[i].w /* radius */))
      {

        // If the cylinder has a z-range check if we consider that cylinder
        // to be in range. https://github.com/fiedl/hole-ice-study/issues/34
        //
        if ((cylinderPositionsAndRadii[i].z == 0) || ((cylinderPositionsAndRadii[i].z != 0) && !(((photonPosAndTime.z < cylinderPositionsAndRadii[i].z - 0.5) && (photonPosAndTime.z + photonRange * photonDirAndWlen.z < cylinderPositionsAndRadii[i].z - 0.5)) || ((photonPosAndTime.z > cylinderPositionsAndRadii[i].z + 0.5) && (photonPosAndTime.z + photonRange * photonDirAndWlen.z > cylinderPositionsAndRadii[i].z + 0.5)))))
        {
          indices_of_cylinders_in_range[j] = i;
          j += 1;
        }
      }
    }
  }

  // Now loop over all cylinders in range and calculate corrections
  // for `*distancePropagated` and `*distanceToAbsorption`.
  //
  for (unsigned int j = 0; j < numberOfCylinders; j++) {
    const int i = indices_of_cylinders_in_range[j];
    if (i == -1) {
      break;
    } else {

      IntersectionProblemParameters_t p = {

        // Input values
        photonPosAndTime.x,
        photonPosAndTime.y,
        cylinderPositionsAndRadii[i].x,
        cylinderPositionsAndRadii[i].y,
        cylinderPositionsAndRadii[i].w, // radius
        photonDirAndWlen,
        1.0, // distance used to calculate s1 and s2 relative to

        // Output values (will be calculated)
        0, // discriminant
        0, // s1
        0  // s2

      };

      calculate_intersections(&p);

      //printf("  intersection:\n");
      //printf("    cylinder: i = %i\n", i);
      //printf("    intersection_s1 = %f\n", intersection_s1(p));
      //printf("    intersection_s2 = %f\n", intersection_s2(p));

      if (intersection_discriminant(p) > 0) {
        if ((intersection_s1(p) <= 0) && (intersection_s2(p) >= 0)) {
          // The photon is already within the hole ice.
          local_scattering_lengths[0] = cylinderScatteringLengths[i];
          local_absorption_lengths[0] = cylinderAbsorptionLengths[i];
        } else if (intersection_s1(p) > 0) {
          // The photon enters the hole ice on its way.
          *number_of_medium_changes += 1;
          distances_to_medium_changes[*number_of_medium_changes] = intersection_s1(p);
          local_scattering_lengths[*number_of_medium_changes] = cylinderScatteringLengths[i];
          local_absorption_lengths[*number_of_medium_changes] = cylinderAbsorptionLengths[i];
        }
        if (intersection_s2(p) > 0) {
          // The photon leaves the hole ice on its way.
          *number_of_medium_changes += 1;
          distances_to_medium_changes[*number_of_medium_changes] = intersection_s2(p);
          int currentPhotonLayer = min(max(findLayerForGivenZPos(photonPosAndTime.z + photonDirAndWlen.z * intersection_s2(p)), 0), MEDIUM_LAYERS-1);
          local_scattering_lengths[*number_of_medium_changes] = getScatteringLength(currentPhotonLayer, photonDirAndWlen.w);
          local_absorption_lengths[*number_of_medium_changes] = getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w);
        }
      }

    }
  }

}

#endif
