#ifndef PROPAGATION_THROUGH_MEDIA_C
#define PROPAGATION_THROUGH_MEDIA_C

#include "../intersection/intersection.c"


// PROPAGATION THROUGH DIFFERENT MEDIA 2018: Layers, Cylinders
// -----------------------------------------------------------------------------

// We know how many scattering lengths (`sca_step_left`) and
// absorption lengths (`abs_lens_left`) the photon will
// travel in this step.
//
// Because the mean scattering and absorption lengths are local
// properties, i.e. depend on the ice layer or whether the photon
// is within a hole-ice cylinder, we need to convert `sca_step_left`
// and `abs_lens_left` to geometrical distances in order to determine
// where the next interaction point is, i.e. how far to propagate
// the photon in this step.

inline void apply_propagation_through_different_media(floating4_t photonPosAndTime, floating4_t photonDirAndWlen, unsigned int numberOfCylinders, __constant floating4_t *cylinderPositionsAndRadii, __constant floating_t *cylinderScatteringLengths, __constant floating_t *cylinderAbsorptionLengths, floating_t *sca_step_left, floating_t *abs_lens_left, floating_t *distancePropagated, floating_t *distanceToAbsorption)
{

  int number_of_medium_changes = 0;
  floating_t distances_to_medium_changes[MEDIUM_LAYERS] = {0.0};
  int currentPhotonLayer = min(max(findLayerForGivenZPos(photonPosAndTime.z), 0), MEDIUM_LAYERS-1);
  floating_t local_scattering_lengths[MEDIUM_LAYERS] = {getScatteringLength(currentPhotonLayer, photonDirAndWlen.w)};
  floating_t local_absorption_lengths[MEDIUM_LAYERS] = {getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w)};

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
          sqr(*sca_step_left * local_scattering_lengths[0] + cylinderPositionsAndRadii[i].w /* radius */))
      {

        // If the cylinder has a z-range check if we consider that cylinder
        // to be in range. https://github.com/fiedl/hole-ice-study/issues/34
        //
        if ((cylinderPositionsAndRadii[i].z == 0) || ((cylinderPositionsAndRadii[i].z != 0) && !(((photonPosAndTime.z < cylinderPositionsAndRadii[i].z - 0.5) && (photonPosAndTime.z + *sca_step_left * local_scattering_lengths[0] * photonDirAndWlen.z < cylinderPositionsAndRadii[i].z - 0.5)) || ((photonPosAndTime.z > cylinderPositionsAndRadii[i].z + 0.5) && (photonPosAndTime.z + *sca_step_left * local_scattering_lengths[0] * photonDirAndWlen.z > cylinderPositionsAndRadii[i].z + 0.5)))))
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
          local_scattering_lengths[number_of_medium_changes] = cylinderScatteringLengths[i];
          local_absorption_lengths[number_of_medium_changes] = cylinderAbsorptionLengths[i];
        } else if (intersection_s1(p) > 0) {
          // The photon enters the hole ice on its way.
          number_of_medium_changes += 1;
          distances_to_medium_changes[number_of_medium_changes] = intersection_s1(p);
          local_scattering_lengths[number_of_medium_changes] = cylinderScatteringLengths[i];
          local_absorption_lengths[number_of_medium_changes] = cylinderAbsorptionLengths[i];
        }
        if (intersection_s2(p) > 0) {
          // The photon leaves the hole ice on its way.
          number_of_medium_changes += 1;
          distances_to_medium_changes[number_of_medium_changes] = intersection_s2(p);
          currentPhotonLayer = min(max(findLayerForGivenZPos(photonPosAndTime.z + photonDirAndWlen.z * intersection_s2(p)), 0), MEDIUM_LAYERS-1);
          local_scattering_lengths[number_of_medium_changes] = getScatteringLength(currentPhotonLayer, photonDirAndWlen.w);
          local_absorption_lengths[number_of_medium_changes] = getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w);
        }
      }

    }
  }

  // Sort the arrays `distances_to_medium_changes`, `local_scattering_lengths` and
  // `local_absorption_lengths` by ascending distance to have the medium changes
  // in the right order.
  //
  // https://en.wikiversity.org/wiki/C_Source_Code/Sorting_array_in_ascending_and_descending_order
  //
  for (int k = 0; k < number_of_medium_changes; k++) {
    for (int l = 0; l < number_of_medium_changes; l++) {
      if (distances_to_medium_changes[l] > distances_to_medium_changes[k]) {
        floating_t tmp_distance = distances_to_medium_changes[k];
        floating_t tmp_scattering = local_scattering_lengths[k];
        floating_t tmp_absorption = local_absorption_lengths[k];
        distances_to_medium_changes[k] = distances_to_medium_changes[l];
        local_scattering_lengths[k] = local_scattering_lengths[l];
        local_absorption_lengths[k] = local_absorption_lengths[l];
        distances_to_medium_changes[l] = tmp_distance;
        local_scattering_lengths[l] = tmp_scattering;
        local_absorption_lengths[l] = tmp_absorption;
      }
    }
  }

  //printf("  before medium-changes loop:\n");
  //printf("    number_of_medium_changes = %i\n", number_of_medium_changes);
  //printf("    *sca_step_left = %f\n", *sca_step_left);
  //printf("    *abs_lens_left = %f\n", *abs_lens_left);

  // For each medium, calculate the distance in that medium
  // and modify `*distancePropagated`, `*distanceToAbsorption`,
  // `*sca_step_left` and `*abs_lens_left`.
  for (int j = 0; (j < number_of_medium_changes) && (*sca_step_left > 0); j++) {
    const floating_t max_distance_in_current_medium = distances_to_medium_changes[j+1] - distances_to_medium_changes[j];
    if (*sca_step_left * local_scattering_lengths[j] > max_distance_in_current_medium) {
      *sca_step_left -= my_divide(max_distance_in_current_medium, local_scattering_lengths[j]);
      *distancePropagated += max_distance_in_current_medium;
    } else {
      *distancePropagated += *sca_step_left * local_scattering_lengths[j];
      *sca_step_left = 0;
    }
    if (*abs_lens_left * local_absorption_lengths[j] > max_distance_in_current_medium) {
      *abs_lens_left -= my_divide(max_distance_in_current_medium, local_absorption_lengths[j]);
      *distanceToAbsorption += max_distance_in_current_medium;
    } else {
      *distanceToAbsorption += *abs_lens_left * local_absorption_lengths[j];
      *abs_lens_left = 0;
    }
    //printf("  within:\n");
    //printf("    j = %i\n", j);
    //printf("    local_scattering_length = %f\n", local_scattering_lengths[j]);
    //printf("    local_absorption_lengths = %f\n", local_absorption_lengths[j]);
    //printf("    *sca_step_left = %f\n", *sca_step_left);
    //printf("    *abs_lens_left = %f\n", *abs_lens_left);
    //printf("    *distancePropagated = %f\n", *distancePropagated);
    //printf("    *distanceToAbsorption = %f\n", *distanceToAbsorption);
  }

  // Spend the rest of the budget with the last medium properties.
  *distancePropagated += *sca_step_left * local_scattering_lengths[number_of_medium_changes];
  *distanceToAbsorption += *abs_lens_left * local_absorption_lengths[number_of_medium_changes];

  if (*distanceToAbsorption < *distancePropagated) {
    *distancePropagated = *distanceToAbsorption;
    *distanceToAbsorption = ZERO;
    *abs_lens_left = ZERO;
  } else {
    *abs_lens_left -= my_divide(*distancePropagated, local_absorption_lengths[number_of_medium_changes]);
  }

  //printf("  after:\n");
  //printf("    *distancePropagated = %f\n", *distancePropagated);
  //printf("    *distanceToAbsorption = %f\n", *distanceToAbsorption);
  //printf("    *sca_step_left = %f\n", *sca_step_left);
  //printf("    *abs_lens_left = %f\n", *abs_lens_left);
}

#endif