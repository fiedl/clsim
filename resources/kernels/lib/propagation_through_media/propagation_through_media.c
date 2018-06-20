#ifndef PROPAGATION_THROUGH_MEDIA_C
#define PROPAGATION_THROUGH_MEDIA_C

#include "propagation_through_media.h"
#include "../ice_layers/ice_layers.c"
#ifdef HOLE_ICE
  #include "../hole_ice/hole_ice.c"
#endif


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

inline void apply_propagation_through_different_media(
  floating4_t photonPosAndTime, floating4_t photonDirAndWlen,
  #ifdef HOLE_ICE
    unsigned int numberOfCylinders, __constant floating4_t *cylinderPositionsAndRadii,
    __constant floating_t *cylinderScatteringLengths, __constant floating_t *cylinderAbsorptionLengths,
  #endif
  floating_t *distances_to_medium_changes, floating_t *local_scattering_lengths, floating_t *local_absorption_lengths,
  floating_t *sca_step_left, floating_t *abs_lens_left,
  floating_t *distancePropagated, floating_t *distanceToAbsorption)
{

  printf("Hole-ice code: THIS IS WORK-IN-PROGRESS. Do not use this, yet! See https://github.com/fiedl/hole-ice-study.\n");

  //clock_t t0 = clock();

  int number_of_medium_changes = 0;
  //clock_t t01 = clock();
  distances_to_medium_changes[0] = 0.0;
  //clock_t t02 = clock();
  int currentPhotonLayer = min(max(findLayerForGivenZPos(photonPosAndTime.z), 0), MEDIUM_LAYERS-1);
  //clock_t t03 = clock();
  local_scattering_lengths[0] = getScatteringLength(currentPhotonLayer, photonDirAndWlen.w);
  //clock_t t04 = clock();
  local_absorption_lengths[0] = getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w);
  //clock_t t05 = clock();

  // int number_of_medium_changes = 1;
  // floating_t distances_to_medium_changes[MEDIUM_LAYERS] = {0.0, 1.0};
  // int currentPhotonLayer = min(max(findLayerForGivenZPos(photonPosAndTime.z), 0), MEDIUM_LAYERS-1);
  // floating_t local_scattering_lengths[MEDIUM_LAYERS] = {getScatteringLength(currentPhotonLayer, photonDirAndWlen.w), 100.0};
  // floating_t local_absorption_lengths[MEDIUM_LAYERS] = {getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w), 100.0};

  //floating_t local_scattering_lengths[MEDIUM_LAYERS] = {1000.0};
  //floating_t local_absorption_lengths[MEDIUM_LAYERS] = {1000.0};


  // To check which medium boundaries are in range, we need to estimate
  // how far the photon can travel in this step.
  //
  const floating_t photonRange = *sca_step_left * local_scattering_lengths[0];

  //printf("MEDIUM CHANGES 2018 DEBUG:\n");
  //printf("  photonRange = %f\n", photonRange);
  //printf("  distancePropagated = %f\n", *distancePropagated);
  //printf("  distanceToAbsorption = %f\n", *distanceToAbsorption);
  //printf("  sca_step_left = %f\n", *sca_step_left);
  //printf("  abs_lens_left = %f\n", *abs_lens_left);

  //clock_t t1 = clock();
  add_ice_layers_on_photon_path_to_medium_changes(
    photonPosAndTime,
    photonDirAndWlen,
    photonRange,

    // These values will be updates within this function:
    &number_of_medium_changes,
    distances_to_medium_changes,
    local_scattering_lengths,
    local_absorption_lengths
  );

  //clock_t t2 = clock();
  #ifdef HOLE_ICE
    add_hole_ice_cylinders_on_photon_path_to_medium_changes(
      photonPosAndTime,
      photonDirAndWlen,
      photonRange,
      numberOfCylinders,
      cylinderPositionsAndRadii,

      // These values will be updates within this function:
      &number_of_medium_changes,
      distances_to_medium_changes,
      local_scattering_lengths,
      local_absorption_lengths
    );
  #endif

  // number_of_medium_changes = 1;
  // distances_to_medium_changes[0] = 0.0;
  // distances_to_medium_changes[1] = 1.0;
  // currentPhotonLayer = min(max(findLayerForGivenZPos(photonPosAndTime.z), 0), MEDIUM_LAYERS-1);
  // local_scattering_lengths[0] = getScatteringLength(currentPhotonLayer, photonDirAndWlen.w);
  // local_scattering_lengths[1] = 100.0;
  // local_absorption_lengths[0] = getAbsorptionLength(currentPhotonLayer, photonDirAndWlen.w);
  // local_absorption_lengths[1] = 100.0;

  //clock_t t3 = clock();
  sort_medium_changes_by_ascending_distance(
    number_of_medium_changes,

    // These values will be updates within this function:
    distances_to_medium_changes,
    local_scattering_lengths,
    local_absorption_lengths
  );

  //clock_t t4 = clock();
  loop_over_media_and_calculate_geometrical_distances_up_to_the_next_scattering_point(
    number_of_medium_changes,
    distances_to_medium_changes,
    local_scattering_lengths,
    local_absorption_lengths,

    // These values will be updates within this function:
    sca_step_left,
    abs_lens_left,
    distancePropagated,
    distanceToAbsorption
  );
  //clock_t t5 = clock();

  //printf("  after:\n");
  //printf("    distancePropagated = %f\n", *distancePropagated);
  //printf("    distanceToAbsorption = %f\n", *distanceToAbsorption);
  //printf("    sca_step_left = %f\n", *sca_step_left);
  //printf("    abs_lens_left = %f\n", *abs_lens_left);
  //printf("    number_of_medium_chnages = %i\n", number_of_medium_changes);
  //for (int i = 0; i <= number_of_medium_changes; i++)
  //{
  //  printf("    i = %i\n", i);
  //  printf("      distance = %f\n", distances_to_medium_changes[i]);
  //  printf("      scattering length = %f\n", local_scattering_lengths[i]);
  //  printf("      absorption length = %f\n", local_absorption_lengths[i]);
  //}

  //printf("  after:\n");
  //printf("    *distancePropagated = %f\n", *distancePropagated);
  //printf("    *distanceToAbsorption = %f\n", *distanceToAbsorption);
  //printf("    *sca_step_left = %f\n", *sca_step_left);
  //printf("    *abs_lens_left = %f\n", *abs_lens_left);

  //printf("PROFILING start %lu\n", t1 - t0);
  //printf("PROFILING start_01 %lu\n", t01 - t0);
  //printf("PROFILING start_02 %lu\n", t02 - t01);
  //printf("PROFILING start_03 %lu\n", t03 - t02);
  //printf("PROFILING start_04 %lu\n", t04 - t03);
  //printf("PROFILING start_05 %lu\n", t05 - t04);
  //printf("PROFILING start_1 %lu\n", t1 - t05);
  //printf("PROFILING add_ice_layers %lu\n", t2 - t1);
  //printf("PROFILING add_hole_ice %lu\n", t3 - t2);
  //printf("PROFILING sort %lu\n", t4 - t3);
  //printf("PROFILING loop_over_media %lu\n", t5 - t4);
  //printf("PROFILING total %lu\n", t5 - t0);

}

inline void sort_medium_changes_by_ascending_distance(int number_of_medium_changes, floating_t *distances_to_medium_changes, floating_t *local_scattering_lengths, floating_t *local_absorption_lengths)
{
  // Sort the arrays `distances_to_medium_changes`, `local_scattering_lengths` and
  // `local_absorption_lengths` by ascending distance to have the medium changes
  // in the right order.
  //
  // https://en.wikiversity.org/wiki/C_Source_Code/Sorting_array_in_ascending_and_descending_order
  //
  for (int k = 0; k <= number_of_medium_changes; k++) {
    for (int l = 0; l <= number_of_medium_changes; l++) {
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
}

inline void loop_over_media_and_calculate_geometrical_distances_up_to_the_next_scattering_point(int number_of_medium_changes, floating_t *distances_to_medium_changes, floating_t *local_scattering_lengths, floating_t *local_absorption_lengths, floating_t *sca_step_left, floating_t *abs_lens_left, floating_t *distancePropagated, floating_t *distanceToAbsorption)
{
  // We know how many scattering lengths (`sca_step_left`) and how many
  // absorption lengths (`abs_lens_left`) we may spend when propagating
  // through the different media.
  //
  // Convert these into the geometrical distances `distancePropagated` (scattering)
  // and `distanceToAbsorption` (absorption) and decrease `sca_step_left` and
  // `abs_lens_left` accordingly.
  //
  // Abort when the next scattering point is reached, i.e. `sca_step_left == 0`.
  // At this point, `abs_lens_left` may still be greater than zero, because
  // the photon may be scattered several times until it is absorbed.
  //
  for (int j = 0; (j < number_of_medium_changes) && (*sca_step_left > 0); j++) {
    floating_t max_distance_in_current_medium = distances_to_medium_changes[j+1] - distances_to_medium_changes[j];

    //printf("  j = %i\n", j);
    //printf("    sca: max_distance_in_current_medium = %f\n", max_distance_in_current_medium);

    if (*sca_step_left * local_scattering_lengths[j] > max_distance_in_current_medium) {
      //printf("    The photon scatters after leaving this medium.\n");
      // The photon scatters after leaving this medium.
      *sca_step_left -= my_divide(max_distance_in_current_medium, local_scattering_lengths[j]);
      *distancePropagated += max_distance_in_current_medium;
    } else {
      //printf("    The photon scatters within this medium.\n");
      // The photon scatters within this medium.
      max_distance_in_current_medium = *sca_step_left * local_scattering_lengths[j];
      *distancePropagated += max_distance_in_current_medium;
      *sca_step_left = 0;
    }

    //printf("    abs: max_distance_in_current_medium = %f\n", max_distance_in_current_medium);
    if (*abs_lens_left * local_absorption_lengths[j] > max_distance_in_current_medium) {
      //printf("    The photon is absorbed after leaving this medium.\n");
      // The photon is absorbed after leaving this medium.
      *abs_lens_left -= my_divide(max_distance_in_current_medium, local_absorption_lengths[j]);
      *distanceToAbsorption += max_distance_in_current_medium;
    } else {
      //printf("    The photon is absorbed within this medium.\n");
      // The photon is absorbed within this medium.
      *distanceToAbsorption += *abs_lens_left * local_absorption_lengths[j];
      *abs_lens_left = 0;
    }
  }

  // Spend the rest of the budget with the last medium properties.
  if (*sca_step_left > 0) {
    *distancePropagated += *sca_step_left * local_scattering_lengths[number_of_medium_changes];
    *distanceToAbsorption += *abs_lens_left * local_absorption_lengths[number_of_medium_changes];
    *abs_lens_left -= my_divide(*distancePropagated, local_absorption_lengths[number_of_medium_changes]);
  }

  // If the photon is absorbed, only propagate up to the absorption point.
  if (*distanceToAbsorption < *distancePropagated) {
    *distancePropagated = *distanceToAbsorption;
    *distanceToAbsorption = ZERO;
    *abs_lens_left = ZERO;
  }
}

#endif
