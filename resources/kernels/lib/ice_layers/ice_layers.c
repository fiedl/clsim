#ifndef ICE_LAYERS_C
#define ICE_LAYERS_C

#include "ice_layers.h"

inline void add_ice_layers_on_photon_path_to_medium_changes(floating4_t photonPosAndTime, floating4_t photonDirAndWlen, floating_t photonRange, int *number_of_medium_changes, floating_t *distances_to_medium_changes, floating_t *local_scattering_lengths, floating_t *local_absorption_lengths)
{

  // The closest ice layer is special, because we need to check how far
  // it is away from the photon. After that, all photon layers are equidistant.
  //
  floating_t z_of_closest_ice_layer_boundary =
      mediumLayerBoundary(photon_layer(photonPosAndTime.z));
  if (photonDirAndWlen.z > ZERO) z_of_closest_ice_layer_boundary +=
      (floating_t)MEDIUM_LAYER_THICKNESS;

  *number_of_medium_changes += 1;
  distances_to_medium_changes[*number_of_medium_changes] =
      my_divide(z_of_closest_ice_layer_boundary - photonPosAndTime.z, photonDirAndWlen.z);
  int next_photon_layer =
      photon_layer(z_of_closest_ice_layer_boundary + photonDirAndWlen.z);
  local_scattering_lengths[*number_of_medium_changes] =
      10000.0;
      //getScatteringLength(next_photon_layer, photonDirAndWlen.w);
  local_absorption_lengths[*number_of_medium_changes] =
      getAbsorptionLength(next_photon_layer, photonDirAndWlen.w);

  // printf("ICE LAYER DEBUG\n");
  // printf("  photonPosAndTime.z = %f\n", photonPosAndTime.z);
  // printf("  photonDirAndWlen.z = %f\n", photonDirAndWlen.z);
  // printf("  z_of_closest_ice_layer_boundary = %f\n", z_of_closest_ice_layer_boundary);
  // printf("  distance to boundary = %f\n", distances_to_medium_changes[*number_of_medium_changes]);
  // printf("  current photon layer = %i\n", photon_layer(photonPosAndTime.z));
  // printf("  next_photon_layer = %i\n", next_photon_layer);

  // Now loop through the equidistant layers in range.
  //
  const floating_t max_trajectory_length_between_two_layers =
      my_divide((floating_t)MEDIUM_LAYER_THICKNESS, my_fabs(photonDirAndWlen.z));
  while (distances_to_medium_changes[*number_of_medium_changes] + max_trajectory_length_between_two_layers < photonRange)
  {
    *number_of_medium_changes += 1;
    distances_to_medium_changes[*number_of_medium_changes] =
        distances_to_medium_changes[*number_of_medium_changes - 1]
        + max_trajectory_length_between_two_layers;
    next_photon_layer = photon_layer(photonPosAndTime.z
        + (distances_to_medium_changes[*number_of_medium_changes] + 0.01) * photonDirAndWlen.z);
    local_scattering_lengths[*number_of_medium_changes] =
        10000.0;
        //getScatteringLength(next_photon_layer, photonDirAndWlen.w);
    local_absorption_lengths[*number_of_medium_changes] =
        getAbsorptionLength(next_photon_layer, photonDirAndWlen.w);
  }

  // printf("  *number_of_medium_chnages = %i\n", *number_of_medium_changes);
  // for (int i = 0; i <= *number_of_medium_changes; i++)
  // {
  //   printf("  i = %i\n", i);
  //   printf("    distance = %f\n", distances_to_medium_changes[i]);
  //   printf("    scattering length = %f\n", local_scattering_lengths[i]);
  //   printf("    absorption length = %f\n", local_absorption_lengths[i]);
  // }

}

inline int photon_layer(floating_t z)
{
  return min(max(findLayerForGivenZPos(z), 0), MEDIUM_LAYERS-1);
}

#endif