#ifndef ICE_LAYERS_H
#define ICE_LAYERS_H

inline void add_ice_layers_on_photon_path_to_medium_changes(floating4_t photonPosAndTime, floating4_t photonDirAndWlen, floating_t photonRange, int *number_of_medium_changes, floating_t *distances_to_medium_changes, floating_t *local_scattering_lengths, floating_t *local_absorption_lengths);

inline int photon_layer(floating_t z);

#endif