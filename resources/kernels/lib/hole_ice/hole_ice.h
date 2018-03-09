#ifndef HOLE_ICE_H
#define HOLE_ICE_H

inline void add_hole_ice_cylinders_on_photon_path_to_medium_changes(floating4_t photonPosAndTime, floating4_t photonDirAndWlen, floating_t photonRange, unsigned int numberOfCylinders, __constant floating4_t *cylinderPositionsAndRadii, int *number_of_medium_changes, floating_t *distances_to_medium_changes, floating_t *local_scattering_lengths, floating_t *local_absorption_lengths);

#endif