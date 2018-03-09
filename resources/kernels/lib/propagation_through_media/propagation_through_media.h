#ifndef PROPAGATION_THROUGH_MEDIA_H
#define PROPAGATION_THROUGH_MEDIA_H

inline void apply_propagation_through_different_media(floating4_t photonPosAndTime, floating4_t photonDirAndWlen, unsigned int numberOfCylinders, __constant floating4_t *cylinderPositionsAndRadii, __constant floating_t *cylinderScatteringLengths, __constant floating_t *cylinderAbsorptionLengths, floating_t *sca_step_left, floating_t *abs_lens_left, floating_t *distancePropagated, floating_t *distanceToAbsorption);

inline void sort_medium_changes_by_ascending_distance(int number_of_medium_changes, floating_t *distances_to_medium_changes, floating_t *local_scattering_lengths, floating_t *local_absorption_lengths);

#endif
