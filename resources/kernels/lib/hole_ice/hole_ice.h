#ifndef HOLE_ICE_H
#define HOLE_ICE_H

typedef struct HoleIceProblemParameters {
  floating_t distance;
  floating_t interaction_length_factor;
  floating_t entry_point_ratio;
  floating_t termination_point_ratio;
  bool starts_within_hole_ice;
  unsigned int number_of_medium_changes; // will be calculated
} HoleIceProblemParameters_t;

#endif