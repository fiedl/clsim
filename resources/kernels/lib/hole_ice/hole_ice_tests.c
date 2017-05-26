#include <stdio.h>
#include <stdbool.h>
#include "math.h"
#include "hole_ice_tests.h"
#include "hole_ice.c"

inline floating_t my_sqrt(floating_t a) {return sqrt(a);}
inline floating_t sqr(floating_t a) {return a * a;}
inline floating_t my_nan() { return NAN; }
inline bool my_is_nan(floating_t a) { return isnan(a); }


int test_trajectory_outside_the_circle()
{
  IntersectionProblemParameters_t p = {
    -1.0, -1.0, // A
     4.0, -1.0, // B
     1.0,  1.0, // M
     1.0        // r
  };
  floating_t trajectory_length_from_scattering = 10.0;
  floating_t trajectory_length_from_absorption = 100.0;
  floating_t number_of_scattering_lengths = 1.0;
  floating_t number_of_absorption_lengths = 10.0;
  floating_t currentScaLen = trajectory_length_from_scattering / number_of_scattering_lengths;
  floating_t currentAbsLen = trajectory_length_from_absorption / number_of_absorption_lengths;
  floating_t holeIceScatteringLengthFactor = 0.5;
  floating_t holeIceAbsorptionLengthFactor = 0.5;

  floating_t expected_trajectory_length_from_scattering = 10.0;

  perform_hole_ice_correction_for(
    &trajectory_length_from_scattering,
    &trajectory_length_from_absorption,
    &number_of_scattering_lengths,
    &number_of_absorption_lengths,
    currentScaLen,
    currentAbsLen,
    holeIceScatteringLengthFactor,
    holeIceAbsorptionLengthFactor,
    p
  );

  if (trajectory_length_from_scattering == expected_trajectory_length_from_scattering) return 0; else return 1;
}

// In this test, the photon travels right through the circle.
// The trajectory within the circle is 2 units long. The scattering
// and absorption lengths are scaled with 0.5, i.e. are half the
// usual length. Therefore, corresponding to the increased probability
// to being scattered or absorbed within the circle, the propagated
// distance is reduced.
//
int test_trajectory_through_the_circle()
{
  IntersectionProblemParameters_t p = {
    -1.0,  1.0, // A
     4.0,  1.0, // B
     1.0,  1.0, // M
     1.0        // r
  };
  floating_t trajectory_length_from_scattering = 10.0;
  floating_t trajectory_length_from_absorption = 100.0;
  floating_t number_of_scattering_lengths = 1.0;
  floating_t number_of_absorption_lengths = 10.0;
  floating_t currentScaLen = trajectory_length_from_scattering / number_of_scattering_lengths;
  floating_t currentAbsLen = trajectory_length_from_absorption / number_of_absorption_lengths;
  floating_t holeIceScatteringLengthFactor = 0.5;
  floating_t holeIceAbsorptionLengthFactor = 0.5;

  floating_t expected_trajectory_length_from_scattering = 10.0 - 2.0 * (1.0 / 0.5 - 1.0);

  perform_hole_ice_correction_for(
    &trajectory_length_from_scattering,
    &trajectory_length_from_absorption,
    &number_of_scattering_lengths,
    &number_of_absorption_lengths,
    currentScaLen,
    currentAbsLen,
    holeIceScatteringLengthFactor,
    holeIceAbsorptionLengthFactor,
    p
  );

  if (trajectory_length_from_scattering == expected_trajectory_length_from_scattering) return 0; else return 1;
}

// In this test, the photon only travels within the circle.
// The whole trajectory is scaled therefore.
//
int test_trajectory_inside_the_circle()
{
  IntersectionProblemParameters_t p = {
     0.5,  1.0, // A
     1.5,  1.0, // B
     1.0,  1.0, // M
     1.0        // r
  };
  floating_t trajectory_length_from_scattering = 10.0;
  floating_t trajectory_length_from_absorption = 100.0;
  floating_t number_of_scattering_lengths = 1.0;
  floating_t number_of_absorption_lengths = 10.0;
  floating_t currentScaLen = trajectory_length_from_scattering / number_of_scattering_lengths;
  floating_t currentAbsLen = trajectory_length_from_absorption / number_of_absorption_lengths;
  floating_t holeIceScatteringLengthFactor = 0.5;
  floating_t holeIceAbsorptionLengthFactor = 0.5;

  floating_t expected_trajectory_length_from_scattering = 10.0 - holeIceScatteringLengthFactor * 1.0;

  perform_hole_ice_correction_for(
    &trajectory_length_from_scattering,
    &trajectory_length_from_absorption,
    &number_of_scattering_lengths,
    &number_of_absorption_lengths,
    currentScaLen,
    currentAbsLen,
    holeIceScatteringLengthFactor,
    holeIceAbsorptionLengthFactor,
    p
  );

  if (trajectory_length_from_scattering == expected_trajectory_length_from_scattering) return 0; else return 1;
}


int main()
{
  int numberOfFailedTests =
    test_trajectory_outside_the_circle() +
    test_trajectory_through_the_circle() +
    test_trajectory_inside_the_circle() +
    0;

  printf("Tests complete. %i tests have failed.\n", numberOfFailedTests);
  return numberOfFailedTests;
}

