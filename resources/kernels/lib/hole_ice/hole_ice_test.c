#include <stdio.h>
#define PRINTF_ENABLED

#include "hole_ice_test.h"
#include "hole_ice.c"
#include "gtest/gtest.h"
#include "math.h"

inline floating_t my_sqrt(floating_t a) {return sqrt(a);}
inline floating_t sqr(floating_t a) {return a * a;}
inline floating_t my_nan() { return NAN; }
inline bool my_is_nan(floating_t a) { return (a != a); }
inline floating_t min(floating_t a, floating_t b) { return fmin(a, b); }

const floating_t desired_numeric_accuracy = 0.001;

IntersectionProblemParameters_t p = {
  0.0, 0.0, // A
  0.0, 0.0, // B
  0.0, 0.0, // M
  10.0      // r
};

HoleIceProblemParameters_t hip = {
  0.0,   // distance
  0.0,   // interaction_length_factor
  0.0,   // entry_point_ratio
  0.0,   // termination_point_ratio
  false, // starts_within_hole_ice
  0      // number_of_medium_changes (will be calculated)
};

namespace {
  floating_t extremeInteractionFactor = 0.0;

  TEST(ExtremeDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.bx = 20.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, extremeInteractionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.bx = 5.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, extremeInteractionFactor, p), -10.0, desired_numeric_accuracy);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.bx = 0.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, extremeInteractionFactor, p), -10.0, desired_numeric_accuracy);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = 15.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, extremeInteractionFactor, p), -15.0, desired_numeric_accuracy);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -15.0; p.bx = 15.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, extremeInteractionFactor, p), -25.0, desired_numeric_accuracy);
  }
}

namespace {
  floating_t interactionFactor = 0.5;

  TEST(DistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.bx = 20.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, interactionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.bx = 5.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, interactionFactor, p), -5.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.bx = 0.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, interactionFactor, p), -5.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = 30.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, interactionFactor, p), -10.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, interactionFactor, p), -6.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -35.0; p.bx = 35.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, interactionFactor, p), -20.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, interactionFactor, p), -11.0, desired_numeric_accuracy);
  }
}

namespace {
  floating_t asymInteractionFactor = 0.25;

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.bx = 20.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, asymInteractionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.bx = 5.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, asymInteractionFactor, p), -7.5, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.bx = 0.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, asymInteractionFactor, p), -7.5, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = 100.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, asymInteractionFactor, p), -30.0, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, asymInteractionFactor, p), -9.0, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -100.0; p.bx = 100.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, asymInteractionFactor, p), -60.0, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, asymInteractionFactor, p), -16.5, desired_numeric_accuracy);
  }
}

namespace {
  floating_t rtlInteractionFactor = 0.5;

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = -15.0; p.bx = -20.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, rtlInteractionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = 5.0; p.bx = -5.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, rtlInteractionFactor, p), -5.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = 15.0; p.bx = 0.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, rtlInteractionFactor, p), -5.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = -30.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, rtlInteractionFactor, p), -10.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.bx = -12.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, rtlInteractionFactor, p), -6.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = 35.0; p.bx = -35.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, rtlInteractionFactor, p), -20.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = 20.0; p.bx = -12.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, rtlInteractionFactor, p), -11.0, desired_numeric_accuracy);
  }
}

namespace {
  floating_t threeDInteractionFactor = 0.5;
  floating_t threeDScalingFactor = 2.0;

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.bx = 20.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.bx = 5.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-5.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.bx = 0.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-5.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = 30.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-10.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.bx = 12.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-6.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -35.0; p.bx = 35.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-20.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.bx = 12.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-11.0) * threeDScalingFactor, desired_numeric_accuracy);
  }
}

namespace {
  // TODO: Hole ice scenario with absorption and scattering.
  // Test apply_hole_ice_correction.

  floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
  floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
  unsigned int numberOfCylinders = 1;
  floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};

  TEST(ApplyHoleIceCorrection, ScatterBeforeHoleIce) {
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 5;
    floating_t distanceToAbsorption = 400;

    apply_hole_ice_correction(
      photonPosAndTime,
      photonDirAndWlen,
      numberOfCylinders,
      cylinderPositionsAndRadii,
      holeIceScatteringLengthFactor,
      holeIceAbsorptionLengthFactor,
      &distancePropagated,
      &distanceToAbsorption
    );

    EXPECT_NEAR(distancePropagated, 5.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, 400.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, ScatterWithinHoleIce) {
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 40;
    floating_t distanceToAbsorption = 400;

    apply_hole_ice_correction(
      photonPosAndTime,
      photonDirAndWlen,
      numberOfCylinders,
      cylinderPositionsAndRadii,
      holeIceScatteringLengthFactor,
      holeIceAbsorptionLengthFactor,
      &distancePropagated,
      &distanceToAbsorption
    );

    EXPECT_NEAR(distancePropagated, 40.0 - 0.5 * 10.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, (400.0 + 5.0 * (1 - 1.0 / 0.8)), desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, ScatterAfterHoleIce) {
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 80;
    floating_t distanceToAbsorption = 400;

    apply_hole_ice_correction(
      photonPosAndTime,
      photonDirAndWlen,
      numberOfCylinders,
      cylinderPositionsAndRadii,
      holeIceScatteringLengthFactor,
      holeIceAbsorptionLengthFactor,
      &distancePropagated,
      &distanceToAbsorption
    );

    EXPECT_NEAR(distancePropagated, (80 + 20 * (1 - 1 / 0.5)), desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, (400.0 + 20.0 * (1 - 1 / 0.8)), desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, ImmediateAbsorptionInHoleIce) {
    floating_t holeIceScatteringLengthFactor = 1.0;
    floating_t holeIceAbsorptionLengthFactor = 0.0;
    floating_t distancePropagated = 60;
    floating_t distanceToAbsorption = 400;

    apply_hole_ice_correction(
      photonPosAndTime,
      photonDirAndWlen,
      numberOfCylinders,
      cylinderPositionsAndRadii,
      holeIceScatteringLengthFactor,
      holeIceAbsorptionLengthFactor,
      &distancePropagated,
      &distanceToAbsorption
    );

    EXPECT_NEAR(distancePropagated, 60.0, desired_numeric_accuracy); // Will not be corrected, because this case is already dealt with after the hole ice code in the propagation kernel.
    EXPECT_NEAR(distanceToAbsorption, 20.0 + 10.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, PhotonStartsOnRightCylinderBoarder) {
    floating_t holeIceScatteringLengthFactor = 1.0;
    floating_t holeIceAbsorptionLengthFactor = 1.0;
    floating_t distancePropagated = 40;
    floating_t distanceToAbsorption = 400;

    // // set above:
    // floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    // floating4_t photonPosAndTime = {0.0, 10.0, 1.0, 0.0};
    photonPosAndTime.x = 30.0;

    apply_hole_ice_correction(
      photonPosAndTime,
      photonDirAndWlen,
      numberOfCylinders,
      cylinderPositionsAndRadii,
      holeIceScatteringLengthFactor,
      holeIceAbsorptionLengthFactor,
      &distancePropagated,
      &distanceToAbsorption
    );

    EXPECT_NEAR(distancePropagated, 40.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, 400.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, PhotonStartsOnLeftCylinderBoarder) {
    floating_t holeIceScatteringLengthFactor = 1.0;
    floating_t holeIceAbsorptionLengthFactor = 1.0;
    floating_t distancePropagated = 40;
    floating_t distanceToAbsorption = 400;

    // // set above:
    // floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    // floating4_t photonPosAndTime = {0.0, 10.0, 1.0, 0.0};
    photonPosAndTime.x = 10.0;

    apply_hole_ice_correction(
      photonPosAndTime,
      photonDirAndWlen,
      numberOfCylinders,
      cylinderPositionsAndRadii,
      holeIceScatteringLengthFactor,
      holeIceAbsorptionLengthFactor,
      &distancePropagated,
      &distanceToAbsorption
    );

    EXPECT_NEAR(distancePropagated, 40.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, 400.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, NanIssue14) {
    // This is a reproduction of this issue:
    // https://github.com/fiedl/hole-ice-study/issues/14
    //
    // Example dataset from the logs:
    // NAN DEBUG: scaCorrection=0.000000, absCorrection=nan, photonPosAndTime=(-255.680984,-521.281982,499.060303,.), photonDirAndWlen=(-0.352114,-0.008777,0.935916,.), cylinderPositionsAndRadii={{-256.023010,-521.281982,0.000000,0.300000}}, holeIceScatteringLengthFactor=1.000000, holeIceAbsorptionLengthFactor=1.000000, distancePropagatedBeforeCorrection=0.485262, distanceToAbsorptionBeforeCorrection=59.835110

    photonPosAndTime.x = -255.680984; photonPosAndTime.y = -521.281982; photonPosAndTime.z = 499.060303;
    photonDirAndWlen.x = -0.352114; photonDirAndWlen.y = -0.008777; photonDirAndWlen.z = 0.935916;
    numberOfCylinders = 1;
    cylinderPositionsAndRadii[0].x = -256.023010;
    cylinderPositionsAndRadii[0].y = -521.281982;
    cylinderPositionsAndRadii[0].w = 0.300000;
    floating_t holeIceScatteringLengthFactor = 1.0;
    floating_t holeIceAbsorptionLengthFactor = 1.0;
    floating_t distancePropagated = 0.485262;
    floating_t distanceToAbsorption = 59.835110;

    apply_hole_ice_correction(
      photonPosAndTime,
      photonDirAndWlen,
      numberOfCylinders,
      cylinderPositionsAndRadii,
      holeIceScatteringLengthFactor,
      holeIceAbsorptionLengthFactor,
      &distancePropagated,
      &distanceToAbsorption
    );

    EXPECT_NEAR(distancePropagated, 0.485262, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, 59.835110, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, NanIssue14InstantAbsorption) {
    photonPosAndTime.x = -255.680984; photonPosAndTime.y = -521.281982; photonPosAndTime.z = 499.060303;
    photonDirAndWlen.x = -0.352114; photonDirAndWlen.y = -0.008777; photonDirAndWlen.z = 0.935916;
    numberOfCylinders = 1;
    cylinderPositionsAndRadii[0].x = -256.023010;
    cylinderPositionsAndRadii[0].y = -521.281982;
    cylinderPositionsAndRadii[0].w = 0.300000;
    floating_t holeIceScatteringLengthFactor = 1.0;
    floating_t holeIceAbsorptionLengthFactor = 0.0;
    floating_t distancePropagated = 0.485262;
    floating_t distanceToAbsorption = 59.835110;

    const floating_t xyProjectionFactor = my_sqrt(1 - sqr(photonDirAndWlen.z));
    IntersectionProblemParameters_t intersection_problem = {
      -255.680984, -521.281982, // A
      -255.680984 + (-0.352114 * distancePropagated * xyProjectionFactor), -521.281982 + (-0.008777 * distancePropagated * xyProjectionFactor), // B
      -256.023010, -521.281982, // M
      0.300000 // r
    };
    const floating_t distance_to_first_intersection_point = intersection_s1(intersection_problem) * distancePropagated;

    apply_hole_ice_correction(
      photonPosAndTime,
      photonDirAndWlen,
      numberOfCylinders,
      cylinderPositionsAndRadii,
      holeIceScatteringLengthFactor,
      holeIceAbsorptionLengthFactor,
      &distancePropagated,
      &distanceToAbsorption
    );

    EXPECT_NEAR(distancePropagated, 0.485262, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, distance_to_first_intersection_point, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, SignIssue17) {
    // In this scenario the photon is scattered away before reaching the hole ice.
    // Absorption and scattering should not be corrected in this case.
    // https://github.com/fiedl/hole-ice-study/issues/17

    photonPosAndTime.x = -255.523010; photonPosAndTime.y = -521.281982; photonPosAndTime.z = 499.133972;
    photonDirAndWlen.x = -0.496493; photonDirAndWlen.y = 0.001049; photonDirAndWlen.z = 0.868040;
    numberOfCylinders = 1;
    cylinderPositionsAndRadii[0].x = -256.023010;
    cylinderPositionsAndRadii[0].y = -521.281982;
    cylinderPositionsAndRadii[0].w = 0.300000;
    floating_t holeIceScatteringLengthFactor = 0.1;
    floating_t holeIceAbsorptionLengthFactor = 0.1;
    floating_t distancePropagatedBeforeCorrection = 0.648981;
    floating_t distanceToAbsorptionBeforeCorrection = 31.514143;
    floating_t distancePropagated = distancePropagatedBeforeCorrection;
    floating_t distanceToAbsorption = distanceToAbsorptionBeforeCorrection;

    apply_hole_ice_correction(
      photonPosAndTime,
      photonDirAndWlen,
      numberOfCylinders,
      cylinderPositionsAndRadii,
      holeIceScatteringLengthFactor,
      holeIceAbsorptionLengthFactor,
      &distancePropagated,
      &distanceToAbsorption
    );

    EXPECT_TRUE(distancePropagated == distancePropagatedBeforeCorrection);
    EXPECT_TRUE(distanceToAbsorption == distanceToAbsorptionBeforeCorrection);
  }

}