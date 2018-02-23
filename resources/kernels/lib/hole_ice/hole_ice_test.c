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
inline floating_t dot(floating4_t a, floating4_t b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

const floating_t desired_numeric_accuracy = 0.008;

IntersectionProblemParameters_t p = {
  0.0, 0.0,             // A
  0.0, 0.0,             // M
  10.0,                 // r
  {1.0, 0.0, 0.0, 0.0}, // direction
  0.0,                  // distance
  0.0,                  // discriminant
  0.0,                  // s1
  0.0,                  // s2
};

HoleIceProblemParameters_t hip = {
  0.0,   // distance
  0.0,   // interaction_length_factor
  0.0,   // entry_point_ratio
  0.0,   // termination_point_ratio
  false  // starts_within_hole_ice
};

namespace {
  floating_t extremeInteractionFactor = 0.0;

  TEST(ExtremeDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.distance = 20.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, extremeInteractionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.distance = 5.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, extremeInteractionFactor, p), -10.0, desired_numeric_accuracy);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.distance = 0.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, extremeInteractionFactor, p), -10.0, desired_numeric_accuracy);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.distance = 15.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, extremeInteractionFactor, p), -15.0, desired_numeric_accuracy);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -15.0; p.distance = 15.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, extremeInteractionFactor, p), -25.0, desired_numeric_accuracy);
  }
}

namespace {
  floating_t interactionFactor = 0.5;

  TEST(DistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.distance = 20.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, interactionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.distance = 5.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, interactionFactor, p), -5.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.distance = 0.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, interactionFactor, p), -5.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.distance = 30.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, interactionFactor, p), -10.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.distance = 12.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, interactionFactor, p), -6.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -35.0; p.distance = 35.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, interactionFactor, p), -20.0, desired_numeric_accuracy);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.distance = 12.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, interactionFactor, p), -11.0, desired_numeric_accuracy);
  }
}

namespace {
  floating_t asymInteractionFactor = 0.25;

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.distance = 20.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, asymInteractionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.distance = 5.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, asymInteractionFactor, p), -7.5, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.distance = 0.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, asymInteractionFactor, p), -7.5, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.distance = 100.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, asymInteractionFactor, p), -30.0, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.distance = 12.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, asymInteractionFactor, p), -9.0, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -100.0; p.distance = 100.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, asymInteractionFactor, p), -60.0, desired_numeric_accuracy);
  }

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.distance = 12.0 - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, asymInteractionFactor, p), -16.5, desired_numeric_accuracy);
  }
}

namespace {
  floating_t rtlInteractionFactor = 0.5;

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = -15.0; p.distance = -(-20.0 - p.ax); p.direction.x = -1;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, rtlInteractionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = 5.0; p.distance = -(-5.0 - p.ax); p.direction.x = -1;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, rtlInteractionFactor, p), -5.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = 15.0; p.distance = -(0.0 - p.ax); p.direction.x = -1;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, rtlInteractionFactor, p), -5.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.distance = -(-30.0 - p.ax); p.direction.x = -1;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, rtlInteractionFactor, p), -10.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.distance = -(-12.0 - p.ax); p.direction.x = -1;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, rtlInteractionFactor, p), -6.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = 35.0; p.distance = -(-35.0 - p.ax); p.direction.x = -1;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, rtlInteractionFactor, p), -20.0, desired_numeric_accuracy);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = 20.0; p.distance = -(-12.0 - p.ax); p.direction.x = -1;
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(p.distance, rtlInteractionFactor, p), -11.0, desired_numeric_accuracy);
  }
}

namespace {
  floating_t threeDInteractionFactor = 0.5;
  floating_t threeDScalingFactor = 2.0;

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.distance = 20.0 - p.ax; p.direction.x = 1.0;
    const floating_t dst = threeDScalingFactor * (p.distance);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), 0.0, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.distance = 5.0 - p.ax; p.direction.x = 1.0;
    const floating_t dst = threeDScalingFactor * (p.distance);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-5.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.distance = 0.0 - p.ax; p.direction.x = 1.0;
    const floating_t dst = threeDScalingFactor * (p.distance);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-5.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.distance = 30.0 - p.ax; p.direction.x = 1.0;
    const floating_t dst = threeDScalingFactor * (p.distance);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-10.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.distance = 12.0 - p.ax; p.direction.x = 1.0;
    const floating_t dst = threeDScalingFactor * (p.distance);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-6.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -35.0; p.distance = 35.0 - p.ax; p.direction.x = 1.0;
    const floating_t dst = threeDScalingFactor * (p.distance);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-20.0) * threeDScalingFactor, desired_numeric_accuracy);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.distance = 12.0 - p.ax; p.direction.x = 1.0;
    const floating_t dst = threeDScalingFactor * (p.distance);
    EXPECT_NEAR(hole_ice_distance_correction_for_intersection_problem(dst, threeDInteractionFactor, p), (-11.0) * threeDScalingFactor, desired_numeric_accuracy);
  }
}

namespace {
  TEST(ApplyHoleIceCorrection, NoIntersectionAtAll) {
    // This is case (a) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 50.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.1;
    floating_t holeIceAbsorptionLengthFactor = 0.1;
    floating_t distancePropagatedBeforeCorrection = 100.0;
    floating_t distanceToAbsorptionBeforeCorrection = 100.0;
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

  TEST(ApplyHoleIceCorrection, StartingAtTangentPoint) {
    // This is case (b) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {20.0, 20.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {-1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 100;
    floating_t distanceToAbsorption = 100;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

  TEST(ApplyHoleIceCorrection, TangentFlyingAway) {
    // This is case (c) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 20.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {-1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 100;
    floating_t distanceToAbsorption = 100;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

  TEST(ApplyHoleIceCorrection, TangentFlyingTowards) {
    // This is case (d) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 20.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 100;
    floating_t distanceToAbsorption = 100;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

  TEST(ApplyHoleIceCorrection, StartInsideAndEndInside) {
    // This is case (e.i) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {15.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.2;
    floating_t distancePropagated = 10;
    floating_t distanceToAbsorption = 11;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, 10 * 0.5, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, 11 * 0.2, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsOutsideAndFaceAway) {
    // This is case (e.a) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {-1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.2;
    floating_t distancePropagated = 10;
    floating_t distanceToAbsorption = 10;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

  TEST(ApplyHoleIceCorrection, StartsTowardsButScatterAndAbsorbBefore) {
    // This is case (e.t) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.2;
    floating_t distancePropagated = 10;
    floating_t distanceToAbsorption = 10;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

  TEST(ApplyHoleIceCorrection, StartsInsideScatterBeforeLeaving) {
    // This is case (f.i) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {15.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.2;
    floating_t distancePropagated = 10;
    floating_t distanceToAbsorption = 100;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, 10 - 10 * 0.5, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, 100 + (10 * 0.5) * (1 - 1/0.2), desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsTowardsScatterBeforeAndAbsorbInside) {
    // This is case (f.t) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.2;
    floating_t distancePropagated = 10;
    floating_t distanceToAbsorption = 40;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

  TEST(ApplyHoleIceCorrection, StartsTowardsScatterBeforeAbsorbAfterHoleIce) {
    // This is case (g) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
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

  TEST(ApplyHoleIceCorrection, StartsInsideScatterOutsideAbsorbInside) {
    // This is case (h.i) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {15.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.2;
    floating_t distancePropagated = 100;
    floating_t distanceToAbsorption = 10;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, distancePropagatedBeforeCorrection + (1 - 1/0.5) * 15.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, distanceToAbsorptionBeforeCorrection * 0.2, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsTowardsScatterInsideButAbsorbBeforeHoleIce) {
    // This is case (h.t) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.2;
    floating_t distancePropagated = 40;
    floating_t distanceToAbsorption = 10;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, distancePropagatedBeforeCorrection - 10 * 0.5, desired_numeric_accuracy); // Will be calculated anyway, but the photon will not scatter again, because it is absorbed before that.
    EXPECT_NEAR(distanceToAbsorption, distanceToAbsorptionBeforeCorrection, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsInsideScatterAndAbsorbOutside) {
    // This is case (i.i) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {15.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.2;
    floating_t distancePropagated = 100;
    floating_t distanceToAbsorption = 200;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, distancePropagatedBeforeCorrection + (1 - 1/0.5) * 15.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, distanceToAbsorptionBeforeCorrection + (1 - 1/0.2) * 15.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsTowardsScatterAndAbsorbInside) {
    // This is case (i.t) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 40;
    floating_t distanceToAbsorption = 45;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, distancePropagatedBeforeCorrection - (1 - 0.5) * 10.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, distanceToAbsorptionBeforeCorrection + (1 - 1/0.8) * (1 - 0.5) * 10.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsTowardsScatterWithinHoleIce) {
    // This is case (j) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
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

  TEST(ApplyHoleIceCorrection, StartsTowardsScatterAfterButAbsorbBeforeHoleIce) {
    // This is case (k) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 400;
    floating_t distanceToAbsorption = 10;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, distancePropagatedBeforeCorrection + (1 - 1/0.5) * 20.0, desired_numeric_accuracy); // Will be calculated anyway, but the photon will not scatter again, because it is absorbed before that.
    EXPECT_NEAR(distanceToAbsorption, distanceToAbsorptionBeforeCorrection, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsTowardsScatterAfterButAbsorbInsideHoleIce) {
    // This is case (l) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 400;
    floating_t distanceToAbsorption = 40;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, distancePropagatedBeforeCorrection + (1 - 1/0.5) * 20.0, desired_numeric_accuracy); // Will be calculated anyway, but the photon will not scatter again, because it is absorbed before that.
    EXPECT_NEAR(distanceToAbsorption, distanceToAbsorptionBeforeCorrection - (1 - 0.8) * 10.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsTowardsScatterAndAbsorbAfterHoleIce) {
    // This is case (m) from https://github.com/fiedl/hole-ice-study/issues/20.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
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
    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
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
    floating4_t photonPosAndTime = {30.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 1.0;
    floating_t holeIceAbsorptionLengthFactor = 1.0;
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

    EXPECT_NEAR(distancePropagated, 40.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, 400.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, PhotonStartsOnLeftCylinderBoarder) {
    floating4_t photonPosAndTime = {10.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 1.0;
    floating_t holeIceAbsorptionLengthFactor = 1.0;
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

    EXPECT_NEAR(distancePropagated, 40.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, 400.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, NanIssue14) {
    // This is a reproduction of this issue:
    // https://github.com/fiedl/hole-ice-study/issues/14
    //
    // Example dataset from the logs:
    // NAN DEBUG: scaCorrection=0.000000, absCorrection=nan, photonPosAndTime=(-255.680984,-521.281982,499.060303,.), photonDirAndWlen=(-0.352114,-0.008777,0.935916,.), cylinderPositionsAndRadii={{-256.023010,-521.281982,0.000000,0.300000}}, holeIceScatteringLengthFactor=1.000000, holeIceAbsorptionLengthFactor=1.000000, distancePropagatedBeforeCorrection=0.485262, distanceToAbsorptionBeforeCorrection=59.835110

    floating4_t photonPosAndTime = {-255.680984, -521.281982, 499.060303, 0.0};
    floating4_t photonDirAndWlen = {-0.352114, -0.008777, 0.935916, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{-256.023010, -521.281982, 0.0, 0.300000}};
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
    floating4_t photonPosAndTime = {-255.680984, -521.281982, 499.060303, 0.0};
    floating4_t photonDirAndWlen = {-0.352114, -0.008777, 0.935916, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{-256.023010, -521.281982, 0.0, 0.300000}};
    floating_t holeIceScatteringLengthFactor = 1.0;
    floating_t holeIceAbsorptionLengthFactor = 0.0;
    floating_t distancePropagated = 0.485262;
    floating_t distanceToAbsorption = 59.835110;

    IntersectionProblemParameters_t intersection_problem = {
      -255.680984, -521.281982, // A
      -256.023010, -521.281982, // M
      0.300000,                 // r
      photonDirAndWlen,         // direction
      distancePropagated,       // distance
      0,0,0                     // (will be filled)
    };
    calculate_intersections(&intersection_problem);
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

  TEST(ApplyHoleIceCorrection, NewNanIssue14) {
    // This is a reproduction of this issue:
    // https://github.com/fiedl/hole-ice-study/issues/14#issuecomment-363432459
    //
    // Example dataset from the logs:
    // NAN DEBUG: scaCorrection=0.000000, absCorrection=nan, photonPosAndTime=(-256.090363,-523.684937,503.091675,.), photonDirAndWlen=(0.841316,0.526975,-0.120350,.), cylinderPositionsAndRadii={{-256.023010,-521.281982,0.000000,0.300000}}, holeIceScatteringLengthFactor=0.100000, holeIceAbsorptionLengthFactor=0.100000, distancePropagatedBeforeCorrection=2.708774, distanceToAbsorptionBeforeCorrection=63.779461

    floating4_t photonPosAndTime = {-256.090363, -523.684937, 503.091675, 0.0};
    floating4_t photonDirAndWlen = {0.841316, 0.526975, -0.120350, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{-256.023010, -521.281982, 0.0, 0.300000}};
    floating_t holeIceScatteringLengthFactor = 0.1;
    floating_t holeIceAbsorptionLengthFactor = 0.1;
    floating_t distancePropagatedBeforeCorrection = 2.708774;
    floating_t distanceToAbsorptionBeforeCorrection = 63.779461;
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

  TEST(ApplyHoleIceCorrection, AnotherNanIssue14) {
    // This is a reproduction of this issue:
    // https://github.com/fiedl/hole-ice-study/issues/14#issuecomment-363432459
    //
    // Example dataset from the logs:
    // NAN DEBUG: scaCorrection=0.000000, absCorrection=nan, photonPosAndTime=(-255.095627,-521.281860,500.002075,.), photonDirAndWlen=(-0.898932,0.426471,0.100219,.), cylinderPositionsAndRadii={{-256.023010,-521.281982,0.000000,0.300000}}, holeIceScatteringLengthFactor=0.100000, holeIceAbsorptionLengthFactor=0.100000, distancePropagatedBeforeCorrection=0.636836, distanceToAbsorptionBeforeCorrection=155.735428

    floating4_t photonPosAndTime = {-255.095627, -521.281860, 500.002075, 0.0};
    floating4_t photonDirAndWlen = {-0.898932, 0.426471, 0.100219, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{-256.023010, -521.281982, 0.0, 0.300000}};
    floating_t holeIceScatteringLengthFactor = 0.1;
    floating_t holeIceAbsorptionLengthFactor = 0.1;
    floating_t distancePropagatedBeforeCorrection = 0.636836;
    floating_t distanceToAbsorptionBeforeCorrection = 155.735428;
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

  TEST(ApplyHoleIceCorrection, SignIssue17) {
    // In this scenario the photon is scattered away before reaching the hole ice.
    // Absorption and scattering should not be corrected in this case.
    // https://github.com/fiedl/hole-ice-study/issues/17

    floating4_t photonPosAndTime = {-255.523010, -521.281982, 499.133972, 0.0};
    floating4_t photonDirAndWlen = {-0.496493, 0.001049, 0.868040, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{-256.023010, -521.281982, 0.0, 0.300000}};
    floating_t holeIceScatteringLengthFactor = 0.1;
    floating_t holeIceAbsorptionLengthFactor = 0.1;
    floating_t distancePropagatedBeforeCorrection = 0.348981;
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

  TEST(ApplyHoleIceCorrection, StartsTowardsScattersInsideScatteringFactorGreaterOne) {
    // In this scenario, the photon would scatter within the hole ice.
    //
    // But the scattering factor is greater than one, i.e. the scattering propability
    // is smaller within the hole ice, resulting in a hole-ice correction that makes
    // the photon scatter after leaving the hole ice.
    //
    // This is "Case 3b" from
    // https://github.com/fiedl/hole-ice-study/issues/2.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 100.0;
    floating_t holeIceAbsorptionLengthFactor = 0.5;
    floating_t distancePropagated = 40.0;
    floating_t distanceToAbsorption = 400.0;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, distancePropagatedBeforeCorrection + (1 - 1/100.0) * 20.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, distanceToAbsorptionBeforeCorrection + (1 - 1/0.5) * 20.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsTowardsScattersInsideAbsorptionFactorGreaterOne) {
    // In this scenario, the photon scatters within the hole ice and would also
    // be absorbed within the hole ice.
    //
    // But the absorption factor is greater than one, i.e. the absorption propability
    // is smaller within the hole ice, resulting in a hole-ice correction that makes
    // the photon be absorbed after leaving the hole ice.
    //
    // https://github.com/fiedl/hole-ice-study/issues/2.

    floating4_t photonPosAndTime = {-20.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 100.0;
    floating_t distancePropagated = 40.0;
    floating_t distanceToAbsorption = 45.0;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, distancePropagatedBeforeCorrection - 10.0 * 0.5 , desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, distanceToAbsorptionBeforeCorrection + (1 - 1/100.0) * (10.0 * 0.5), desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, StartsInsideScattersInsideScatteringFactorGreaterOne) {
    // In this scenario, the photon starts within the hole ice and
    // would scatter within the hole ice.
    //
    // But the scattering factor is greater than one, i.e. the scattering propability
    // is smaller within the hole ice, resulting in a hole-ice correction that makes
    // the photon scatter after leaving the hole ice.
    //
    // This is "Case 2b" from
    // https://github.com/fiedl/hole-ice-study/issues/2.

    floating4_t photonPosAndTime = {15.0, 10.0, 1.0, 0.0};
    floating4_t photonDirAndWlen = {1.0,  0.0, 0.0, 700e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{20.0, 10.0, 0.0, 10.0}};
    floating_t holeIceScatteringLengthFactor = 100.0;
    floating_t holeIceAbsorptionLengthFactor = 0.5;
    floating_t distancePropagated = 5.0;
    floating_t distanceToAbsorption = 400.0;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, distancePropagatedBeforeCorrection + (1 - 1/100.0) * 15.0, desired_numeric_accuracy);
    EXPECT_NEAR(distanceToAbsorption, distanceToAbsorptionBeforeCorrection + (1 - 1/0.5) * 15.0, desired_numeric_accuracy);
  }

  TEST(ApplyHoleIceCorrection, GeoGebraTestIssue25) {
    // https://github.com/fiedl/hole-ice-study/issues/25

    floating4_t photonPosAndTime = {10.0, 15.0, 1.2, 0.0};
    floating4_t photonDirAndWlen = {(8.12 - 10.0) / 17.9, (28.91 - 15.0) / 17.9, (12.3 - 1.2) / 17.9, 73542800e-9};
    unsigned int numberOfCylinders = 1;
    floating4_t cylinderPositionsAndRadii[] = {{10.2, 20.5, 0.0, 2.76}};
    floating_t holeIceScatteringLengthFactor = 0.0;
    floating_t holeIceAbsorptionLengthFactor = 0.0;
    floating_t distancePropagated = 17.9;
    floating_t distanceToAbsorption = 17.9;
    floating_t distancePropagatedBeforeCorrection = distancePropagated;
    floating_t distanceToAbsorptionBeforeCorrection = distanceToAbsorption;

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

    EXPECT_NEAR(distancePropagated, 3.6, desired_numeric_accuracy);
  }
}