#include <stdio.h>
#define PRINTF_ENABLED

#include "hole_ice_test.h"
#include "hole_ice.c"
#include "gtest/gtest.h"
#include "math.h"

inline floating_t my_sqrt(floating_t a) {return sqrt(a);}
inline floating_t sqr(floating_t a) {return a * a;}
inline floating_t my_nan() { return NAN; }
inline bool my_is_nan(floating_t a) { return isnan(a); }
inline floating_t min(floating_t a, floating_t b) { return fmin(a, b); }

IntersectionProblemParameters_t p = {
  0.0, 0.0, // A
  0.0, 0.0, // B
  0.0, 0.0, // M
  10.0      // r
};

namespace {
  floating_t extremeInteractionFactor = 0.0;

  TEST(ExtremeDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.bx = 20.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, extremeInteractionFactor, p), 0.0, 0.001);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.bx = 5.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, extremeInteractionFactor, p), -10.0, 0.001);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.bx = 0.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, extremeInteractionFactor, p), -10.0, 0.001);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = 15.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, extremeInteractionFactor, p), -15.0, 0.001);
  }

  TEST(ExtremeDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -15.0; p.bx = 15.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, extremeInteractionFactor, p), -25.0, 0.001);
  }
}

namespace {
  floating_t interactionFactor = 0.5;

  TEST(DistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.bx = 20.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), 0.0, 0.001);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.bx = 5.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -5.0, 0.001);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.bx = 0.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -5.0, 0.001);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = 30.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -10.0, 0.001);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -6.0, 0.001);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -35.0; p.bx = 35.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -20.0, 0.001);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -11.0, 0.001);
  }
}

namespace {
  floating_t asymInteractionFactor = 0.25;

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.bx = 20.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, asymInteractionFactor, p), 0.0, 0.001);
  }

  TEST(AsymDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.bx = 5.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, asymInteractionFactor, p), -7.5, 0.001);
  }

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.bx = 0.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, asymInteractionFactor, p), -7.5, 0.001);
  }

  TEST(AsymDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = 100.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, asymInteractionFactor, p), -30.0, 0.001);
  }

  TEST(AsymDistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, asymInteractionFactor, p), -9.0, 0.001);
  }

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -100.0; p.bx = 100.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, asymInteractionFactor, p), -60.0, 0.001);
  }

  TEST(AsymDistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, asymInteractionFactor, p), -16.5, 0.001);
  }
}

namespace {
  floating_t rtlInteractionFactor = 0.5;

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = -15.0; p.bx = -20.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, rtlInteractionFactor, p), 0.0, 0.001);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = 5.0; p.bx = -5.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, rtlInteractionFactor, p), -5.0, 0.001);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = 15.0; p.bx = 0.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, rtlInteractionFactor, p), -5.0, 0.001);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = -30.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, rtlInteractionFactor, p), -10.0, 0.001);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.bx = -12.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, rtlInteractionFactor, p), -6.0, 0.001);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = 35.0; p.bx = -35.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, rtlInteractionFactor, p), -20.0, 0.001);
  }

  TEST(RightToLeftDistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = 20.0; p.bx = -12.0;
    const floating_t dst = -(p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, rtlInteractionFactor, p), -11.0, 0.001);
  }
}

namespace {
  floating_t threeDInteractionFactor = 0.5;
  floating_t threeDScalingFactor = 2.0;

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithoutIntersections) {
    p.ax = 15.0; p.bx = 20.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, threeDInteractionFactor, p), 0.0, 0.001);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsInsideWithoutIntersections) {
    p.ax = -5.0; p.bx = 5.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, threeDInteractionFactor, p), (-5.0) * threeDScalingFactor, 0.001);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithOneIntersection) {
    p.ax = -15.0; p.bx = 0.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, threeDInteractionFactor, p), (-5.0) * threeDScalingFactor, 0.001);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsInsideWithOneIntersection) {
    p.ax = 0.0; p.bx = 30.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, threeDInteractionFactor, p), (-10.0) * threeDScalingFactor, 0.001);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.bx = 12.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, threeDInteractionFactor, p), (-6.0) * threeDScalingFactor, 0.001);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -35.0; p.bx = 35.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, threeDInteractionFactor, p), (-20.0) * threeDScalingFactor, 0.001);
  }

  TEST(ThreeDDistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.bx = 12.0;
    const floating_t dst = threeDScalingFactor * (p.bx - p.ax);
    EXPECT_NEAR(hole_ice_distance_correction(dst, threeDInteractionFactor, p), (-11.0) * threeDScalingFactor, 0.001);
  }
}

namespace {
  // TODO: Hole ice scenario with absorption and scattering.
  // Test apply_hole_ice_correction.
  
  floating4_t photonPosAndTime = {0.0, 10.0, 1.0, 0.0};
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

    EXPECT_NEAR(distancePropagated, 5.0, 0.001);
    EXPECT_NEAR(distanceToAbsorption, 400.0, 0.001);
  }
  
  TEST(ApplyHoleIceCorrection, ScatterWithinHoleIce) {
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
    floating_t distancePropagated = 20;
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

    EXPECT_NEAR(distancePropagated, 15.0, 0.001);
    EXPECT_NEAR(distanceToAbsorption, (400.0 + 5.0 * (0.8 - 1)), 0.001);
  }
  
  TEST(ApplyHoleIceCorrection, ScatterAfterHoleIce) {
    floating_t holeIceScatteringLengthFactor = 0.5;
    floating_t holeIceAbsorptionLengthFactor = 0.8;
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

    EXPECT_NEAR(distancePropagated, 40.0, 0.001);
    EXPECT_NEAR(distanceToAbsorption, (400.0 + 20.0 * (0.8 - 1)), 0.001);
  }
  
  TEST(ApplyHoleIceCorrection, ImmediateAbsorptionInHoleIce) {
    floating_t holeIceScatteringLengthFactor = 1.0;
    floating_t holeIceAbsorptionLengthFactor = 0.0;
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

    EXPECT_NEAR(distancePropagated, 20.0, 0.001);
    EXPECT_NEAR(distanceToAbsorption, 20.0, 0.001);
  }
}