#include "hole_ice_test.h"
#include "hole_ice.c"
#include "gtest/gtest.h"
#include "math.h"

inline floating_t my_sqrt(floating_t a) {return sqrt(a);}
inline floating_t sqr(floating_t a) {return a * a;}
inline floating_t my_nan() { return NAN; }
inline bool my_is_nan(floating_t a) { return isnan(a); }

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
    p.ax = 0.0; p.bx = 20.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -5.0, 0.001);
  }

  TEST(DistanceCorrectionTest, BeginsInsideWithOneIntersectionButNoIntersectionAfterScaling) {
    p.ax = 0.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -6.0, 0.001);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithTwoIntersections) {
    p.ax = -20.0; p.bx = 20.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -10.0, 0.001);
  }

  TEST(DistanceCorrectionTest, BeginsOutsideWithTwoIntersectionsButOnlyOneIntersectionAfterScaling) {
    p.ax = -20.0; p.bx = 12.0;
    const floating_t dst = p.bx - p.ax;
    EXPECT_NEAR(hole_ice_distance_correction(dst, interactionFactor, p), -11.0, 0.001);
  }
}