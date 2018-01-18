// #include <stdio.h>
// #define PRINTF_ENABLED

#include "intersection_test.h"
#include "intersection.c"
#include "gtest/gtest.h"
#include "math.h"

inline floating_t my_sqrt(floating_t a) {return sqrt(a);}
inline floating_t sqr(floating_t a) {return a * a;}
inline floating_t my_nan() { return NAN; }
inline bool my_is_nan(floating_t a) { return isnan(a); }

IntersectionProblemParameters_t parameters_for_tangent = {
  0.0, 0.5, // A
  5.0, 0.5, // B
  1.0, 1.0, // M
  0.5       // r
};

IntersectionProblemParameters_t parameters_for_two_intersections = {
  0.0, 1.0, // A
  5.0, 1.0, // B
  1.0, 1.0, // M
  0.5       // r
};

IntersectionProblemParameters_t parameters_for_no_intersection = {
  0.0, 0.5, // A
  5.0, 0.5, // B
  1.0, 1.0, // M
  0.2       // r
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_inside = {
  1.0, 1.0, // A
  5.0, 1.0, // B
  1.0, 1.0, // M
  1.0       // r
};

IntersectionProblemParameters_t parameters_for_two_intersections_from_right_to_left = {
  2.0, 1.0, // A
 -3.0, 1.0, // B
  1.0, 1.0, // M
  0.5       // r
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_on_border_outwards = {
  30.0, 10.0,  // A
  100.0, 10.0, // B
  20.0, 10.0,  // M
  10.0         // r
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_on_border_inwards = {
  30.0, 10.0, // A
  20.0, 10.0, // B
  20.0, 10.0, // M
  10.0        // r
};

IntersectionProblemParameters_t parameters_for_trajectory_ending_on_border = {
  20.0, 10.0, // A
  30.0, 10.0, // B
  20.0, 10.0, // M
  10.0        // r
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_outside_no_intersection_outwards = {
  15.0, 10.0, // A
  20.0, 10.0, // B
  0.0, 10.0,  // M
  10.0        // r
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_outside_no_intersection_towards = {
  20.0, 10.0, // A
  15.0, 10.0, // B
  0.0, 10.0,  // M
  10.0        // r
};

TEST(IntersectionPointsTest, TangentPoint) {
  EXPECT_EQ(intersection_x1(parameters_for_tangent), 1.0);
  EXPECT_EQ(intersection_y1(parameters_for_tangent), 0.5);
  EXPECT_EQ(intersection_x2(parameters_for_tangent), 1.0);
  EXPECT_EQ(intersection_y2(parameters_for_tangent), 0.5);
}
TEST(IntersectionPointsTest, TwoIntersectionPoints) {
  EXPECT_EQ(intersection_x1(parameters_for_two_intersections), 0.5);
  EXPECT_EQ(intersection_y1(parameters_for_two_intersections), 1.0);
  EXPECT_EQ(intersection_x2(parameters_for_two_intersections), 1.5);
  EXPECT_EQ(intersection_y2(parameters_for_two_intersections), 1.0);
}
TEST(IntersectionPointsTest, NoIntersectionPoint) {
  EXPECT_TRUE(my_is_nan(intersection_x1(parameters_for_no_intersection)));
  EXPECT_TRUE(my_is_nan(intersection_y1(parameters_for_no_intersection)));
  EXPECT_TRUE(my_is_nan(intersection_x2(parameters_for_no_intersection)));
  EXPECT_TRUE(my_is_nan(intersection_y2(parameters_for_no_intersection)));
}
TEST(IntersectionPointsTest, StartingInside) {
  EXPECT_TRUE(my_is_nan(intersection_x1(parameters_for_trajectory_starting_inside)));
  EXPECT_TRUE(my_is_nan(intersection_y1(parameters_for_trajectory_starting_inside)));
  EXPECT_EQ(intersection_x2(parameters_for_trajectory_starting_inside), 2.0);
  EXPECT_EQ(intersection_y2(parameters_for_trajectory_starting_inside), 1.0);
}

TEST(IntersectionRatioTest, TangentPoint) {
  EXPECT_EQ(intersection_s1(parameters_for_tangent), 0.2);
  EXPECT_EQ(intersection_s2(parameters_for_tangent), 0.2);
}
TEST(IntersectionRatioTest, TwoIntersectionPoints) {
  EXPECT_EQ(intersection_s1(parameters_for_two_intersections), 0.1);
  EXPECT_EQ(intersection_s2(parameters_for_two_intersections), 0.3);
}
TEST(IntersectionRatioTest, NoIntersectionPoint) {
  EXPECT_TRUE(my_is_nan(intersection_s1(parameters_for_no_intersection)));
  EXPECT_TRUE(my_is_nan(intersection_s2(parameters_for_no_intersection)));
}
TEST(IntersectionRatioTest, StartingInside) {
  EXPECT_TRUE(my_is_nan(intersection_s1(parameters_for_trajectory_starting_inside)));
  EXPECT_EQ(intersection_s2(parameters_for_trajectory_starting_inside), 0.25);
}
TEST(IntersectionRatioTest, TwoIntersectionPointsRtlDirection) {
  EXPECT_EQ(intersection_s1(parameters_for_two_intersections_from_right_to_left), 0.1);
  EXPECT_EQ(intersection_s2(parameters_for_two_intersections_from_right_to_left), 0.3);
}

TEST(TrajectoryStartsInside, TangentPoint) {
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_tangent));
}
TEST(TrajectoryStartsInside, TwoIntersectionPoints) {
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_two_intersections));
}
TEST(TrajectoryStartsInside, NoIntersectionPoint) {
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_no_intersection));
}
TEST(TrajectoryStartsInside, StartingInside) {
  EXPECT_TRUE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_inside));
}
TEST(TrajectoryStartsInside, StartingOnBorderOutwards) {
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_on_border_outwards));
}
TEST(TrajectoryStartsInside, StartingOnBorderInwards) {
  EXPECT_TRUE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_on_border_inwards));
}
TEST(TrajectoryStartsInside, StartingOutsideNoIntersectionOutwards) {
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_outside_no_intersection_outwards));
}
TEST(TrajectoryStartsInside, StartingOutsideNoIntersectionTowards) {
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_outside_no_intersection_towards));
}

TEST(TrajectoryEndsInside, EndingOnBorder) {
  // printf("TEST DEBUG:  s1 = %f\n", intersection_s1(parameters_for_trajectory_ending_on_border));
  // printf("TEST DEBUG:  s2 = %f\n", intersection_s2(parameters_for_trajectory_ending_on_border));
  // printf("TEST DEBUG:  d  = %f\n", intersection_discriminant(parameters_for_trajectory_ending_on_border));
  EXPECT_TRUE(intersecting_trajectory_ends_inside(parameters_for_trajectory_ending_on_border));
}

