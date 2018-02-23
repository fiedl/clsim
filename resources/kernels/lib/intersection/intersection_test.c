// #include <stdio.h>
// #define PRINTF_ENABLED

#include "intersection_test.h"
#include "intersection.c"
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

IntersectionProblemParameters_t parameters_for_tangent = {
  0.0, 0.5,              // A
  1.0, 1.0,              // M
  0.5,                   // r
  {1.0, 0.0, 0.0, 0.0},  // direction
  5.0                    // distance
};

IntersectionProblemParameters_t parameters_for_two_intersections = {
  0.0, 1.0,              // A
  1.0, 1.0,              // M
  0.5,                   // r
  {1.0, 0.0, 0.0, 0.0},  // direction
  5.0                    // distance
};

IntersectionProblemParameters_t parameters_for_no_intersection = {
  0.0, 0.5,              // A
  1.0, 1.0,              // M
  0.2,                   // r
  {1.0, 0.0, 0.0, 0.0},  // direction
  5.0                    // distance
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_inside = {
  1.0, 1.0,              // A
  1.0, 1.0,              // M
  1.0,                   // r
  {1.0, 0.0, 0.0, 0.0},  // direction
  4.0                    // distance
};

IntersectionProblemParameters_t parameters_for_two_intersections_from_right_to_left = {
  2.0, 1.0,              // A
  1.0, 1.0,              // M
  0.5,                   // r
  {-1.0, 0.0, 0.0, 0.0}, // direction
  5.0                    // distance
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_on_border_outwards = {
  30.0, 10.0,            // A
  20.0, 10.0,            // M
  10.0,                  // r
  {1.0, 0.0, 0.0, 0.0},  // direction
  70.0                   // distance
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_on_border_inwards = {
  30.0, 10.0,            // A
  20.0, 10.0,            // M
  10.0,                  // r
  {-1.0, 0.0, 0.0, 0.0}, // direction
  10.0                   // distance
};

IntersectionProblemParameters_t parameters_for_trajectory_ending_on_border = {
  20.0, 10.0,            // A
  20.0, 10.0,            // M
  10.0,                  // r
  {1.0, 0.0, 0.0, 0.0},  // direction
  10.0                   // distance
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_outside_no_intersection_outwards = {
  15.0, 10.0,            // A
  0.0, 10.0,             // M
  10.0,                  // r
  {1.0, 0.0, 0.0, 0.0},  // direction
  5.0                    // distance
};

IntersectionProblemParameters_t parameters_for_trajectory_starting_outside_no_intersection_towards = {
  20.0, 10.0,            // A
  0.0, 10.0,             // M
  10.0,                  // r
  {-1.0, 0.0, 0.0, 0.0}, // direction
  5.0                    // distance
};

TEST(IntersectionPointsTest, TangentPoint) {
  calculate_intersections(&parameters_for_tangent);
  EXPECT_EQ(intersection_x1(parameters_for_tangent), 1.0);
  EXPECT_EQ(intersection_y1(parameters_for_tangent), 0.5);
  EXPECT_EQ(intersection_x2(parameters_for_tangent), 1.0);
  EXPECT_EQ(intersection_y2(parameters_for_tangent), 0.5);
}
TEST(IntersectionPointsTest, TwoIntersectionPoints) {
  calculate_intersections(&parameters_for_two_intersections);
  EXPECT_EQ(intersection_x1(parameters_for_two_intersections), 0.5);
  EXPECT_EQ(intersection_y1(parameters_for_two_intersections), 1.0);
  EXPECT_EQ(intersection_x2(parameters_for_two_intersections), 1.5);
  EXPECT_EQ(intersection_y2(parameters_for_two_intersections), 1.0);
}
TEST(IntersectionPointsTest, NoIntersectionPoint) {
  calculate_intersections(&parameters_for_no_intersection);
  EXPECT_TRUE(my_is_nan(intersection_x1(parameters_for_no_intersection)));
  EXPECT_TRUE(my_is_nan(intersection_y1(parameters_for_no_intersection)));
  EXPECT_TRUE(my_is_nan(intersection_x2(parameters_for_no_intersection)));
  EXPECT_TRUE(my_is_nan(intersection_y2(parameters_for_no_intersection)));
}
TEST(IntersectionPointsTest, StartingInside) {
  calculate_intersections(&parameters_for_trajectory_starting_inside);
  EXPECT_TRUE(my_is_nan(intersection_x1(parameters_for_trajectory_starting_inside)));
  EXPECT_TRUE(my_is_nan(intersection_y1(parameters_for_trajectory_starting_inside)));
  EXPECT_EQ(intersection_x2(parameters_for_trajectory_starting_inside), 2.0);
  EXPECT_EQ(intersection_y2(parameters_for_trajectory_starting_inside), 1.0);
}

TEST(IntersectionRatioTest, TangentPoint) {
  calculate_intersections(&parameters_for_tangent);
  EXPECT_EQ(intersection_s1(parameters_for_tangent), 0.2);
  EXPECT_EQ(intersection_s2(parameters_for_tangent), 0.2);
}
TEST(IntersectionRatioTest, TwoIntersectionPoints) {
  calculate_intersections(&parameters_for_two_intersections);
  EXPECT_EQ(intersection_s1(parameters_for_two_intersections), 0.1);
  EXPECT_EQ(intersection_s2(parameters_for_two_intersections), 0.3);
}
TEST(IntersectionRatioTest, NoIntersectionPoint) {
  calculate_intersections(&parameters_for_no_intersection);
  EXPECT_TRUE(my_is_nan(intersection_s1(parameters_for_no_intersection)));
  EXPECT_TRUE(my_is_nan(intersection_s2(parameters_for_no_intersection)));
}
TEST(IntersectionRatioTest, StartingInside) {
  calculate_intersections(&parameters_for_trajectory_starting_inside);
  EXPECT_TRUE(intersection_s1(parameters_for_trajectory_starting_inside) < 0);
  EXPECT_EQ(intersection_s2(parameters_for_trajectory_starting_inside), 0.25);
}
TEST(IntersectionRatioTest, TwoIntersectionPointsRtlDirection) {
  calculate_intersections(&parameters_for_two_intersections_from_right_to_left);
  EXPECT_EQ(intersection_s1(parameters_for_two_intersections_from_right_to_left), 0.1);
  EXPECT_EQ(intersection_s2(parameters_for_two_intersections_from_right_to_left), 0.3);
}

TEST(TrajectoryStartsInside, TangentPoint) {
  calculate_intersections(&parameters_for_tangent);
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_tangent));
}
TEST(TrajectoryStartsInside, TwoIntersectionPoints) {
  calculate_intersections(&parameters_for_two_intersections);
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_two_intersections));
}
TEST(TrajectoryStartsInside, NoIntersectionPoint) {
  calculate_intersections(&parameters_for_no_intersection);
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_no_intersection));
}
TEST(TrajectoryStartsInside, StartingInside) {
  calculate_intersections(&parameters_for_trajectory_starting_inside);
  EXPECT_TRUE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_inside));
}
TEST(TrajectoryStartsInside, StartingOnBorderOutwards) {
  calculate_intersections(&parameters_for_trajectory_starting_on_border_outwards);
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_on_border_outwards));
}
TEST(TrajectoryStartsInside, StartingOnBorderInwards) {
  calculate_intersections(&parameters_for_trajectory_starting_on_border_inwards);
  EXPECT_TRUE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_on_border_inwards));
}
TEST(TrajectoryStartsInside, StartingOutsideNoIntersectionOutwards) {
  calculate_intersections(&parameters_for_trajectory_starting_outside_no_intersection_outwards);
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_outside_no_intersection_outwards));
}
TEST(TrajectoryStartsInside, StartingOutsideNoIntersectionTowards) {
  calculate_intersections(&parameters_for_trajectory_starting_outside_no_intersection_towards);
  EXPECT_FALSE(intersecting_trajectory_starts_inside(parameters_for_trajectory_starting_outside_no_intersection_towards));
}

TEST(TrajectoryEndsInside, EndingOnBorder) {
  calculate_intersections(&parameters_for_trajectory_ending_on_border);
  EXPECT_TRUE(intersecting_trajectory_ends_inside(parameters_for_trajectory_ending_on_border));
}

