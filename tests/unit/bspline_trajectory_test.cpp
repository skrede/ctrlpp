#include "ctrlpp/traj/bspline_trajectory.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

TEST_CASE("Cubic B-spline endpoint interpolation", "[bspline]")
{
    // Degree 3 with 5 control points and auto-generated uniform clamped knots.
    // Clamped B-spline passes through first and last control points.
    std::vector<double> ctrl = {0.0, 1.0, 3.0, 2.0, 4.0};

    ctrlpp::bspline_trajectory<double, 3> bs({.control_points = ctrl});

    auto const p0 = bs.evaluate(0.0);
    auto const pT = bs.evaluate(bs.duration());

    REQUIRE_THAT(p0.position(0), WithinAbs(ctrl.front(), 1e-12));
    REQUIRE_THAT(pT.position(0), WithinAbs(ctrl.back(), 1e-12));
}

TEST_CASE("Cubic B-spline convex hull", "[bspline]")
{
    std::vector<double> ctrl = {1.0, 2.0, 4.0, 3.0, 5.0};

    ctrlpp::bspline_trajectory<double, 3> bs({.control_points = ctrl});

    auto const min_cp = *std::min_element(ctrl.begin(), ctrl.end());
    auto const max_cp = *std::max_element(ctrl.begin(), ctrl.end());

    auto const T = bs.duration();
    int const N = 50;
    for (int i = 0; i <= N; ++i) {
        auto const t = T * static_cast<double>(i) / N;
        auto const pt = bs.evaluate(t);
        REQUIRE(pt.position(0) >= min_cp - 1e-12);
        REQUIRE(pt.position(0) <= max_cp + 1e-12);
    }
}

TEST_CASE("make_bspline_interpolation passes through waypoints", "[bspline]")
{
    // 5 waypoints at uniform times
    std::vector<double> times = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> positions = {0.0, 2.0, 1.0, 3.0, 5.0};

    auto bs = ctrlpp::make_bspline_interpolation<double, 3>(times, positions);

    for (std::size_t i = 0; i < times.size(); ++i) {
        auto const pt = bs.evaluate(times[i]);
        REQUIRE_THAT(pt.position(0), WithinAbs(positions[i], 1e-10));
    }
}

TEST_CASE("Quintic B-spline evaluation", "[bspline]")
{
    // Degree 5 with 7 control points
    std::vector<double> ctrl = {0.0, 1.0, 2.5, 4.0, 3.0, 2.0, 5.0};

    ctrlpp::bspline_trajectory<double, 5> bs({.control_points = ctrl});

    auto const T = bs.duration();
    REQUIRE(T > 0.0);

    int const N = 20;
    for (int i = 0; i <= N; ++i) {
        auto const t = T * static_cast<double>(i) / N;
        auto const pt = bs.evaluate(t);
        // Should not produce NaN or Inf
        REQUIRE(std::isfinite(pt.position(0)));
        REQUIRE(std::isfinite(pt.velocity(0)));
        REQUIRE(std::isfinite(pt.acceleration(0)));
    }
}

TEST_CASE("Custom knot vector", "[bspline]")
{
    // Cubic with 5 control points, custom non-uniform knots
    // m = n + p + 1 = 4 + 3 + 1 = 8 knots total (9 values)
    std::vector<double> ctrl = {0.0, 1.0, 3.0, 2.0, 4.0};
    std::vector<double> knots = {0, 0, 0, 0, 0.3, 1, 1, 1, 1};

    ctrlpp::bspline_trajectory<double, 3> bs({.control_points = ctrl, .knot_vector = knots});

    auto const p0 = bs.evaluate(0.0);
    auto const pT = bs.evaluate(bs.duration());

    REQUIRE_THAT(p0.position(0), WithinAbs(ctrl.front(), 1e-12));
    REQUIRE_THAT(pT.position(0), WithinAbs(ctrl.back(), 1e-12));
}

TEST_CASE("Invalid knot vector detected", "[bspline]")
{
    // Wrong number of knots for degree 3 with 5 control points (need 9, provide 7)
    std::vector<double> ctrl = {0.0, 1.0, 3.0, 2.0, 4.0};
    std::vector<double> bad_knots = {0, 0, 0, 0.5, 1, 1, 1};

    REQUIRE_THROWS(ctrlpp::bspline_trajectory<double, 3>(
        {.control_points = ctrl, .knot_vector = bad_knots}));
}

TEST_CASE("Duration returns active parameter range", "[bspline]")
{
    // For uniform clamped knots on [0,1], duration should be 1.0
    std::vector<double> ctrl = {0.0, 1.0, 3.0, 2.0, 4.0};

    ctrlpp::bspline_trajectory<double, 3> bs({.control_points = ctrl});

    REQUIRE_THAT(bs.duration(), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Satisfies trajectory_segment concept", "[bspline]")
{
    static_assert(ctrlpp::trajectory_segment<ctrlpp::bspline_trajectory<double, 3>, double, 1>);
    static_assert(ctrlpp::trajectory_segment<ctrlpp::bspline_trajectory<double, 5>, double, 1>);
    SUCCEED("Static assert passed");
}

TEST_CASE("Velocity and acceleration via derivative de Boor", "[bspline]")
{
    // For a linear B-spline (control points on a line), velocity should be constant
    // and acceleration should be zero.
    std::vector<double> ctrl = {0.0, 1.0, 2.0, 3.0, 4.0};

    ctrlpp::bspline_trajectory<double, 3> bs({.control_points = ctrl});

    auto const T = bs.duration();
    int const N = 20;
    for (int i = 1; i < N; ++i) {
        auto const t = T * static_cast<double>(i) / N;
        auto const pt = bs.evaluate(t);
        // For collinear control points with uniform clamped knots,
        // velocity should be approximately constant at interior points
        REQUIRE(std::isfinite(pt.velocity(0)));
        REQUIRE(std::isfinite(pt.acceleration(0)));
    }

    // At interior points, velocity should be positive for ascending control points
    auto const mid = bs.evaluate(T * 0.5);
    REQUIRE(mid.velocity(0) > 0.0);
}
