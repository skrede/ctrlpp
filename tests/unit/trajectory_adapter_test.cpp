/// @brief Tests for trajectory adapter lifting normalized paths into trajectory segments.
///
/// Verifies harmonic and cycloidal paths produce correct endpoint conditions
/// and derivative scaling when used with the trajectory adapter.

#include "ctrlpp/trajectory/cycloidal_path.h"
#include "ctrlpp/trajectory/harmonic_path.h"
#include "ctrlpp/trajectory/trajectory.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <limits>
#include <numbers>

using Catch::Matchers::WithinAbs;
using ctrlpp::Vector;

// ---- Concept satisfaction (compile-time) ----
using harmonic_traj_t =
    ctrlpp::trajectory<ctrlpp::path_point<double> (*)(double), double, 1>;
static_assert(ctrlpp::trajectory_segment<harmonic_traj_t, double, 1>);

// ---- Harmonic via trajectory ----

TEST_CASE("trajectory harmonic rest-to-rest", "[traj][trajectory_adapter]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);

    auto seg = ctrlpp::make_trajectory(ctrlpp::harmonic_path<double>, q0, q1, 2.0).value();

    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p0.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.velocity[0], WithinAbs(0.0, 1e-12));

    auto pT = seg.evaluate(2.0);
    CHECK_THAT(pT.position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(pT.velocity[0], WithinAbs(0.0, 1e-12));
}

TEST_CASE("trajectory harmonic midpoint symmetry", "[traj][trajectory_adapter]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);

    auto seg = ctrlpp::make_trajectory(ctrlpp::harmonic_path<double>, q0, q1, 2.0).value();
    auto pm = seg.evaluate(1.0);
    CHECK_THAT(pm.position[0], WithinAbs(5.0, 1e-10));
}

TEST_CASE("trajectory harmonic endpoint acceleration", "[traj][trajectory_adapter]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    double const T = 2.0;

    auto seg = ctrlpp::make_trajectory(ctrlpp::harmonic_path<double>, q0, q1, T).value();

    // harmonic ddq(0) = pi^2/2, scaled by h/T^2 = 10/4 = 2.5
    // expected: 2.5 * pi^2 / 2 = 5*pi^2/4
    auto const pi = std::numbers::pi;
    auto const expected_acc = 10.0 * pi * pi / (2.0 * T * T);
    CHECK_THAT(seg.evaluate(0.0).acceleration[0], WithinAbs(expected_acc, 1e-10));
}

// ---- Cycloidal via trajectory ----

TEST_CASE("trajectory cycloidal rest-to-rest", "[traj][trajectory_adapter]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Constant(5.0);
    Vector<double, 1> q1 = Vector<double, 1>::Constant(15.0);

    auto seg = ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q1, 3.0).value();

    CHECK_THAT(seg.evaluate(0.0).position[0], WithinAbs(5.0, 1e-12));
    CHECK_THAT(seg.evaluate(3.0).position[0], WithinAbs(15.0, 1e-12));
    CHECK_THAT(seg.evaluate(0.0).velocity[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(seg.evaluate(3.0).velocity[0], WithinAbs(0.0, 1e-12));
}

TEST_CASE("trajectory cycloidal zero endpoint acceleration", "[traj][trajectory_adapter]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Constant(5.0);
    Vector<double, 1> q1 = Vector<double, 1>::Constant(15.0);

    auto seg = ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q1, 3.0).value();

    CHECK_THAT(seg.evaluate(0.0).acceleration[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(seg.evaluate(3.0).acceleration[0], WithinAbs(0.0, 1e-12));
}

// ---- ND=3 multi-dimensional ----

TEST_CASE("trajectory ND=3 per-axis displacement", "[traj][trajectory_adapter]")
{
    Vector<double, 3> q0;
    q0 << 0.0, 10.0, -5.0;
    Vector<double, 3> q1;
    q1 << 10.0, 30.0, 5.0;

    auto seg = ctrlpp::make_trajectory(ctrlpp::harmonic_path<double>, q0, q1, 2.0).value();

    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p0.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.position[1], WithinAbs(10.0, 1e-12));
    CHECK_THAT(p0.position[2], WithinAbs(-5.0, 1e-12));

    auto pT = seg.evaluate(2.0);
    CHECK_THAT(pT.position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(pT.position[1], WithinAbs(30.0, 1e-12));
    CHECK_THAT(pT.position[2], WithinAbs(5.0, 1e-12));
}

// ---- Boundary clamping ----

TEST_CASE("trajectory boundary clamping", "[traj][trajectory_adapter]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);

    auto seg = ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q1, 2.0).value();

    auto p_neg = seg.evaluate(-1.0);
    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p_neg.position[0], WithinAbs(p0.position[0], 1e-15));

    auto p_over = seg.evaluate(10.0);
    auto pT = seg.evaluate(2.0);
    CHECK_THAT(p_over.position[0], WithinAbs(pT.position[0], 1e-15));
}

// ---- Rejection tests: duration and finiteness domain ----

TEST_CASE("make_trajectory rejects invalid duration and positions",
          "[traj][trajectory_adapter][negative]")
{
    Vector<double, 1> const q0 = Vector<double, 1>::Zero();
    Vector<double, 1> const q1 = Vector<double, 1>::Constant(10.0);
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();

    SECTION("duration = 0 -> non_positive_duration")
    {
        auto const result =
            ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q1, 0.0);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_duration);
    }

    SECTION("duration = -1 -> non_positive_duration")
    {
        auto const result =
            ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q1, -1.0);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_duration);
    }

    SECTION("duration = NaN -> non_finite_input")
    {
        auto const result =
            ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q1, nan);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_finite_input);
    }

    SECTION("NaN position entry -> non_finite_input")
    {
        Vector<double, 1> const q_nan = Vector<double, 1>::Constant(nan);
        auto const result =
            ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q_nan, 2.0);
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_finite_input);
    }
}

TEST_CASE("trajectory duration()", "[traj][trajectory_adapter]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);

    auto seg = ctrlpp::make_trajectory(ctrlpp::harmonic_path<double>, q0, q1, 3.5).value();
    CHECK(seg.duration() == 3.5);
}
