/// @brief Tests for cubic, quintic, and septic polynomial trajectory segments.
///
/// Verifies boundary conditions, midpoint symmetry, clamping, and concept satisfaction.

#include "ctrlpp/traj/cubic_trajectory.h"
#include "ctrlpp/traj/quintic_trajectory.h"
#include "ctrlpp/traj/septic_trajectory.h"
#include "ctrlpp/traj/trajectory_segment.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using ctrlpp::Vector;

// ---- Concept satisfaction (compile-time) ----
static_assert(ctrlpp::trajectory_segment<ctrlpp::cubic_trajectory<double, 1>, double, 1>);
static_assert(ctrlpp::trajectory_segment<ctrlpp::quintic_trajectory<double, 1>, double, 1>);
static_assert(ctrlpp::trajectory_segment<ctrlpp::septic_trajectory<double, 1>, double, 1>);
static_assert(ctrlpp::trajectory_segment<ctrlpp::cubic_trajectory<double, 3>, double, 3>);

// ---- Cubic trajectory tests ----

TEST_CASE("cubic_trajectory rest-to-rest boundary conditions", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> v0 = Vector<double, 1>::Zero();
    Vector<double, 1> v1 = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_cubic_trajectory(q0, q1, v0, v1, 2.0);

    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p0.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.velocity[0], WithinAbs(0.0, 1e-12));

    auto pT = seg.evaluate(2.0);
    CHECK_THAT(pT.position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(pT.velocity[0], WithinAbs(0.0, 1e-12));
}

TEST_CASE("cubic_trajectory midpoint symmetry", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> v0 = Vector<double, 1>::Zero();
    Vector<double, 1> v1 = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_cubic_trajectory(q0, q1, v0, v1, 2.0);
    auto pm = seg.evaluate(1.0);
    CHECK_THAT(pm.position[0], WithinAbs(5.0, 1e-10));
}

TEST_CASE("cubic_trajectory boundary clamping", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> v0 = Vector<double, 1>::Zero();
    Vector<double, 1> v1 = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_cubic_trajectory(q0, q1, v0, v1, 2.0);

    auto p_neg = seg.evaluate(-1.0);
    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p_neg.position[0], WithinAbs(p0.position[0], 1e-15));

    auto p_over = seg.evaluate(10.0);
    auto pT = seg.evaluate(2.0);
    CHECK_THAT(p_over.position[0], WithinAbs(pT.position[0], 1e-15));
}

TEST_CASE("cubic_trajectory non-zero velocity BCs", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> v0 = Vector<double, 1>::Constant(2.0);
    Vector<double, 1> v1 = Vector<double, 1>::Constant(-3.0);

    auto seg = ctrlpp::make_cubic_trajectory(q0, q1, v0, v1, 4.0);

    CHECK_THAT(seg.evaluate(0.0).velocity[0], WithinAbs(2.0, 1e-12));
    CHECK_THAT(seg.evaluate(4.0).velocity[0], WithinAbs(-3.0, 1e-12));
    CHECK_THAT(seg.evaluate(0.0).position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(seg.evaluate(4.0).position[0], WithinAbs(10.0, 1e-12));
}

TEST_CASE("cubic_trajectory duration()", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> v0 = Vector<double, 1>::Zero();
    Vector<double, 1> v1 = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_cubic_trajectory(q0, q1, v0, v1, 2.0);
    CHECK(seg.duration() == 2.0);
}

// ---- Quintic trajectory tests ----

TEST_CASE("quintic_trajectory rest-to-rest boundary conditions", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> v0 = Vector<double, 1>::Zero();
    Vector<double, 1> v1 = Vector<double, 1>::Zero();
    Vector<double, 1> a0 = Vector<double, 1>::Zero();
    Vector<double, 1> a1 = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_quintic_trajectory(q0, q1, v0, v1, a0, a1, 2.0);

    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p0.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.velocity[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.acceleration[0], WithinAbs(0.0, 1e-12));

    auto pT = seg.evaluate(2.0);
    CHECK_THAT(pT.position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(pT.velocity[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(pT.acceleration[0], WithinAbs(0.0, 1e-12));
}

TEST_CASE("quintic_trajectory midpoint symmetry", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> zero = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_quintic_trajectory(q0, q1, zero, zero, zero, zero, 2.0);
    auto pm = seg.evaluate(1.0);
    CHECK_THAT(pm.position[0], WithinAbs(5.0, 1e-10));
}

TEST_CASE("quintic_trajectory boundary clamping", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> zero = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_quintic_trajectory(q0, q1, zero, zero, zero, zero, 2.0);

    auto p_neg = seg.evaluate(-1.0);
    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p_neg.position[0], WithinAbs(p0.position[0], 1e-15));

    auto p_over = seg.evaluate(10.0);
    auto pT = seg.evaluate(2.0);
    CHECK_THAT(p_over.position[0], WithinAbs(pT.position[0], 1e-15));
}

TEST_CASE("quintic_trajectory non-zero BCs", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> v0 = Vector<double, 1>::Constant(1.0);
    Vector<double, 1> v1 = Vector<double, 1>::Constant(-1.0);
    Vector<double, 1> a0 = Vector<double, 1>::Constant(0.5);
    Vector<double, 1> a1 = Vector<double, 1>::Constant(-0.5);

    auto seg = ctrlpp::make_quintic_trajectory(q0, q1, v0, v1, a0, a1, 3.0);

    CHECK_THAT(seg.evaluate(0.0).position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(seg.evaluate(3.0).position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(seg.evaluate(0.0).velocity[0], WithinAbs(1.0, 1e-12));
    CHECK_THAT(seg.evaluate(3.0).velocity[0], WithinAbs(-1.0, 1e-12));
    CHECK_THAT(seg.evaluate(0.0).acceleration[0], WithinAbs(0.5, 1e-12));
    CHECK_THAT(seg.evaluate(3.0).acceleration[0], WithinAbs(-0.5, 1e-12));
}

// ---- Septic trajectory tests ----

TEST_CASE("septic_trajectory rest-to-rest boundary conditions", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> zero = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_septic_trajectory(q0, q1, zero, zero, zero, zero, zero, zero, 2.0);

    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p0.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.velocity[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.acceleration[0], WithinAbs(0.0, 1e-12));

    auto pT = seg.evaluate(2.0);
    CHECK_THAT(pT.position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(pT.velocity[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(pT.acceleration[0], WithinAbs(0.0, 1e-12));
}

TEST_CASE("septic_trajectory midpoint symmetry", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> zero = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_septic_trajectory(q0, q1, zero, zero, zero, zero, zero, zero, 2.0);
    auto pm = seg.evaluate(1.0);
    CHECK_THAT(pm.position[0], WithinAbs(5.0, 1e-10));
}

TEST_CASE("septic_trajectory boundary clamping", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> zero = Vector<double, 1>::Zero();

    auto seg = ctrlpp::make_septic_trajectory(q0, q1, zero, zero, zero, zero, zero, zero, 2.0);

    auto p_neg = seg.evaluate(-1.0);
    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p_neg.position[0], WithinAbs(p0.position[0], 1e-15));

    auto p_over = seg.evaluate(10.0);
    auto pT = seg.evaluate(2.0);
    CHECK_THAT(p_over.position[0], WithinAbs(pT.position[0], 1e-15));
}

TEST_CASE("septic_trajectory non-zero jerk BCs", "[traj][polynomial]")
{
    Vector<double, 1> q0 = Vector<double, 1>::Zero();
    Vector<double, 1> q1 = Vector<double, 1>::Constant(10.0);
    Vector<double, 1> zero = Vector<double, 1>::Zero();
    Vector<double, 1> j0 = Vector<double, 1>::Constant(1.0);
    Vector<double, 1> j1 = Vector<double, 1>::Constant(-1.0);

    auto seg = ctrlpp::make_septic_trajectory(q0, q1, zero, zero, zero, zero, j0, j1, 2.0);

    CHECK_THAT(seg.evaluate(0.0).position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(seg.evaluate(2.0).position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(seg.evaluate(0.0).velocity[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(seg.evaluate(2.0).velocity[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(seg.evaluate(0.0).acceleration[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(seg.evaluate(2.0).acceleration[0], WithinAbs(0.0, 1e-12));
}

// ---- ND=3 multi-dimensional test ----

TEST_CASE("cubic_trajectory ND=3 independent axes", "[traj][polynomial]")
{
    Vector<double, 3> q0;
    q0 << 0.0, 0.0, 0.0;
    Vector<double, 3> q1;
    q1 << 10.0, 20.0, 30.0;
    Vector<double, 3> v0;
    v0 << 1.0, 0.0, -1.0;
    Vector<double, 3> v1;
    v1 << -1.0, 0.0, 1.0;

    auto seg = ctrlpp::make_cubic_trajectory(q0, q1, v0, v1, 2.0);

    auto p0 = seg.evaluate(0.0);
    CHECK_THAT(p0.position[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.position[1], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.position[2], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.velocity[0], WithinAbs(1.0, 1e-12));
    CHECK_THAT(p0.velocity[1], WithinAbs(0.0, 1e-12));
    CHECK_THAT(p0.velocity[2], WithinAbs(-1.0, 1e-12));

    auto pT = seg.evaluate(2.0);
    CHECK_THAT(pT.position[0], WithinAbs(10.0, 1e-12));
    CHECK_THAT(pT.position[1], WithinAbs(20.0, 1e-12));
    CHECK_THAT(pT.position[2], WithinAbs(30.0, 1e-12));
    CHECK_THAT(pT.velocity[0], WithinAbs(-1.0, 1e-12));
    CHECK_THAT(pT.velocity[1], WithinAbs(0.0, 1e-12));
    CHECK_THAT(pT.velocity[2], WithinAbs(1.0, 1e-12));
}
