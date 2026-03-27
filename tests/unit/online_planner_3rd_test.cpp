#include "ctrlpp/traj/online_planner_3rd.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;

// -- Test 1: Step response settles to target -----------------------------------
TEST_CASE("OnlinePlanner3rd: step response settles", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<double> planner({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    planner.update(10.0);

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(10.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(pt.acceleration[0], WithinAbs(0.0, 1e-6));
    REQUIRE(planner.is_settled());
}

// -- Test 2: Velocity never exceeds v_max --------------------------------------
TEST_CASE("OnlinePlanner3rd: velocity constraint", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<double> planner({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    planner.update(10.0);

    double constexpr tol = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = 10.0 * static_cast<double>(i) / 1000.0;
        auto const pt = planner.sample(t);
        REQUIRE(std::abs(pt.velocity[0]) <= 5.0 + tol);
    }
}

// -- Test 3: Acceleration never exceeds a_max ----------------------------------
TEST_CASE("OnlinePlanner3rd: acceleration constraint", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<double> planner({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    planner.update(10.0);

    double constexpr tol = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = 10.0 * static_cast<double>(i) / 1000.0;
        auto const pt = planner.sample(t);
        REQUIRE(std::abs(pt.acceleration[0]) <= 10.0 + tol);
    }
}

// -- Test 4: Jerk never exceeds j_max (numerical differentiation) ---------------
TEST_CASE("OnlinePlanner3rd: jerk constraint", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<double> planner({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    planner.update(10.0);

    double constexpr tol = 1e-3; // numerical differentiation tolerance
    double constexpr small_dt = 1e-5;

    for (int i = 0; i < 1000; ++i) {
        double const t = 10.0 * static_cast<double>(i) / 1000.0;
        auto const pt0 = planner.sample(t);
        auto const pt1 = planner.sample(t + small_dt);
        double const jerk = (pt1.acceleration[0] - pt0.acceleration[0]) / small_dt;
        REQUIRE(std::abs(jerk) <= 50.0 + tol);
    }
}

// -- Test 5: Mid-motion target change ------------------------------------------
TEST_CASE("OnlinePlanner3rd: mid-motion target change", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<double> planner({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    planner.update(10.0);

    // Sample partway
    auto const mid = planner.sample(0.5);
    REQUIRE(mid.position[0] > 0.0);

    // Change target mid-motion
    planner.update(5.0);

    // Eventually settles to new target
    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(5.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));

    // Check constraints after target change
    double constexpr tol = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = 0.5 + 10.0 * static_cast<double>(i) / 1000.0;
        auto const p = planner.sample(t);
        REQUIRE(std::abs(p.velocity[0]) <= 5.0 + tol);
        REQUIRE(std::abs(p.acceleration[0]) <= 10.0 + tol);
    }
}

// -- Test 6: Negative displacement ---------------------------------------------
TEST_CASE("OnlinePlanner3rd: negative displacement", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<double> planner({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    planner.update(-5.0);

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(-5.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));
}

// -- Test 7: is_settled returns correct state ----------------------------------
TEST_CASE("OnlinePlanner3rd: is_settled transitions", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<double> planner({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.is_settled());

    planner.update(10.0);
    planner.sample(0.1);
    REQUIRE_FALSE(planner.is_settled());

    planner.sample(100.0);
    REQUIRE(planner.is_settled());
}

// -- Test 8: Reset functionality -----------------------------------------------
TEST_CASE("OnlinePlanner3rd: reset", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<double> planner({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    planner.update(10.0);
    planner.sample(100.0);

    planner.reset(3.0);
    REQUIRE(planner.is_settled());

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(3.0, 1e-12));
}

// -- Test 9: Profile has smooth acceleration (double-S shape) ------------------
TEST_CASE("OnlinePlanner3rd: smooth acceleration profile", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<double> planner({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    planner.update(10.0);

    // Early: acceleration should be ramping up (jerk > 0)
    auto const early = planner.sample(0.05);
    REQUIRE(early.velocity[0] > 0.0);
    REQUIRE(early.acceleration[0] > 0.0);

    // Acceleration should be smooth (not step-like as in 2nd-order)
    auto const a1 = planner.sample(0.01).acceleration[0];
    auto const a2 = planner.sample(0.02).acceleration[0];
    auto const a3 = planner.sample(0.03).acceleration[0];
    // Acceleration should be monotonically increasing during jerk phase
    REQUIRE(a2 > a1);
    REQUIRE(a3 > a2);
}

// -- Test 10: Float type works -------------------------------------------------
TEST_CASE("OnlinePlanner3rd: float type", "[traj][online_planner_3rd]")
{
    ctrlpp::online_planner_3rd<float> planner({.v_max = 5.0f, .a_max = 10.0f, .j_max = 50.0f});

    planner.update(10.0f);
    auto const pt = planner.sample(100.0f);
    REQUIRE_THAT(pt.position[0], WithinAbs(10.0, 1e-3));
}
