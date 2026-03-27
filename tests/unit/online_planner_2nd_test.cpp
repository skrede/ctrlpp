#include "ctrlpp/traj/online_planner_2nd.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinAbs;

// -- Test 1: Step response settles to target -----------------------------------
TEST_CASE("OnlinePlanner2nd: step response settles", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<double> planner({.v_max = 5.0, .a_max = 10.0});

    planner.update(10.0);

    // Sample far enough in the future that the planner should have settled
    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(10.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));
    REQUIRE(planner.is_settled());
}

// -- Test 2: Velocity never exceeds v_max --------------------------------------
TEST_CASE("OnlinePlanner2nd: velocity constraint", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<double> planner({.v_max = 5.0, .a_max = 10.0});

    planner.update(10.0);

    double constexpr tol = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = 10.0 * static_cast<double>(i) / 1000.0;
        auto const pt = planner.sample(t);
        REQUIRE(std::abs(pt.velocity[0]) <= 5.0 + tol);
    }
}

// -- Test 3: Acceleration never exceeds a_max ----------------------------------
TEST_CASE("OnlinePlanner2nd: acceleration constraint", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<double> planner({.v_max = 5.0, .a_max = 10.0});

    planner.update(10.0);

    double constexpr tol = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = 10.0 * static_cast<double>(i) / 1000.0;
        auto const pt = planner.sample(t);
        REQUIRE(std::abs(pt.acceleration[0]) <= 10.0 + tol);
    }
}

// -- Test 4: Mid-motion target change ------------------------------------------
TEST_CASE("OnlinePlanner2nd: mid-motion target change", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<double> planner({.v_max = 5.0, .a_max = 10.0});

    planner.update(10.0);

    // Sample partway to build up state
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

// -- Test 5: Negative displacement ---------------------------------------------
TEST_CASE("OnlinePlanner2nd: negative displacement", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<double> planner({.v_max = 5.0, .a_max = 10.0});

    planner.update(-5.0);

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(-5.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));

    // Velocity should go negative during motion
    auto const mid = planner.sample(0.5);
    REQUIRE(mid.velocity[0] < 0.0);
}

// -- Test 6: Zero displacement -- immediately settled --------------------------
TEST_CASE("OnlinePlanner2nd: zero displacement", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<double> planner({.v_max = 5.0, .a_max = 10.0});

    planner.update(0.0);

    REQUIRE(planner.is_settled());
    auto const pt = planner.sample(0.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-12));
}

// -- Test 7: is_settled returns correct state ----------------------------------
TEST_CASE("OnlinePlanner2nd: is_settled transitions", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<double> planner({.v_max = 5.0, .a_max = 10.0});

    // Initially at rest at origin -- settled
    REQUIRE(planner.is_settled());

    planner.update(10.0);

    // After update with different target -- not settled (sampling before arrival)
    planner.sample(0.1);
    REQUIRE_FALSE(planner.is_settled());

    // After sufficient time -- settled
    planner.sample(100.0);
    REQUIRE(planner.is_settled());
}

// -- Test 8: Reset functionality -----------------------------------------------
TEST_CASE("OnlinePlanner2nd: reset", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<double> planner({.v_max = 5.0, .a_max = 10.0});

    planner.update(10.0);
    planner.sample(100.0);

    planner.reset(3.0);
    REQUIRE(planner.is_settled());

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(3.0, 1e-12));
}

// -- Test 9: Profile shape is trapezoidal-like ---------------------------------
TEST_CASE("OnlinePlanner2nd: trapezoidal profile shape", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<double> planner({.v_max = 5.0, .a_max = 10.0});

    planner.update(10.0);

    // During acceleration phase: velocity increases, acceleration is positive
    auto const early = planner.sample(0.1);
    REQUIRE(early.velocity[0] > 0.0);
    REQUIRE(early.acceleration[0] > 0.0);

    // During cruise (if long enough): acceleration should be near zero
    // For h=10, v_max=5, a_max=10: T_a=0.5, T_v=1.5, T_d=0.5
    auto const cruise = planner.sample(1.0);
    REQUIRE_THAT(cruise.velocity[0], WithinAbs(5.0, 0.1));
    REQUIRE_THAT(cruise.acceleration[0], WithinAbs(0.0, 0.1));
}

// -- Test 10: Float type works -------------------------------------------------
TEST_CASE("OnlinePlanner2nd: float type", "[traj][online_planner_2nd]")
{
    ctrlpp::online_planner_2nd<float> planner({.v_max = 5.0f, .a_max = 10.0f});

    planner.update(10.0f);
    auto const pt = planner.sample(100.0f);
    REQUIRE_THAT(pt.position[0], WithinAbs(10.0, 1e-3));
}
