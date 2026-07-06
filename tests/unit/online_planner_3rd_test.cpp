#include "ctrlpp/trajectory/online_planner_3rd.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>

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

// -- Test 11: Same-direction retarget carries velocity through (no full-stop dip)
//
// Retargeting farther in the same direction while cruising must NOT brake the
// motion to a full stop and re-accelerate. A time-optimal planner keeps the
// current velocity and simply extends the profile, so once the planner is
// cruising the velocity never rises again after the retarget.
TEST_CASE("OnlinePlanner3rd: same-direction retarget keeps velocity",
          "[traj][online_planner_3rd]")
{
    double constexpr v_max = 5.0;
    ctrlpp::online_planner_3rd<double> planner({.v_max = v_max, .a_max = 10.0, .j_max = 50.0});

    // Command a far target and cruise up to v_max.
    planner.update(100.0);

    double constexpr dt = 0.01;
    double constexpr t_retarget = 2.0;
    double t = 0.0;
    for (; t < t_retarget - 0.5 * dt; t += dt) {
        planner.sample(t);
    }

    // Confirm the planner is cruising at v_max with (near) zero acceleration
    // before the retarget, so the "no velocity rise afterward" check is clean.
    auto const at_retarget = planner.sample(t_retarget);
    REQUIRE_THAT(at_retarget.velocity[0], WithinAbs(v_max, 1e-6));
    REQUIRE_THAT(at_retarget.acceleration[0], WithinAbs(0.0, 1e-6));

    // Retarget farther in the same direction (non-overshoot, same sign).
    planner.update(200.0);

    // Immediately after the retarget the planner should keep cruising, not brake.
    auto const shortly_after = planner.sample(t_retarget + 1.0);
    REQUIRE(shortly_after.velocity[0] > 0.9 * v_max);

    // Sweep to settling: after retarget-at-cruise the velocity is monotonically
    // non-increasing (flat cruise, then a single deceleration to rest). The old
    // brake-to-rest-then-replan strategy would dip to zero and rise again.
    // The Lipschitz bound |dv| <= a_max*dt also guards against the end-of-profile
    // snap masking a profile that fails to bring the velocity to zero at T.
    double constexpr a_max = 10.0;
    double prev_v = at_retarget.velocity[0];
    double constexpr tol = 1e-6;
    for (t = t_retarget + dt; t < 80.0; t += dt) {
        auto const pt = planner.sample(t);
        REQUIRE(pt.velocity[0] <= prev_v + tol);
        REQUIRE(std::abs(pt.velocity[0]) <= v_max + tol);
        REQUIRE(std::abs(pt.velocity[0] - prev_v) <= a_max * dt + 1e-3);
        prev_v = pt.velocity[0];
    }

    auto const settled = planner.sample(200.0);
    REQUIRE_THAT(settled.position[0], WithinAbs(200.0, 1e-6));
    REQUIRE_THAT(settled.velocity[0], WithinAbs(0.0, 1e-6));
}

// -- Test 12: try_create rejects out-of-domain limits ----------------------------
//
// All three limits divide in the planner math (cruise duration h/v_max,
// jerk-phase durations a_max/j_max and |a|/j_max), so the domain of each is
// finite and strictly positive; everything else is rejected with the
// limit-specific enumerator.
TEST_CASE("OnlinePlanner3rd: try_create rejects invalid limits",
          "[traj][online_planner_3rd][negative]")
{
    using planner_t = ctrlpp::online_planner_3rd<double>;
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    SECTION("invalid v_max -> non_positive_velocity_limit")
    {
        for (double const v_max : {0.0, -1.0, nan, inf}) {
            auto const result =
                planner_t::try_create({.v_max = v_max, .a_max = 10.0, .j_max = 50.0});
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error()
                    == ctrlpp::trajectory_error::non_positive_velocity_limit);
        }
    }

    SECTION("invalid a_max -> non_positive_acceleration_limit")
    {
        for (double const a_max : {0.0, -1.0, nan, inf}) {
            auto const result =
                planner_t::try_create({.v_max = 5.0, .a_max = a_max, .j_max = 50.0});
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error()
                    == ctrlpp::trajectory_error::non_positive_acceleration_limit);
        }
    }

    SECTION("invalid j_max -> non_positive_jerk_limit")
    {
        for (double const j_max : {0.0, -1.0, nan, inf}) {
            auto const result =
                planner_t::try_create({.v_max = 5.0, .a_max = 10.0, .j_max = j_max});
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_jerk_limit);
        }
    }
}

// -- Test 13: factory-built and ctor-built planners agree ------------------------
TEST_CASE("OnlinePlanner3rd: try_create and constructor produce identical profiles",
          "[traj][online_planner_3rd]")
{
    auto factory_built = ctrlpp::online_planner_3rd<double>::try_create(
        {.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});
    REQUIRE(factory_built.has_value());
    ctrlpp::online_planner_3rd<double> ctor_built(
        {.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    factory_built->update(10.0);
    ctor_built.update(10.0);

    double constexpr dt = 0.01;
    double t = 0.0;
    for (int i = 0; i < 100; ++i) {
        t += dt;
        auto const a = factory_built->sample(t);
        auto const b = ctor_built.sample(t);
        REQUIRE(a.position[0] == b.position[0]);
        REQUIRE(a.velocity[0] == b.velocity[0]);
        REQUIRE(a.acceleration[0] == b.acceleration[0]);
    }

    // Retarget mid-motion and keep comparing.
    factory_built->update(-5.0);
    ctor_built.update(-5.0);
    for (int i = 0; i < 100; ++i) {
        t += dt;
        auto const a = factory_built->sample(t);
        auto const b = ctor_built.sample(t);
        REQUIRE(a.position[0] == b.position[0]);
        REQUIRE(a.velocity[0] == b.velocity[0]);
        REQUIRE(a.acceleration[0] == b.acceleration[0]);
    }
}
