#include "ctrlpp/trajectory/online_planner_2nd.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <limits>
#include <utility>

using Catch::Matchers::WithinAbs;

namespace
{

// Construction is fallible, so every valid-input site goes through the factory
// and asserts success here; a config that becomes unrealizable fails the test
// instead of quietly skipping it. Rejection cases never use this helper.
template <typename Scalar>
auto make_planner(typename ctrlpp::online_planner_2nd<Scalar>::config const& cfg)
    -> ctrlpp::online_planner_2nd<Scalar>
{
    auto created = ctrlpp::online_planner_2nd<Scalar>::create(cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

}

// -- Test 1: Step response settles to target -----------------------------------
TEST_CASE("OnlinePlanner2nd: step response settles", "[traj][online_planner_2nd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

    planner.update(0.0);

    REQUIRE(planner.is_settled());
    auto const pt = planner.sample(0.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-12));
}

// -- Test 7: is_settled returns correct state ----------------------------------
TEST_CASE("OnlinePlanner2nd: is_settled transitions", "[traj][online_planner_2nd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

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
    auto planner = make_planner<float>({.v_max = 5.0f, .a_max = 10.0f});

    planner.update(10.0f);
    auto const pt = planner.sample(100.0f);
    REQUIRE_THAT(pt.position[0], WithinAbs(10.0, 1e-3));
}

// -- Test 11: Same-direction retarget carries velocity through (no full-stop dip)
//
// Retargeting farther in the same direction while cruising must NOT brake the
// motion to a full stop and re-accelerate. A time-optimal planner keeps the
// current velocity and extends the profile, so once cruising the velocity never
// rises again after the retarget.
TEST_CASE("OnlinePlanner2nd: same-direction retarget keeps velocity",
          "[traj][online_planner_2nd]")
{
    double constexpr v_max = 5.0;
    double constexpr a_max = 10.0;
    auto planner = make_planner<double>({.v_max = v_max, .a_max = a_max});

    // Command a far target and cruise up to v_max.
    planner.update(100.0);

    double constexpr dt = 0.01;
    double constexpr t_retarget = 2.0;
    double t = 0.0;
    for (; t < t_retarget - 0.5 * dt; t += dt) {
        planner.sample(t);
    }

    // Confirm the planner is cruising at v_max with zero acceleration before the
    // retarget, so the "no velocity rise afterward" check is clean.
    auto const at_retarget = planner.sample(t_retarget);
    REQUIRE_THAT(at_retarget.velocity[0], WithinAbs(v_max, 1e-6));
    REQUIRE_THAT(at_retarget.acceleration[0], WithinAbs(0.0, 1e-6));

    // Retarget farther in the same direction (non-overshoot, same sign).
    planner.update(200.0);

    // Immediately after the retarget the planner should keep cruising, not brake.
    auto const shortly_after = planner.sample(t_retarget + 1.0);
    REQUIRE(shortly_after.velocity[0] > 0.9 * v_max);

    // After retarget-at-cruise the velocity is monotonically non-increasing (flat
    // cruise, then a single deceleration to rest), with bounded acceleration.
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

// -- Test 12: create rejects out-of-domain limits ----------------------------
//
// Both limits divide in the planner math (stopping distance v^2/(2*a_max),
// phase durations v_v/a_max), so the domain of each is finite and strictly
// positive; everything else is rejected with the limit-specific enumerator.
TEST_CASE("OnlinePlanner2nd: create rejects invalid limits",
          "[traj][online_planner_2nd][negative]")
{
    using planner_t = ctrlpp::online_planner_2nd<double>;
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    SECTION("invalid v_max -> non_positive_velocity_limit")
    {
        for (double const v_max : {0.0, -1.0, nan, inf}) {
            auto const result = planner_t::create({.v_max = v_max, .a_max = 10.0});
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error()
                    == ctrlpp::trajectory_error::non_positive_velocity_limit);
        }
    }

    SECTION("invalid a_max -> non_positive_acceleration_limit")
    {
        for (double const a_max : {0.0, -1.0, nan, inf}) {
            auto const result = planner_t::create({.v_max = 5.0, .a_max = a_max});
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error()
                    == ctrlpp::trajectory_error::non_positive_acceleration_limit);
        }
    }
}
