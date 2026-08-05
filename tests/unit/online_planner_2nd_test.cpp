#include "ctrlpp/trajectory/online_planner_2nd.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
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

    REQUIRE(planner.update(10.0).has_value());

    // Sample far enough in the future that the planner should have settled
    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(10.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));
    REQUIRE(planner.is_settled());
}

TEST_CASE("OnlinePlanner2nd: non-finite targets are rejected without mutation",
          "[traj][online_planner_2nd][negative]")
{
    auto planner = make_planner<double>({.v_max = 1.0, .a_max = 1.0});

    auto const rejected =
        planner.update(std::numeric_limits<double>::quiet_NaN());
    REQUIRE_FALSE(rejected.has_value());
    CHECK(rejected.error() == ctrlpp::trajectory_error::non_finite_input);

    auto const unchanged = planner.sample(0.0);
    CHECK(unchanged.position(0) == 0.0);
    CHECK(unchanged.velocity(0) == 0.0);
    CHECK(unchanged.acceleration(0) == 0.0);
    CHECK(planner.is_settled());

    REQUIRE(planner.update(1.0).has_value());
    CHECK(planner.sample(10.0).position(0) == 1.0);
}

TEST_CASE("OnlinePlanner2nd: nonzero sub-picometer moves remain bounded",
          "[traj][online_planner_2nd][precision]")
{
    auto planner = make_planner<double>({.v_max = 1.0, .a_max = 1.0});
    constexpr double target = 5e-13;

    REQUIRE(planner.update(target).has_value());
    auto const duration = planner.diagnostics().planned_duration;
    REQUIRE(duration > 0.0);

    auto const start = planner.sample(0.0);
    CHECK(start.position(0) == 0.0);
    auto const middle = planner.sample(duration / 2.0);
    CHECK(middle.position(0) > 0.0);
    CHECK(middle.position(0) < target);
    CHECK(std::abs(middle.velocity(0)) <= 1.0);
    CHECK(std::abs(middle.acceleration(0)) <= 1.0);

    auto const end = planner.sample(duration);
    CHECK(end.position(0) == target);
    CHECK(end.velocity(0) == 0.0);
}

// -- Test 2: Velocity never exceeds v_max --------------------------------------
TEST_CASE("OnlinePlanner2nd: velocity constraint", "[traj][online_planner_2nd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

    REQUIRE(planner.update(10.0).has_value());

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

    REQUIRE(planner.update(10.0).has_value());

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

    REQUIRE(planner.update(10.0).has_value());

    // Sample partway to build up state
    auto const mid = planner.sample(0.5);
    REQUIRE(mid.position[0] > 0.0);

    // Change target mid-motion
    REQUIRE(planner.update(5.0).has_value());

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

    REQUIRE(planner.update(-5.0).has_value());

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

    REQUIRE(planner.update(0.0).has_value());

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

    REQUIRE(planner.update(10.0).has_value());

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

    REQUIRE(planner.update(10.0).has_value());
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

    REQUIRE(planner.update(10.0).has_value());

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

    REQUIRE(planner.update(10.0f).has_value());
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
    REQUIRE(planner.update(100.0).has_value());

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
    REQUIRE(planner.update(200.0).has_value());

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

// -- Test 13: the brake-then-replan substitution is REPORTED --------------------
//
// This planner already stored the outcome of its own decision in members; what
// it never did was let the caller read it. The motion is limit-respecting either
// way, so the assertions below are on the report and the motion is checked
// second.
TEST_CASE("OnlinePlanner2nd: brake-then-replan substitution is reported",
          "[traj][online_planner_2nd][diagnostics]")
{
    double constexpr v_max = 5.0;
    double constexpr a_max = 10.0;

    auto planner = make_planner<double>({.v_max = v_max, .a_max = a_max});

    SECTION("a target behind the motion reports a reversal")
    {
        REQUIRE(planner.update(100.0).has_value());

        double t = 0.0;
        for (int i = 0; i < 200; ++i) {
            t = 0.01 * static_cast<double>(i);
            planner.sample(t);
        }
        auto const cruising = planner.sample(t);
        REQUIRE_THAT(cruising.velocity[0], WithinAbs(v_max, 1e-9));

        double const q0 = cruising.position[0];
        double const target = q0 - 10.0;
        CAPTURE(q0, cruising.velocity[0], target);

        REQUIRE(planner.update(target).has_value());

        auto const& diag = planner.diagnostics();
        REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::braked_and_replanned);
        REQUIRE(diag.substitution_reason
                == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
        REQUIRE(diag.brake_duration > 0.0);
        REQUIRE(diag.commanded_target == target);
        REQUIRE(diag.initial_velocity == cruising.velocity[0]);

        // Commanded against planned: the replan starts at the stopping point,
        // which lies on the far side of the commanded start from the target.
        REQUIRE(diag.replan_start_position > q0);
        REQUIRE(diag.planned_duration > diag.brake_duration);

        // Braking from v0 at a_max sweeps v0^2 / (2 * a_max) before rest.
        double const stop_dist = cruising.velocity[0] * cruising.velocity[0] / (2.0 * a_max);
        REQUIRE_THAT(diag.replan_start_position, WithinAbs(q0 + stop_dist, 1e-9));

        double constexpr tol = 1e-6;
        double const t_end = t + diag.planned_duration;
        for (double s = t; s < t_end; s += 0.001) {
            auto const pt = planner.sample(s);
            REQUIRE(std::abs(pt.velocity[0]) <= v_max + tol);
        }

        auto const settled = planner.sample(t_end + 1.0);
        REQUIRE_THAT(settled.position[0], WithinAbs(target, 1e-6));
        REQUIRE_THAT(settled.velocity[0], WithinAbs(0.0, 1e-6));
    }

    SECTION("a target inside the stopping distance reports an overshoot")
    {
        REQUIRE(planner.update(100.0).has_value());

        double t = 0.0;
        for (int i = 0; i < 200; ++i) {
            t = 0.01 * static_cast<double>(i);
            planner.sample(t);
        }
        auto const cruising = planner.sample(t);

        double const q0 = cruising.position[0];
        double const v0 = cruising.velocity[0];
        // Half the stopping distance: the planner cannot come to rest at the
        // target without passing it, so it brakes and comes back.
        double const stop_dist = v0 * v0 / (2.0 * a_max);
        double const target = q0 + stop_dist / 2.0;
        CAPTURE(q0, v0, stop_dist, target);

        REQUIRE(planner.update(target).has_value());

        auto const& diag = planner.diagnostics();
        REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::braked_and_replanned);
        REQUIRE(diag.substitution_reason
                == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
        REQUIRE(diag.brake_duration > 0.0);
        REQUIRE(diag.replan_start_position > target);
    }
}

// -- Test 14: an unsubstituted plan reports itself as such ---------------------
TEST_CASE("OnlinePlanner2nd: an ordinary move reports the commanded profile",
          "[traj][online_planner_2nd][diagnostics]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});
    planner.reset(0.0);

    REQUIRE(planner.update(10.0).has_value());

    auto const& diag = planner.diagnostics();
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::commanded_profile);
    REQUIRE(diag.substitution_reason == ctrlpp::online_planner_substitution_reason::none);
    REQUIRE(diag.brake_duration == 0.0);
    REQUIRE(diag.replan_start_position == 0.0);
    REQUIRE(diag.commanded_target == 10.0);
    REQUIRE(diag.planned_duration > 0.0);
}

// -- Test 15: a move commanded from inside the settle tolerance reports settled -
TEST_CASE("OnlinePlanner2nd: a settled command reports a zero-duration plan",
          "[traj][online_planner_2nd][diagnostics]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});
    planner.reset(3.0);

    REQUIRE(planner.update(3.0).has_value());

    auto const& diag = planner.diagnostics();
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::settled);
    REQUIRE(diag.substitution_reason == ctrlpp::online_planner_substitution_reason::none);
    REQUIRE(diag.planned_duration == 0.0);
    REQUIRE(diag.brake_duration == 0.0);
    REQUIRE(diag.replan_start_position == 3.0);
}

// -- Test 16: the carry-velocity reason is unreachable from this planner -------
//
// The two planners share one diagnostics type, and one of its substitution
// reasons belongs to a shape only the jerk-bounded planner builds. That is
// recorded here rather than papered over with a second type: sweep the same
// command families that drive the third-order planner through all three of its
// branches, and the reason never appears.
TEST_CASE("OnlinePlanner2nd: never reports the carry-velocity reason",
          "[traj][online_planner_2nd][diagnostics]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

    for (int i = 0; i < 200; ++i) {
        double const t = 0.05 * static_cast<double>(i);
        // Alternate far, near, behind, and already-reached targets so every
        // branch of compute_profile is commanded at some point in the sweep.
        double const q = planner.sample(t).position[0];
        switch (i % 4) {
        case 0: REQUIRE(planner.update(q + 50.0).has_value()); break;
        case 1: REQUIRE(planner.update(q + 0.05).has_value()); break;
        case 2: REQUIRE(planner.update(q - 20.0).has_value()); break;
        default: REQUIRE(planner.update(q).has_value()); break;
        }
        REQUIRE(planner.diagnostics().substitution_reason
                != ctrlpp::online_planner_substitution_reason::carry_velocity_shape_unavailable);
    }
}

// -- Test 17: reset clears a substitution report -------------------------------
TEST_CASE("OnlinePlanner2nd: reset clears the substitution report",
          "[traj][online_planner_2nd][diagnostics]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0});

    REQUIRE(planner.update(100.0).has_value());
    double t = 0.0;
    for (int i = 0; i < 200; ++i) {
        t = 0.01 * static_cast<double>(i);
        planner.sample(t);
    }
    auto const cruising = planner.sample(t);
    REQUIRE(planner.update(cruising.position[0] - 10.0).has_value());
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::braked_and_replanned);

    planner.reset(0.0);

    auto const& diag = planner.diagnostics();
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::settled);
    REQUIRE(diag.substitution_reason == ctrlpp::online_planner_substitution_reason::none);
    REQUIRE(diag.brake_duration == 0.0);
    REQUIRE(diag.planned_duration == 0.0);
}

namespace
{

struct planner_limits
{
    double v_max;
    double a_max;
};

// Seven pairs spanning just over seven decades of velocity limit. The two limits
// climb together so every pair plans the same trapezoidal shape on the same time
// scale (v_max / a_max is constant along the ladder), which keeps the sweep a
// sweep of SCALE rather than a sweep of branches. Powers of two, so nothing in
// the input is rounded on the way in.
constexpr std::array<planner_limits, 7> bit_identity_limits{{
    {0x1p-10, 0x1p-9},
    {0x1p-6, 0x1p-5},
    {0x1p-2, 0x1p-1},
    {0x1p+2, 0x1p+3},
    {0x1p+6, 0x1p+7},
    {0x1p+10, 0x1p+11},
    {0x1p+14, 0x1p+15},
}};

}

// -- Test 18: omitting the settle tolerances reproduces the previous behavior ---
//
// The two tolerance fields carry, as defaults, the value the planner compared
// both residuals against when the distance was a literal inside sample(). The
// claim is that OMITTING them is bit-identical to naming them, for every choice
// of limits. That is why the limits are swept over seven decades instead of one
// pair being tested, and why the comparison is exact equality: a tolerance here
// would be testing a different and weaker statement.
//
// The 1e-9 written out below is deliberately a literal, and the only settle
// literal in this file. It is the pre-knob value; if a future change moves the
// defaults, the omitting planner and the naming planner part company on the
// settle ladder and this case fails, which is what it is for.
TEST_CASE("OnlinePlanner2nd: omitting the settle tolerances is bit-identical to naming them",
          "[traj][online_planner_2nd][settle_tolerance]")
{
    constexpr int drive_steps = 64;
    constexpr int retarget_steps = 64;
    constexpr int ladder_points = 21;
    constexpr int points_per_pair = drive_steps + retarget_steps + ladder_points;

    int compared = 0;

    for (auto const& lim : bit_identity_limits) {
        auto omitted = make_planner<double>({.v_max = lim.v_max, .a_max = lim.a_max});
        auto named = make_planner<double>({
            .v_max = lim.v_max,
            .a_max = lim.a_max,
            .position_settle_tol = 1e-9,
            .velocity_settle_tol = 1e-9,
        });

        // Eight cruise-lengths of travel: long enough that the profile reaches
        // its cruise phase at every point on the ladder.
        double const target = 8.0 * lim.v_max;

        REQUIRE(omitted.update(target).has_value());
        REQUIRE(named.update(target).has_value());

        for (int i = 0; i < drive_steps; ++i) {
            double const t = 12.0 * static_cast<double>(i) / drive_steps;
            auto const from_omitted = omitted.sample(t);
            auto const from_named = named.sample(t);
            CAPTURE(lim.v_max, i, t);
            REQUIRE(from_omitted.position(0) == from_named.position(0));
            REQUIRE(from_omitted.velocity(0) == from_named.velocity(0));
            REQUIRE(from_omitted.acceleration(0) == from_named.acceleration(0));
            REQUIRE(omitted.is_settled() == named.is_settled());
            ++compared;
        }

        // Retarget across the origin from wherever the drive left the planners:
        // a brake-and-replan exercises the braking branch a rest-to-rest move
        // does not reach.
        REQUIRE(omitted.update(-target).has_value());
        REQUIRE(named.update(-target).has_value());

        for (int i = 0; i < retarget_steps; ++i) {
            double const t = 12.0 + 24.0 * static_cast<double>(i) / retarget_steps;
            auto const from_omitted = omitted.sample(t);
            auto const from_named = named.sample(t);
            CAPTURE(lim.v_max, i, t);
            REQUIRE(from_omitted.position(0) == from_named.position(0));
            REQUIRE(from_omitted.velocity(0) == from_named.velocity(0));
            REQUIRE(from_omitted.acceleration(0) == from_named.acceleration(0));
            REQUIRE(omitted.is_settled() == named.is_settled());
            ++compared;
        }

        // The settle ladder. reset() followed by a single update to a small
        // offset puts the planner one sample away from a state whose position
        // residual is EXACTLY that offset and whose velocity is exactly zero,
        // because sampling at the reference time integrates no phase at all. The
        // offsets bracket the default by three decades either side, so a default
        // that moved separates the two planners here whatever the limits are.
        for (int i = 0; i < ladder_points; ++i) {
            double const offset = std::ldexp(1.0, -40 + i);
            omitted.reset(0.0);
            named.reset(0.0);
            REQUIRE(omitted.update(offset).has_value());
            REQUIRE(named.update(offset).has_value());

            auto const from_omitted = omitted.sample(0.0);
            auto const from_named = named.sample(0.0);
            CAPTURE(lim.v_max, i, offset);
            REQUIRE(from_omitted.position(0) == 0.0);
            REQUIRE(from_omitted.velocity(0) == 0.0);
            REQUIRE(from_omitted.position(0) == from_named.position(0));
            REQUIRE(from_omitted.velocity(0) == from_named.velocity(0));
            REQUIRE(from_omitted.acceleration(0) == from_named.acceleration(0));
            REQUIRE(omitted.is_settled() == named.is_settled());
            ++compared;
        }
    }

    // The exact number of compared states, not a bound no loop can fail.
    REQUIRE(compared
            == points_per_pair * static_cast<int>(bit_identity_limits.size()));
}

// -- Test 19: each settle tolerance governs its own dimension and no other -----
//
// Sampling a fixed time-to-go before the end of a move leaves both residuals
// nonzero and three decades apart: the position error falls as tau^2 and the
// speed as tau. That is the state a single shared tolerance cannot describe,
// which is why the fields are separate.
//
// The bracketing values are read off the sampled state rather than written as
// literals, so the case moves with the knob instead of duplicating it. There is
// no acceleration field to bracket: this planner carries no acceleration state
// to settle.
TEST_CASE("OnlinePlanner2nd: each settle tolerance governs its own dimension",
          "[traj][online_planner_2nd][settle_tolerance]")
{
    constexpr double v_max = 0x1p+0;
    constexpr double a_max = 0x1p+1;
    constexpr double target = 0x1p+3;
    constexpr double tau = 0x1p-9;  // time-to-go at the probe, inside the decel ramp

    double t_probe = 0.0;
    double residual_q = 0.0;
    double residual_v = 0.0;

    {
        auto probe = make_planner<double>({.v_max = v_max, .a_max = a_max});
        REQUIRE(probe.update(target).has_value());
        t_probe = probe.diagnostics().planned_duration - tau;
        REQUIRE(t_probe > 0.0);

        auto const pt = probe.sample(t_probe);
        residual_q = std::abs(pt.position(0) - target);
        residual_v = std::abs(pt.velocity(0));

        // Both residuals are outside the default, so the default-configured
        // planner is not settled here and each bracket below is a real move.
        CHECK_FALSE(probe.is_settled());
    }

    // CAPTURE base-ten exponents rather than the values: Catch2 stringifies a
    // double in fixed notation, which renders a 1e-6 position residual as "0.0"
    // and would tell a reader diagnosing a failure the one thing that is untrue
    // of it.
    double const residual_q_decades = std::log10(residual_q);
    double const residual_v_decades = std::log10(residual_v);
    CAPTURE(residual_q_decades, residual_v_decades);

    REQUIRE(residual_q > 0.0);
    REQUIRE(residual_v > 0.0);

    // Two quantities, two units, two magnitudes.
    REQUIRE(residual_q < residual_v);

    auto settled_under = [&](double pos_tol, double vel_tol) {
        auto planner = make_planner<double>({
            .v_max = v_max,
            .a_max = a_max,
            .position_settle_tol = pos_tol,
            .velocity_settle_tol = vel_tol,
        });
        REQUIRE(planner.update(target).has_value());
        auto const pt = planner.sample(t_probe);

        // The knob does not move the motion: the sampled state is the one the
        // default-configured probe produced, bit for bit.
        REQUIRE(std::abs(pt.position(0) - target) == residual_q);
        REQUIRE(std::abs(pt.velocity(0)) == residual_v);
        return planner.is_settled();
    };

    double const above_q = 2.0 * residual_q;
    double const above_v = 2.0 * residual_v;
    double const below_q = 0.5 * residual_q;
    double const below_v = 0.5 * residual_v;

    // Raising both past their residuals settles a state the default does not:
    // the position transition moved to a larger position error, and the velocity
    // transition moved with its own quantity.
    CHECK(settled_under(above_q, above_v));

    // One field below its residual withholds the verdict on its own, with the
    // other raised. Each dimension is therefore decided by its own field and the
    // other cannot rescue it.
    CHECK_FALSE(settled_under(below_q, above_v));
    CHECK_FALSE(settled_under(above_q, below_v));
}

// -- Test 20: the settle policy does not reach the profile computation ---------
//
// The planner asks two different questions about distance. sample() asks "is the
// motion done", which is the policy the fields above carry. The profile
// computation asks "is this command a numerical no-op", and compares the target
// against the current position by EXACT equality. A policy tolerance must not be
// able to answer the second, or it would decide what the planner is allowed to
// compute.
TEST_CASE("OnlinePlanner2nd: a wide settle tolerance still plans a far smaller move",
          "[traj][online_planner_2nd][settle_tolerance]")
{
    auto planner = make_planner<double>({
        .v_max = 0x1p+0,
        .a_max = 0x1p+1,
        .position_settle_tol = 0x1p+0,
        .velocity_settle_tol = 0x1p+0,
    });

    // Twelve decades below the tolerance the planner was given.
    constexpr double tiny_target = 0x1p-40;

    REQUIRE(planner.update(tiny_target).has_value());

    auto const& diag = planner.diagnostics();
    CHECK(diag.disposition == ctrlpp::online_planner_disposition::commanded_profile);
    CHECK(diag.planned_duration > 0.0);

    // The move is realized in full rather than swallowed by the tolerance.
    CHECK(planner.sample(diag.planned_duration).position(0) == tiny_target);
}

// -- Test 21: the numerical-no-op velocity floor is derived, and it scales ------
//
// A separate question from the settle policy above, with a separate owner. The
// policy says when the application considers the axis arrived; this says when
// the arithmetic can no longer tell the commanded state from the one already in
// force. This planner carries no acceleration state, so its short-circuit is
// decided by two clauses alone: an EXACT comparison of the commanded
// displacement against zero, and the speed against one unit in the last place at
// the scale of the velocity limit. The floor moves with the limit rather than
// sitting at an absolute speed borrowed from an axis nobody named.
//
// The count is one because this planner performs no arithmetic on the snapshot
// before the test. That is why it is not the jerk-limited planner's seven: the
// chains differ, and copying that number across would be an unexplained constant
// wearing a different name.
TEST_CASE("OnlinePlanner2nd: the numerical-no-op velocity floor scales with the limit",
          "[traj][online_planner_2nd][no_op_floor]")
{
    // Three limit sets spanning nearly three decades, ratios held so the sweep is
    // a sweep of scale and not of branch structure.
    constexpr std::array<std::array<double, 2>, 3> limit_sets{{
        {5.0, 10.0},
        {0.25, 0.5},
        {40.0, 200.0},
    }};

    constexpr int command_velocity_rounding_ops = 1;

    int straddles = 0;
    for (auto const& limits : limit_sets) {
        double const v_max = limits[0];
        double const a_max = limits[1];
        double const velocity_floor = static_cast<double>(command_velocity_rounding_ops)
                                      * std::numeric_limits<double>::epsilon() * v_max;

        for (double const half_or_double : {0.5, 2.0}) {
            auto planner = make_planner<double>({.v_max = v_max, .a_max = a_max});
            planner.reset(0.0);
            REQUIRE(planner.update(100.0 * v_max).has_value());

            // A time s into the opening acceleration ramp leaves a speed a*s and
            // a position a*s^2/2, so s places the speed at half the floor and at
            // twice it.
            double const s = half_or_double * velocity_floor / a_max;
            auto const probe = planner.sample(s);

            // Commanding the position just sampled makes the displacement clause
            // an exact zero, so the speed clause is what decides.
            REQUIRE(planner.update(probe.position[0]).has_value());

            CAPTURE(v_max, a_max, half_or_double, std::log10(std::abs(probe.velocity[0])),
                    std::log10(velocity_floor), std::log10(std::abs(probe.position[0])));

            if (half_or_double < 1.0) {
                REQUIRE(std::abs(probe.velocity[0]) < velocity_floor);
                REQUIRE(planner.diagnostics().disposition
                        == ctrlpp::online_planner_disposition::settled);
                REQUIRE(planner.diagnostics().planned_duration == 0.0);
                REQUIRE(planner.is_settled());
            } else {
                REQUIRE(std::abs(probe.velocity[0]) > velocity_floor);
                REQUIRE(planner.diagnostics().disposition
                        != ctrlpp::online_planner_disposition::settled);
                REQUIRE_FALSE(planner.is_settled());
            }
        }
        ++straddles;
    }

    // An exact count, not a bound the loop cannot fail.
    REQUIRE(straddles == 3);
}

// -- Test 22: the overshoot verdict is relative, so it is scale-invariant -------
//
// The overshoot test compares two lengths the planner has already computed --
// the stopping distance and the remaining distance -- with a slack that is a
// counted multiple of epsilon rather than a distance. Its verdict is therefore a
// property of the RATIO of the two lengths and nothing else, and the same
// relative shortfall must be called an overshoot on an axis whose stopping
// distance is metres and on one whose stopping distance is picometres.
//
// The lower rungs of the ladder are where an absolute slack could not follow: at
// the smallest limit set the shortfall the axis is asked to resolve is five
// decades below the absolute slack this comparison used to carry, so that
// comparison would have reported no overshoot there.
TEST_CASE("OnlinePlanner2nd: the overshoot verdict is relative and scale-invariant",
          "[traj][online_planner_2nd][overshoot]")
{
    constexpr std::array<std::array<double, 2>, 4> limit_sets{{
        {5.0, 10.0},
        {1e-2, 5.0},
        {1e-4, 5.0},
        {1e-5, 5.0},
    }};

    // A relative shortfall many decades above the counted slack, so what is being
    // tested is the relativity and not the width.
    constexpr double relative_shortfall = 1e-6;

    int rungs = 0;
    for (auto const& limits : limit_sets) {
        double const v_max = limits[0];
        double const a_max = limits[1];

        auto planner = make_planner<double>({.v_max = v_max, .a_max = a_max});
        planner.reset(0.0);
        REQUIRE(planner.update(1000.0 * v_max).has_value());

        double const dt = 4.0 * v_max / a_max / 100.0;
        double t = 0.0;
        for (int i = 0; i < 400; ++i) {
            t = dt * static_cast<double>(i);
            planner.sample(t);
        }
        auto const cruising = planner.sample(t);
        double const q0 = cruising.position[0];
        double const v0 = cruising.velocity[0];

        // The planner's own stopping-distance expression, spelled the same way.
        double const stop_dist = v0 * v0 / (2.0 * a_max);
        REQUIRE(stop_dist > 0.0);

        REQUIRE(planner.update(q0 + stop_dist * (1.0 - relative_shortfall)).has_value());
        auto const& diag = planner.diagnostics();

        CAPTURE(v_max, std::log10(v0), std::log10(stop_dist),
                std::log10(stop_dist * relative_shortfall));

        REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::braked_and_replanned);
        REQUIRE(diag.substitution_reason
                == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
        ++rungs;
    }
    REQUIRE(rungs == 4);
}
