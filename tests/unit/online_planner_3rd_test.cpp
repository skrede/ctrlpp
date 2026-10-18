#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <bit>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <cstdint>
#include <utility>
#include <type_traits>

using Catch::Matchers::WithinAbs;

namespace
{

// Construction is fallible, so every valid-input site goes through the factory
// and asserts success here; a config that becomes unrealizable fails the test
// instead of quietly skipping it. Rejection cases never use this helper.
template <typename Scalar>
auto make_planner(typename ctrlpp::online_planner_3rd<Scalar>::config const& cfg)
    -> ctrlpp::online_planner_3rd<Scalar>
{
    auto created = ctrlpp::online_planner_3rd<Scalar>::create(cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

// Signed distance a jerk-limited stop from speed v (with zero acceleration)
// sweeps: a jerk ramp into deceleration, an optional constant-deceleration
// stretch once a_max is reached, and a jerk ramp back out of it. The same phase
// algebra the planner brakes with, so the two agree to rounding.
//
// @cite biagiotti2009 -- Sec. 3.4.3
auto jerk_limited_stop_distance(double v, double a_max, double j_max) -> double
{
    double const abs_v = std::abs(v);
    double const sign_v = (v >= 0.0) ? 1.0 : -1.0;

    if (abs_v <= a_max * a_max / j_max) {
        // Triangular deceleration: two jerk ramps, no constant stretch. The mean
        // speed over the symmetric ramp pair is abs_v / 2 across 2 * T_j.
        double const T_j = std::sqrt(abs_v / j_max);
        return sign_v * abs_v * T_j;
    }

    double const T_j = a_max / j_max;
    double const T_c = abs_v / a_max - T_j;
    double const v1 = abs_v - j_max * T_j * T_j / 2.0;
    double const q1 = abs_v * T_j - j_max * T_j * T_j * T_j / 6.0;
    double const v2 = v1 - a_max * T_c;
    double const q2 = q1 + v1 * T_c - a_max * T_c * T_c / 2.0;
    double const q3 = q2 + v2 * T_j - a_max * T_j * T_j / 2.0 + j_max * T_j * T_j * T_j / 6.0;
    return sign_v * q3;
}

// Drives the planner up to its cruise velocity on a far target and returns the
// state it is cruising at. Acceleration is zero there, so the planner's next
// update plans straight from this state with no acceleration-nulling phase in
// between and the numbers below are the ones it actually reasons about.
struct cruising_state
{
    double t{};
    double q{};
    double v{};
};

auto cruise_up(ctrlpp::online_planner_3rd<double>& planner, double far_target,
               double dt, int steps) -> cruising_state
{
    REQUIRE(planner.update(far_target).has_value());

    double t = 0.0;
    for (int i = 0; i < steps; ++i) {
        t = dt * static_cast<double>(i);
        planner.sample(t);
    }

    auto const pt = planner.sample(t);
    return {t, pt.position[0], pt.velocity[0]};
}

}

// -- Test 1: Step response settles to target -----------------------------------
TEST_CASE("OnlinePlanner3rd: step response settles", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.update(10.0).has_value());

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(10.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(pt.acceleration[0], WithinAbs(0.0, 1e-6));
    REQUIRE(planner.is_settled());
}

TEST_CASE("OnlinePlanner3rd: non-finite targets are rejected without mutation",
          "[trajectory][online_planner_3rd][negative]")
{
    auto planner =
        make_planner<double>({.v_max = 1.0, .a_max = 1.0, .j_max = 1.0});

    auto const rejected =
        planner.update(std::numeric_limits<double>::infinity());
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

TEST_CASE("OnlinePlanner3rd: nonzero sub-picometer moves remain bounded",
          "[trajectory][online_planner_3rd][precision]")
{
    auto planner =
        make_planner<double>({.v_max = 1.0, .a_max = 1.0, .j_max = 1.0});
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
    CHECK(end.acceleration(0) == 0.0);
}

// -- Test 2: Velocity never exceeds v_max --------------------------------------
TEST_CASE("OnlinePlanner3rd: velocity constraint", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.update(10.0).has_value());

    double constexpr tol = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = 10.0 * static_cast<double>(i) / 1000.0;
        auto const pt = planner.sample(t);
        REQUIRE(std::abs(pt.velocity[0]) <= 5.0 + tol);
    }
}

// -- Test 3: Acceleration never exceeds a_max ----------------------------------
TEST_CASE("OnlinePlanner3rd: acceleration constraint", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.update(10.0).has_value());

    double constexpr tol = 1e-6;
    for (int i = 0; i <= 1000; ++i) {
        double const t = 10.0 * static_cast<double>(i) / 1000.0;
        auto const pt = planner.sample(t);
        REQUIRE(std::abs(pt.acceleration[0]) <= 10.0 + tol);
    }
}

// -- Test 4: Jerk never exceeds j_max (numerical differentiation) ---------------
TEST_CASE("OnlinePlanner3rd: jerk constraint", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.update(10.0).has_value());

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
TEST_CASE("OnlinePlanner3rd: mid-motion target change", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.update(10.0).has_value());

    // Sample partway
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

// -- Test 6: Negative displacement ---------------------------------------------
TEST_CASE("OnlinePlanner3rd: negative displacement", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.update(-5.0).has_value());

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(-5.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));
}

// -- Test 7: is_settled returns correct state ----------------------------------
TEST_CASE("OnlinePlanner3rd: is_settled transitions", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.is_settled());

    REQUIRE(planner.update(10.0).has_value());
    planner.sample(0.1);
    REQUIRE_FALSE(planner.is_settled());

    planner.sample(100.0);
    REQUIRE(planner.is_settled());
}

// -- Test 8: Reset functionality -----------------------------------------------
TEST_CASE("OnlinePlanner3rd: reset", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.update(10.0).has_value());
    planner.sample(100.0);

    planner.reset(3.0);
    REQUIRE(planner.is_settled());

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(3.0, 1e-12));
}

// -- Test 9: Profile has smooth acceleration (double-S shape) ------------------
TEST_CASE("OnlinePlanner3rd: smooth acceleration profile", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.update(10.0).has_value());

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
TEST_CASE("OnlinePlanner3rd: float type", "[trajectory][online_planner_3rd]")
{
    auto planner = make_planner<float>({.v_max = 5.0f, .a_max = 10.0f, .j_max = 50.0f});

    REQUIRE(planner.update(10.0f).has_value());
    auto const pt = planner.sample(100.0f);
    REQUIRE_THAT(static_cast<double>(pt.position[0]), WithinAbs(10.0, 1e-3));
}

// -- Test 11: Same-direction retarget carries velocity through (no full-stop dip)
//
// Retargeting farther in the same direction while cruising must NOT brake the
// motion to a full stop and re-accelerate. A time-optimal planner keeps the
// current velocity and simply extends the profile, so once the planner is
// cruising the velocity never rises again after the retarget.
TEST_CASE("OnlinePlanner3rd: same-direction retarget keeps velocity",
          "[trajectory][online_planner_3rd]")
{
    double constexpr v_max = 5.0;
    auto planner = make_planner<double>({.v_max = v_max, .a_max = 10.0, .j_max = 50.0});

    // Command a far target and cruise up to v_max.
    REQUIRE(planner.update(100.0).has_value());

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
    REQUIRE(planner.update(200.0).has_value());

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

// -- Test 12: create rejects out-of-domain limits ----------------------------
//
// All three limits divide in the planner math (cruise duration h/v_max,
// jerk-phase durations a_max/j_max and |a|/j_max), so the domain of each is
// finite and strictly positive; everything else is rejected with the
// limit-specific enumerator.
TEST_CASE("OnlinePlanner3rd: create rejects invalid limits",
          "[trajectory][online_planner_3rd][negative]")
{
    using planner_t = ctrlpp::online_planner_3rd<double>;
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    SECTION("invalid v_max -> non_positive_velocity_limit")
    {
        for (double const v_max : {0.0, -1.0, nan, inf}) {
            auto const result =
                planner_t::create({.v_max = v_max, .a_max = 10.0, .j_max = 50.0});
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error()
                    == ctrlpp::trajectory_error::non_positive_velocity_limit);
        }
    }

    SECTION("invalid a_max -> non_positive_acceleration_limit")
    {
        for (double const a_max : {0.0, -1.0, nan, inf}) {
            auto const result =
                planner_t::create({.v_max = 5.0, .a_max = a_max, .j_max = 50.0});
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error()
                    == ctrlpp::trajectory_error::non_positive_acceleration_limit);
        }
    }

    SECTION("invalid j_max -> non_positive_jerk_limit")
    {
        for (double const j_max : {0.0, -1.0, nan, inf}) {
            auto const result =
                planner_t::create({.v_max = 5.0, .a_max = 10.0, .j_max = j_max});
            REQUIRE_FALSE(result.has_value());
            REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_jerk_limit);
        }
    }
}

// -- Test 13: the overshoot verdict and the carry-velocity shape agree ----------
//
// The planner brakes to rest when the remaining distance is shorter than the
// carried velocity sweeps on its way to rest, and the carry-velocity double-S
// refuses exactly that command. The boundary is located from the shape's own
// verdict, and every probe around it is a fresh planner driven to the same
// cruise, so the report at each target can be held against the shape's answer
// for the same operands. A target the shape refuses must be reported as a
// predicted overshoot, never as an unpredicted shape refusal.
TEST_CASE("OnlinePlanner3rd: the overshoot verdict agrees with the carry-velocity shape at its boundary",
          "[trajectory][online_planner_3rd][diagnostics]")
{
    double constexpr v_max = 5.0;
    double constexpr a_max = 10.0;
    double constexpr j_max = 50.0;
    int constexpr probes_each_side = 8;

    auto const cruise = [&](ctrlpp::online_planner_3rd<double>& planner) {
        return cruise_up(planner, 100.0, 0.01, 200);
    };
    auto reference = make_planner<double>({.v_max = v_max, .a_max = a_max, .j_max = j_max});
    auto const cruising = cruise(reference);
    REQUIRE(cruising.v > 0.0);
    REQUIRE(cruising.v <= v_max);

    auto const shape = [&](double target) {
        return ctrlpp::double_s_trajectory<double>::create({
            .q0 = cruising.q,
            .q1 = target,
            .v_max = v_max,
            .a_max = a_max,
            .j_max = j_max,
            .v0 = cruising.v,
            .v1 = 0.0,
        });
    };
    auto const commanded = [&](double target) {
        auto planner = make_planner<double>({.v_max = v_max, .a_max = a_max, .j_max = j_max});
        auto const state = cruise(planner);
        REQUIRE(state.q == cruising.q);
        REQUIRE(state.v == cruising.v);
        REQUIRE(planner.update(target).has_value());
        return planner;
    };

    double const stop_dist = jerk_limited_stop_distance(cruising.v, a_max, j_max);
    double refused = cruising.q;
    double admitted = cruising.q + 2.0 * stop_dist;
    REQUIRE_FALSE(shape(refused).has_value());
    REQUIRE(shape(admitted).has_value());
    while (std::nextafter(refused, admitted) != admitted) {
        double const mid = refused + (admitted - refused) / 2.0;
        (shape(mid).has_value() ? admitted : refused) = mid;
    }
    CAPTURE(cruising.q, cruising.v, stop_dist, refused, admitted);

    auto const refusal = shape(refused);
    REQUIRE_FALSE(refusal.has_value());
    REQUIRE(refusal.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity);
    CHECK(commanded(admitted).diagnostics().substitution_reason
          == ctrlpp::online_planner_substitution_reason::none);
    CHECK(commanded(refused).diagnostics().substitution_reason
          == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);

    double target = admitted;
    for (int i = 0; i < probes_each_side; ++i) {
        target = std::nextafter(target, cruising.q);
    }
    for (int i = 0; i <= 2 * probes_each_side; ++i, target = std::nextafter(target, admitted + stop_dist)) {
        auto const reason = commanded(target).diagnostics().substitution_reason;
        auto const expected = shape(target).has_value()
                                  ? ctrlpp::online_planner_substitution_reason::none
                                  : ctrlpp::online_planner_substitution_reason::reversal_or_overshoot;
        CAPTURE(i, target);
        CHECK(reason == expected);
        CHECK(reason != ctrlpp::online_planner_substitution_reason::carry_velocity_shape_unavailable);
    }

    auto planner = commanded(refused);
    auto const& diag = planner.diagnostics();
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(diag.brake_duration > 0.0);
    REQUIRE(diag.commanded_target == refused);
    REQUIRE(diag.initial_velocity == cruising.v);
    REQUIRE(diag.replan_start_position != cruising.q);
    REQUIRE_THAT(diag.replan_start_position, WithinAbs(cruising.q + stop_dist, 1e-9));
    REQUIRE(diag.planned_duration > diag.brake_duration);
    REQUIRE(std::isfinite(diag.planned_duration));

    // The motion runs from the commanded start to the stopping point, never
    // turning back: the interval it stays inside is the one the report names.
    double constexpr tol = 1e-6;
    double const t_end = cruising.t + diag.planned_duration;
    for (double t = cruising.t; t < t_end; t += 0.001) {
        auto const pt = planner.sample(t);
        REQUIRE(std::abs(pt.velocity[0]) <= v_max + tol);
        REQUIRE(std::abs(pt.acceleration[0]) <= a_max + tol);
        REQUIRE(pt.position[0] >= cruising.q - tol);
        REQUIRE(pt.position[0] <= diag.replan_start_position + tol);
    }

    auto const settled = planner.sample(t_end + 1.0);
    REQUIRE_THAT(settled.position[0], WithinAbs(refused, 1e-6));
    REQUIRE_THAT(settled.velocity[0], WithinAbs(0.0, 1e-6));
}

// -- Test 14: the reversal substitution is reported, and is DISTINGUISHABLE -----
//
// Same disposition, different cause. A caller that only learned "something was
// substituted" could not choose between axes or explain the extra time; the
// reason enumerator is what makes the two causes separable.
TEST_CASE("OnlinePlanner3rd: reversal substitution reports a distinct reason",
          "[trajectory][online_planner_3rd][diagnostics]")
{
    double constexpr v_max = 5.0;
    double constexpr a_max = 10.0;
    double constexpr j_max = 50.0;

    auto planner = make_planner<double>({.v_max = v_max, .a_max = a_max, .j_max = j_max});

    auto const cruising = cruise_up(planner, 100.0, 0.01, 200);
    REQUIRE(cruising.v > 0.0);

    // Command a target behind the planner while it moves forward at speed.
    double const target = cruising.q - 10.0;
    CAPTURE(cruising.q, cruising.v, target);

    REQUIRE(planner.update(target).has_value());

    auto const& diag = planner.diagnostics();
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(diag.substitution_reason
            == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
    REQUIRE(diag.substitution_reason
            != ctrlpp::online_planner_substitution_reason::carry_velocity_shape_unavailable);
    REQUIRE(diag.brake_duration > 0.0);

    // Here the plan is strictly longer than its braking: the stopping point lies
    // on the far side of the start from the target, so a real move follows.
    REQUIRE(diag.planned_duration > diag.brake_duration);
    REQUIRE(diag.replan_start_position > cruising.q);
    REQUIRE(diag.replan_start_position > target);

    double constexpr tol = 1e-6;
    double const t_end = cruising.t + diag.planned_duration;
    for (double t = cruising.t; t < t_end; t += 0.001) {
        auto const pt = planner.sample(t);
        REQUIRE(std::abs(pt.velocity[0]) <= v_max + tol);
        REQUIRE(std::abs(pt.acceleration[0]) <= a_max + tol);
    }

    auto const settled = planner.sample(t_end + 1.0);
    REQUIRE_THAT(settled.position[0], WithinAbs(target, 1e-6));
    REQUIRE_THAT(settled.velocity[0], WithinAbs(0.0, 1e-6));
}

// -- Test 15: an unsubstituted plan reports itself as such ---------------------
TEST_CASE("OnlinePlanner3rd: an ordinary move reports the commanded profile",
          "[trajectory][online_planner_3rd][diagnostics]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});
    planner.reset(0.0);

    REQUIRE(planner.update(10.0).has_value());

    auto const& diag = planner.diagnostics();
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::commanded_profile);
    REQUIRE(diag.substitution_reason == ctrlpp::online_planner_substitution_reason::none);
    REQUIRE(diag.brake_duration == 0.0);
    REQUIRE(diag.replan_start_position == 0.0);
    REQUIRE(diag.commanded_target == 10.0);
    REQUIRE(diag.planned_duration > 0.0);

    // A same-direction retarget from cruise keeps the carry-velocity shape, so
    // it too is the commanded profile and not a substitution.
    auto const cruising = cruise_up(planner, 100.0, 0.01, 200);
    REQUIRE(planner.update(200.0).has_value());
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::commanded_profile);
    REQUIRE(diag.substitution_reason == ctrlpp::online_planner_substitution_reason::none);
    REQUIRE(diag.brake_duration == 0.0);
    REQUIRE(diag.initial_velocity == cruising.v);
}

// -- Test 16: a command at the target from rest reports settled -----------------
TEST_CASE("OnlinePlanner3rd: a settled command reports a zero-duration plan",
          "[trajectory][online_planner_3rd][diagnostics]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});
    planner.reset(3.0);

    REQUIRE(planner.update(3.0).has_value());

    auto const& diag = planner.diagnostics();
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::settled);
    REQUIRE(diag.substitution_reason == ctrlpp::online_planner_substitution_reason::none);
    REQUIRE(diag.planned_duration == 0.0);
    REQUIRE(diag.brake_duration == 0.0);
    REQUIRE(diag.commanded_target == 3.0);
    REQUIRE(diag.replan_start_position == 3.0);
}

// -- Test 17: reset clears a substitution report -------------------------------
//
// A report that outlived the plan it describes would be the same defect in a
// different place.
TEST_CASE("OnlinePlanner3rd: reset clears the substitution report",
          "[trajectory][online_planner_3rd][diagnostics]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    auto const cruising = cruise_up(planner, 100.0, 0.01, 200);
    REQUIRE(planner.update(cruising.q - 10.0).has_value());
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
    double j_max;
};

// Seven triples spanning just over seven decades of velocity limit. The three
// limits climb together so every triple plans the same seven-phase shape on the
// same time scale (a_max / j_max and v_max / a_max are constant along the
// ladder), which keeps the sweep a sweep of SCALE rather than a sweep of
// branches. Powers of two, so nothing in the input is rounded on the way in.
constexpr std::array<planner_limits, 7> bit_identity_limits{{
    {0x1p-10, 0x1p-9, 0x1p-8},
    {0x1p-6, 0x1p-5, 0x1p-4},
    {0x1p-2, 0x1p-1, 0x1p+0},
    {0x1p+2, 0x1p+3, 0x1p+4},
    {0x1p+6, 0x1p+7, 0x1p+8},
    {0x1p+10, 0x1p+11, 0x1p+12},
    {0x1p+14, 0x1p+15, 0x1p+16},
}};

}

// -- Test 18: omitting the settle tolerances reproduces the previous behavior ---
//
// The three tolerance fields carry, as defaults, the value the planner compared
// all three residuals against when the distance was a literal inside sample().
// The claim is that OMITTING them is bit-identical to naming them, for every
// choice of limits. That is why the limits are swept over seven decades instead
// of one triple being tested, and why the comparison is exact equality: a
// tolerance here would be testing a different and weaker statement.
//
// The 1e-9 written out below is deliberately a literal, and the only settle
// literal in this file. It is the pre-knob value; if a future change moves the
// defaults, the omitting planner and the naming planner part company on the
// settle ladder and this case fails, which is what it is for.
TEST_CASE("OnlinePlanner3rd: omitting the settle tolerances is bit-identical to naming them",
          "[trajectory][online_planner_3rd][settle_tolerance]")
{
    constexpr int drive_steps = 64;
    constexpr int retarget_steps = 64;
    constexpr int ladder_points = 21;
    constexpr int points_per_triple = drive_steps + retarget_steps + ladder_points;

    int compared = 0;

    for (auto const& lim : bit_identity_limits) {
        auto omitted = make_planner<double>(
            {.v_max = lim.v_max, .a_max = lim.a_max, .j_max = lim.j_max});
        auto named = make_planner<double>({
            .v_max = lim.v_max,
            .a_max = lim.a_max,
            .j_max = lim.j_max,
            .position_settle_tol = 1e-9,
            .velocity_settle_tol = 1e-9,
            .acceleration_settle_tol = 1e-9,
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
        // a brake-and-replan exercises the phase array a rest-to-rest move does
        // not reach.
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
        // residual is EXACTLY that offset and whose velocity and acceleration
        // are exactly zero, because sampling at the reference time integrates no
        // phase at all. The offsets bracket the default by three decades either
        // side, so a default that moved separates the two planners here whatever
        // the limits are.
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
            REQUIRE(from_omitted.acceleration(0) == 0.0);
            REQUIRE(from_omitted.position(0) == from_named.position(0));
            REQUIRE(from_omitted.velocity(0) == from_named.velocity(0));
            REQUIRE(from_omitted.acceleration(0) == from_named.acceleration(0));
            REQUIRE(omitted.is_settled() == named.is_settled());
            ++compared;
        }
    }

    // The exact number of compared states, not a bound no loop can fail.
    REQUIRE(compared
            == points_per_triple * static_cast<int>(bit_identity_limits.size()));
}

// -- Test 19: each settle tolerance governs its own dimension and no other -----
//
// Sampling one jerk-ramp fraction before the end of a move leaves all three
// residuals nonzero and decades apart: the position error falls as tau^3, the
// speed as tau^2 and the acceleration as tau. That is the state a single shared
// tolerance cannot describe, which is why the fields are separate.
//
// The bracketing values are read off the sampled state rather than written as
// literals, so the case moves with the knob instead of duplicating it.
TEST_CASE("OnlinePlanner3rd: each settle tolerance governs its own dimension",
          "[trajectory][online_planner_3rd][settle_tolerance]")
{
    constexpr double v_max = 0x1p+0;
    constexpr double a_max = 0x1p+1;
    constexpr double j_max = 0x1p+2;
    constexpr double target = 0x1p+3;
    constexpr double tau = 0x1p-9;  // time-to-go at the probe, inside the last ramp

    double t_probe = 0.0;
    double residual_q = 0.0;
    double residual_v = 0.0;
    double residual_a = 0.0;

    {
        auto probe = make_planner<double>(
            {.v_max = v_max, .a_max = a_max, .j_max = j_max});
        REQUIRE(probe.update(target).has_value());
        t_probe = probe.diagnostics().planned_duration - tau;
        REQUIRE(t_probe > 0.0);

        auto const pt = probe.sample(t_probe);
        residual_q = std::abs(pt.position(0) - target);
        residual_v = std::abs(pt.velocity(0));
        residual_a = std::abs(pt.acceleration(0));

        // Every residual is outside the default, so the default-configured
        // planner is not settled here and each bracket below is a real move.
        CHECK_FALSE(probe.is_settled());
    }

    // CAPTURE base-ten exponents rather than the values: Catch2 stringifies a
    // double in fixed notation, which renders a 1e-9 position residual as "0.0"
    // and would tell a reader diagnosing a failure the one thing that is untrue
    // of it.
    double const residual_q_decades = std::log10(residual_q);
    double const residual_v_decades = std::log10(residual_v);
    double const residual_a_decades = std::log10(residual_a);
    CAPTURE(residual_q_decades, residual_v_decades, residual_a_decades);

    REQUIRE(residual_q > 0.0);
    REQUIRE(residual_v > 0.0);
    REQUIRE(residual_a > 0.0);

    // Three quantities, three units, three magnitudes.
    REQUIRE(residual_q < residual_v);
    REQUIRE(residual_v < residual_a);

    auto settled_under = [&](double pos_tol, double vel_tol, double acc_tol) {
        auto planner = make_planner<double>({
            .v_max = v_max,
            .a_max = a_max,
            .j_max = j_max,
            .position_settle_tol = pos_tol,
            .velocity_settle_tol = vel_tol,
            .acceleration_settle_tol = acc_tol,
        });
        REQUIRE(planner.update(target).has_value());
        auto const pt = planner.sample(t_probe);

        // The knob does not move the motion: the sampled state is the one the
        // default-configured probe produced, bit for bit.
        REQUIRE(std::abs(pt.position(0) - target) == residual_q);
        REQUIRE(std::abs(pt.velocity(0)) == residual_v);
        REQUIRE(std::abs(pt.acceleration(0)) == residual_a);
        return planner.is_settled();
    };

    double const above_q = 2.0 * residual_q;
    double const above_v = 2.0 * residual_v;
    double const above_a = 2.0 * residual_a;
    double const below_q = 0.5 * residual_q;
    double const below_v = 0.5 * residual_v;
    double const below_a = 0.5 * residual_a;

    // Raising all three past their residuals settles a state the default does
    // not: the position transition moved to a larger position error, and the
    // other two moved with their own quantities.
    CHECK(settled_under(above_q, above_v, above_a));

    // One field below its residual withholds the verdict on its own, with the
    // other two raised. Each dimension is therefore decided by its own field and
    // nothing else can rescue it.
    CHECK_FALSE(settled_under(below_q, above_v, above_a));
    CHECK_FALSE(settled_under(above_q, below_v, above_a));
    CHECK_FALSE(settled_under(above_q, above_v, below_a));
}

// -- Test 20: the settle policy does not reach the profile computation ---------
//
// The planner asks two different questions about distance. sample() asks "is
// the motion done", which is the policy the fields above carry. The profile
// computation asks "is this command a numerical no-op", and compares the target
// against the current position by EXACT equality. A policy tolerance must not be
// able to answer the second, or it would decide what the planner is allowed to
// compute.
TEST_CASE("OnlinePlanner3rd: a wide settle tolerance still plans a far smaller move",
          "[trajectory][online_planner_3rd][settle_tolerance]")
{
    auto planner = make_planner<double>({
        .v_max = 0x1p+0,
        .a_max = 0x1p+1,
        .j_max = 0x1p+2,
        .position_settle_tol = 0x1p+0,
        .velocity_settle_tol = 0x1p+0,
        .acceleration_settle_tol = 0x1p+0,
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

// -- Test 21: a command at the current position is a no-op only from rest ------
//
// Removing any nonzero acceleration takes |a| / j_max at the jerk limit, and the planner is never told the sample
// period, so no magnitude exists below which dropping it would respect j_max. A command at the position just sampled
// is therefore planned whenever the sampled velocity or acceleration is nonzero, however small.
//
// The accelerations are taken a time s into the opening jerk ramp, at half and twice eps * a_max; the speed j*s^2/2
// left there is decades smaller still. From cruise the acceleration is zero to one rounding and the speed is the
// cruise speed, so that state is planned too.
TEST_CASE("OnlinePlanner3rd: a command at the current position is a no-op only from rest",
          "[trajectory][online_planner_3rd][no_op_floor]")
{
    // Three limit sets spanning five decades of acceleration limit, with the
    // ratios held so the sweep is a sweep of scale and not of branch structure.
    constexpr std::array<std::array<double, 3>, 3> limit_sets{{
        {5.0, 10.0, 50.0},
        {0.25, 0.5, 4.0},
        {40.0, 200.0, 2000.0},
    }};

    int planned = 0;
    for (auto const& limits : limit_sets) {
        double const v_max = limits[0];
        double const a_max = limits[1];
        double const j_max = limits[2];
        double const unit_acceleration = std::numeric_limits<double>::epsilon() * a_max;

        for (double const half_or_double : {0.5, 2.0}) {
            auto planner = make_planner<double>(
                {.v_max = v_max, .a_max = a_max, .j_max = j_max});
            planner.reset(0.0);
            REQUIRE(planner.update(100.0 * v_max).has_value());

            double const s = half_or_double * unit_acceleration / j_max;
            auto const probe = planner.sample(s);
            REQUIRE(planner.update(probe.position[0]).has_value());

            CAPTURE(v_max, a_max, j_max, half_or_double, probe.acceleration[0], probe.velocity[0]);
            REQUIRE(probe.acceleration[0] != 0.0);
            REQUIRE(std::abs(probe.velocity[0]) < std::numeric_limits<double>::epsilon() * v_max);

            // The plan nulls the acceleration at the jerk limit before anything else, so it lasts at least that long.
            REQUIRE(planner.diagnostics().disposition
                    != ctrlpp::online_planner_disposition::settled);
            REQUIRE(planner.diagnostics().planned_duration >= std::abs(probe.acceleration[0]) / j_max);
            ++planned;
        }
    }
    REQUIRE(planned == 6);

    // The same from cruise. The two jerk ramps that built the cruise add and
    // remove the same product j*T, which cancels exactly when both are rounded,
    // and leaves the first one's rounding when a fused multiply-add removes the
    // exact product: within half a unit in the last place of j*T <= a_max.
    constexpr double cruise_acceleration_rounding_ops = 1.0;
    constexpr double unit_roundoff = 0.5 * std::numeric_limits<double>::epsilon();
    int cruises = 0;
    for (auto const& limits : limit_sets) {
        auto planner = make_planner<double>(
            {.v_max = limits[0], .a_max = limits[1], .j_max = limits[2]});
        planner.reset(0.0);
        auto const cruising = cruise_up(planner, 1000.0 * limits[0], 0.01, 400);
        CAPTURE(limits[0], cruising.q, cruising.v);

        REQUIRE(std::abs(planner.sample(cruising.t).acceleration[0])
                <= cruise_acceleration_rounding_ops * unit_roundoff * limits[1]);
        REQUIRE(std::abs(cruising.v)
                > std::numeric_limits<double>::epsilon() * limits[0]);

        REQUIRE(planner.update(cruising.q).has_value());
        REQUIRE(planner.diagnostics().disposition
                != ctrlpp::online_planner_disposition::settled);
        ++cruises;
    }
    REQUIRE(cruises == 3);

    // From rest the command is the one already in force.
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});
    planner.reset(0x1.8p+1);
    REQUIRE(planner.update(0x1.8p+1).has_value());
    REQUIRE(planner.diagnostics().disposition == ctrlpp::online_planner_disposition::settled);
    REQUIRE(planner.diagnostics().planned_duration == 0.0);
}

// -- Test 22: every nonzero starting acceleration is nulled in its own phase --
//
// The evaluation starts from the snapshot's acceleration, so the plan nulls every nonzero one in its own phase first,
// however short, and is planned from exactly the state it is evaluated from. Sampled at the phase's own duration
// |a_start| / j_max, the acceleration is zero. The starting accelerations straddle the one whose nulling phase moves
// the axis by five roundings of the stopping distance from full speed, so both are far below anything a phase that
// short changes on the velocity.
TEST_CASE("OnlinePlanner3rd: the acceleration-nulling phase is emitted for every nonzero starting acceleration",
          "[trajectory][online_planner_3rd][no_op_floor]")
{
    constexpr std::array<std::array<double, 3>, 3> limit_sets{{
        {5.0, 10.0, 50.0},
        {0.25, 0.5, 4.0},
        {40.0, 200.0, 2000.0},
    }};

    constexpr int displacement_rounding_ops = 5;

    int straddles = 0;
    for (auto const& limits : limit_sets) {
        double const v_max = limits[0];
        double const a_max = limits[1];
        double const j_max = limits[2];

        double const a_reference = static_cast<double>(displacement_rounding_ops)
                                   * std::numeric_limits<double>::epsilon() * v_max * j_max / (2.0 * a_max);

        for (double const below_or_above : {0.5, 4.0}) {
            auto planner = make_planner<double>(
                {.v_max = v_max, .a_max = a_max, .j_max = j_max});
            planner.reset(0.0);
            REQUIRE(planner.update(100.0 * v_max).has_value());

            double const s = below_or_above * a_reference / j_max;
            auto const probe = planner.sample(s);
            double const a_start = probe.acceleration[0];
            REQUIRE(a_start > 0.0);

            REQUIRE(planner.update(100.0 * v_max).has_value());
            double const nulling_duration = a_start / j_max;
            auto const after = planner.sample(s + nulling_duration);

            CAPTURE(v_max, a_max, j_max, below_or_above, std::log10(a_start),
                    std::log10(a_reference), std::log10(std::abs(after.acceleration[0]) + 1e-300));

            REQUIRE(std::abs(after.acceleration[0])
                    < std::numeric_limits<double>::epsilon() * a_max);
        }
        ++straddles;
    }
    REQUIRE(straddles == 3);
}

// A nulling duration |a| / j_max cannot underflow through the public API, because every sampled acceleration is
// j_max times a time. So the groups emptied here are a degenerate plan's cruise and a cruise retarget's ramp.
TEST_CASE("OnlinePlanner3rd: a group whose phases all vanish still re-bases the evaluation at its end",
          "[trajectory][online_planner_3rd]")
{
    // A zero-length phase adds v*0 + a*0*0/2 + j*0*0*0/6 = 0 to the position and re-bases (v, a) to the values the
    // previous pin already left, so an emptied group changes no sample. The degenerate command has h / (2 j_max) = 1,
    // whose cube root every libm returns exactly, so its jerk time is 1 and every speed, acceleration and sample time
    // below is a short dyadic rational: every product the evaluation forms is exact, and the rows are the same bits
    // whether or not a compiler fuses a multiply-add. Each speed and acceleration is the exact profile's, each position
    // within one rounding of it.
    using sample_row = std::array<double, 4>;
    constexpr std::array<sample_row, 4> degenerate{{
        {0x1.8p+0, 0x1.0aaaaaaaaaaabp-1, 0x1.cp-1, 0x1p-1},
        {0x1.4p+1, 0x1.7aaaaaaaaaaaap+0, 0x1.cp-1, -0x1p-1},
        {0x1.8p+1, 0x1.d555555555555p+0, 0x1p-1, -0x1p+0},
        {0x1.cp+1, 0x1.faaaaaaaaaaabp+0, 0x1p-3, -0x1p-1},
    }};
    constexpr std::array<sample_row, 2> cruise_retarget{{
        {0x1.5p+2, 0x1.1p+2, 0x1p+0, 0x0p+0},
        {0x1.98p+3, 0x1.75cp+3, 0x1.7p-1, -0x1.8p-1},
    }};
    auto planner = make_planner<double>({.v_max = 2.0, .a_max = 2.0, .j_max = 1.0});
    planner.reset(0.0);
    REQUIRE(planner.update(2.0).has_value());
    REQUIRE(planner.diagnostics().planned_duration == 4.0);
    for (auto const& row : degenerate) {
        auto const point = planner.sample(row[0]);
        CAPTURE(row[0]);
        CHECK(point.position[0] == row[1]);
        CHECK(point.velocity[0] == row[2]);
        CHECK(point.acceleration[0] == row[3]);
    }
    auto retargeting = make_planner<double>({.v_max = 1.0, .a_max = 1.0, .j_max = 1.0});
    retargeting.reset(0.0);
    REQUIRE(retargeting.update(10.0).has_value());
    REQUIRE(retargeting.sample(4.0).velocity[0] == 1.0);
    REQUIRE(retargeting.update(12.0).has_value());
    for (auto const& row : cruise_retarget) {
        auto const point = retargeting.sample(row[0]);
        CAPTURE(row[0]);
        CHECK(point.position[0] == row[1]);
        CHECK(point.velocity[0] == row[2]);
        CHECK(point.acceleration[0] == row[3]);
    }
}

// -- Test 23: the zero-displacement floor scales with the planner's own limits --
//
// The headline consequence of scaling the floor rather than fixing it. The SAME
// absolute displacement is a real move on one axis and nothing at all on
// another, because the two axes resolve different lengths. The scale is the
// planner's own stopping distance from full speed, v_max^2 / (2 a_max), which is
// intrinsic to the limits it was given and needs no knowledge of the sample
// period -- which the planner is never told, and which is why an absolute
// constant was never the right answer here.
TEST_CASE("OnlinePlanner3rd: the zero-displacement floor scales with the limits",
          "[trajectory][online_planner_3rd][no_op_floor]")
{
    // One displacement, commanded on both axes.
    constexpr double displacement = 1e-15;

    constexpr int length_rounding_ops = 12 + 1 + 3;

    // The coarse axis resolves 1.25 m of stopping distance; the fine one
    // resolves 5 mm. The commanded displacement straddles their two floors.
    constexpr std::array<std::array<double, 3>, 2> limit_sets{{
        {5.0, 10.0, 50.0},
        {1.0, 100.0, 1000.0},
    }};

    int coarse_axes = 0;
    int fine_axes = 0;
    for (auto const& limits : limit_sets) {
        double const v_max = limits[0];
        double const a_max = limits[1];
        double const j_max = limits[2];

        double const length_scale = v_max * v_max / (2.0 * a_max);
        double const length_floor = static_cast<double>(length_rounding_ops)
                                    * std::numeric_limits<double>::epsilon() * length_scale;

        auto planner = make_planner<double>({.v_max = v_max, .a_max = a_max, .j_max = j_max});
        planner.reset(0.0);
        REQUIRE(planner.update(100.0 * v_max).has_value());

        // A state carrying a speed comfortably above the planner's own speed
        // floor and a position still at the origin, so the commanded
        // displacement is not swallowed by the rounding of the position itself.
        // A time s into the opening jerk ramp leaves a speed j*s^2/2.
        constexpr int speed_rounding_ops = 7;
        double const speed_wanted = 20.0 * static_cast<double>(speed_rounding_ops)
                                    * std::numeric_limits<double>::epsilon() * v_max;
        double const s = std::sqrt(2.0 * speed_wanted / j_max);
        auto const probe = planner.sample(s);

        REQUIRE(planner.update(probe.position[0] + displacement).has_value());
        auto const& diag = planner.diagnostics();

        CAPTURE(v_max, a_max, std::log10(length_scale), std::log10(length_floor),
                std::log10(displacement), std::log10(std::abs(probe.position[0])),
                std::log10(std::abs(probe.velocity[0])));

        if (displacement < length_floor) {
            // Below the floor the commanded displacement is zero to the precision
            // available, so the planner treats the command as a direction it
            // cannot resolve and brakes to rest instead of carrying the velocity.
            REQUIRE(diag.disposition
                    == ctrlpp::online_planner_disposition::braked_and_replanned);
            REQUIRE(diag.substitution_reason
                    == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
            ++coarse_axes;
        } else {
            // Above it the same displacement is a move, and is planned as one.
            REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::commanded_profile);
            REQUIRE(diag.substitution_reason
                    == ctrlpp::online_planner_substitution_reason::none);
            ++fine_axes;
        }
    }

    // One of each, asserted exactly. A displacement that fell on the same side of
    // both floors would leave the scaling untested and this is what says so.
    REQUIRE(coarse_axes == 1);
    REQUIRE(fine_axes == 1);
}

// -- Test 24: the overshoot verdict is relative, so it is scale-invariant -------
//
// The overshoot test is the carry-velocity shape's own minimum-displacement
// test: the remaining distance against what the transition from the current
// velocity to rest sweeps, both lengths the planner computes from its own state
// and limits, with no slack and no external distance. Its verdict is therefore a
// property of the RATIO of the two lengths and nothing else, and the same
// relative shortfall must be called an overshoot on an axis whose stopping
// distance is meters and on one whose stopping distance is picometers.
//
// The lower rungs of the ladder are where an absolute slack could not follow: at
// the smallest limit set the shortfall the axis is asked to resolve is some
// three decades below the absolute slack this comparison used to carry, so that
// comparison would have reported no overshoot there.
TEST_CASE("OnlinePlanner3rd: the overshoot verdict is relative and scale-invariant",
          "[trajectory][online_planner_3rd][overshoot]")
{
    constexpr std::array<std::array<double, 3>, 4> limit_sets{{
        {5.0, 10.0, 50.0},
        {1e-2, 5.0, 50.0},
        {1e-4, 5.0, 50.0},
        {1e-5, 5.0, 50.0},
    }};

    // A relative shortfall many decades above the rounding of either length, so
    // what is being tested is the relativity and not the boundary.
    constexpr double relative_shortfall = 1e-6;

    // Cruise leaves the acceleration within one rounding of zero, half a unit in
    // the last place of the ramp product j*T <= a_max, as the no-op case
    // derives; a fused multiply-add is what leaves it nonzero.
    constexpr double cruise_acceleration_rounding_ops = 1.0;
    constexpr double unit_roundoff = 0.5 * std::numeric_limits<double>::epsilon();

    int rungs = 0;
    for (auto const& limits : limit_sets) {
        double const v_max = limits[0];
        double const a_max = limits[1];
        double const j_max = limits[2];

        auto planner = make_planner<double>({.v_max = v_max, .a_max = a_max, .j_max = j_max});
        planner.reset(0.0);

        // The sampling step spans four times the whole acceleration stretch
        // a_max / j_max + v_max / a_max over the drive, so the axis is genuinely
        // CRUISING at every rung rather than still climbing the opening ramp.
        // That matters: a state still on the ramp carries an acceleration, which
        // sends the command through the acceleration-nulling phase first and
        // moves the start of the plan, and the verdict would then be about that
        // phase rather than about the overshoot comparison this case exists to
        // exercise.
        double const step = 4.0 * (a_max / j_max + v_max / a_max) / 100.0;
        auto const cruising = cruise_up(planner, 1000.0 * v_max, step, 400);
        REQUIRE_THAT(cruising.v, WithinAbs(v_max, v_max * 1e-9));
        REQUIRE(std::abs(planner.sample(cruising.t).acceleration[0])
                <= cruise_acceleration_rounding_ops * unit_roundoff * a_max);

        double const stop_dist = jerk_limited_stop_distance(cruising.v, a_max, j_max);
        REQUIRE(stop_dist > 0.0);

        double const target = cruising.q + stop_dist * (1.0 - relative_shortfall);
        REQUIRE(planner.update(target).has_value());
        auto const& diag = planner.diagnostics();

        CAPTURE(v_max, std::log10(cruising.v), std::log10(stop_dist),
                std::log10(stop_dist * relative_shortfall));

        REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::braked_and_replanned);
        REQUIRE(diag.substitution_reason
                == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
        ++rungs;
    }
    REQUIRE(rungs == 4);
}

// -- What sample() carries forward, pinned to recorded bits --------------------
//
// sample() returns a point and also writes the state the next update() replans
// from and is_settled() reads. Both halves are pinned here, so a change to the
// sample body that alters either one fails rather than moving silently.
//
// The limits, the command and every sample time are chosen so each phase
// duration and each sampled state is a short dyadic rational. Every product the
// evaluation forms is then exact, which makes the recorded bits independent of
// whether a compiler fuses a multiply-add, and lets one table serve both float
// and double.

namespace
{

struct recorded_state
{
    double q;
    double v;
    double a;
};

constexpr std::array<recorded_state, 16> recorded_samples{{
        {0x0p+0, 0x0p+0, 0x0p+0},
        {0x1.bp-5, 0x1.bp-2, 0x1.2p+1},
        {0x1.ap-2, 0x1.8p+0, 0x1.8p+1},
        {0x1.2d8p+0, 0x1.4ap+1, 0x1.2p+1},
        {0x1.2p+1, 0x1.8p+1, 0x0p+0},
        {0x1.bp+1, 0x1.8p+1, 0x0p+0},
        {0x1.2p+2, 0x1.8p+1, 0x0p+0},
        {0x1.68p+2, 0x1.8p+1, 0x0p+0},
        {0x1.bp+2, 0x1.8p+1, 0x0p+0},
        {0x1.f4ap+2, 0x1.4ap+1, -0x1.2p+1},
        {0x1.13p+3, 0x1.8p+0, -0x1.8p+1},
        {0x1.1e5p+3, 0x1.bp-2, -0x1.2p+1},
        {0x1.2p+3, 0x0p+0, 0x0p+0},
        {0x1.2p+3, 0x0p+0, 0x0p+0},
        {0x1.2p+3, 0x0p+0, 0x0p+0},
        {0x1.2p+3, 0x0p+0, 0x0p+0}
}};

// Sampled a quarter second after a replan to the origin issued from each sample
// above, so each row also depends on the time sample() recorded.
constexpr std::array<recorded_state, 16> recorded_replans{{
        {0x0p+0, 0x0p+0, 0x0p+0},
        {0x1.b4p-3, 0x1.98p-1, 0x1.8p-1},
        {0x1.b8p-1, 0x1.08p+1, 0x1.8p+0},
        {0x1.e08p+0, 0x1.7ap+1, 0x1.8p-1},
        {0x1.7ep+1, 0x1.68p+1, -0x1.8p+0},
        {0x1.07p+2, 0x1.68p+1, -0x1.8p+0},
        {0x1.4fp+2, 0x1.68p+1, -0x1.8p+0},
        {0x1.97p+2, 0x1.68p+1, -0x1.8p+0},
        {0x1.dfp+2, 0x1.68p+1, -0x1.8p+0},
        {0x1.0d3p+3, 0x1.1ap+1, -0x1.8p-1},
        {0x1.1c8p+3, 0x1.ep-1, -0x1.8p+0},
        {0x1.1ffp+3, 0x1.8p-5, -0x1.8p-1},
        {0x1.1f8p+3, -0x1.8p-3, -0x1.8p+0},
        {0x1.1f8p+3, -0x1.8p-3, -0x1.8p+0},
        {0x1.1f8p+3, -0x1.8p-3, -0x1.8p+0},
        {0x1.1f8p+3, -0x1.8p-3, -0x1.8p+0}
}};

constexpr std::int32_t recorded_first_settled_sample = 12;

template <typename Scalar>
void require_same_bits(Scalar actual, double recorded)
{
    using bits = std::conditional_t<sizeof(Scalar) == sizeof(std::uint64_t), std::uint64_t, std::uint32_t>;
    REQUIRE(std::bit_cast<bits>(actual) == std::bit_cast<bits>(static_cast<Scalar>(recorded)));
}

template <typename Scalar>
void require_recorded(ctrlpp::trajectory_point<Scalar, 1> const& point, recorded_state const& recorded)
{
    require_same_bits(point.position(0), recorded.q);
    require_same_bits(point.velocity(0), recorded.v);
    require_same_bits(point.acceleration(0), recorded.a);
}

}

TEMPLATE_TEST_CASE("OnlinePlanner3rd: sample returns and carries exactly the recorded state",
                   "[trajectory][online_planner_3rd][sample_contract]", float, double)
{
    using Scalar = TestType;
    auto planner = make_planner<Scalar>({.v_max = Scalar{3}, .a_max = Scalar{3}, .j_max = Scalar{6}});
    REQUIRE(planner.update(Scalar{9}).has_value());
    REQUIRE(planner.diagnostics().planned_duration == Scalar{4.5});

    for (std::int32_t k = 0; k < static_cast<std::int32_t>(recorded_samples.size()); ++k) {
        auto const index = static_cast<std::size_t>(k);
        auto const t = static_cast<Scalar>(k) * Scalar{3} / Scalar{8};
        auto const point = planner.sample(t);
        CAPTURE(k, t);
        require_recorded(point, recorded_samples[index]);
        REQUIRE(planner.is_settled() == (k >= recorded_first_settled_sample));

        // A replan snapshots the carried state, and sampling it at the time
        // sample() recorded integrates no phase, so it hands that state back.
        auto replanned = planner;
        REQUIRE(replanned.update(Scalar{0}).has_value());
        REQUIRE(replanned.diagnostics().initial_velocity == point.velocity(0));
        auto const resumed = replanned.sample(t);
        require_recorded(resumed, recorded_samples[index]);
        require_recorded(replanned.sample(t + Scalar{0.25}), recorded_replans[index]);
    }
}

namespace
{

template <typename Scalar>
struct rest_to_rest_command
{
    Scalar h, v_max, a_max, j_max;
};

// Rest-to-rest commands whose cruise-free peak, as solved before it was formed from the limit side, rounded above
// v_max by one to six ulps: the quadratic root where a_max is reached and the doubly degenerate peak.
template <typename Scalar>
auto rest_to_rest_peak_commands() -> std::vector<rest_to_rest_command<Scalar>>
{
    if constexpr (std::is_same_v<Scalar, double>) {
        return {{0x1.cp+1, 5, 10, 50},
                {0x1.61a58d1b06a78p-1, 0.7, 1.3, 2.9},
                {0x1.5600d03ed0857p+20, 0x1.67f9e8721c5a7p+9, 0x1.7ae695a7dde87p-2, 0x1.4acebc3ac628dp+3},
                {0x1.974b7e1e7faf1p-9, 0x1.54ea2b0c5c43p-8, 0x1.1ff1d3ec394f9p-7, 0x1.a26b26a35ab2ap+0},
                {0x1.ffffffffffffep+0, 1, 1, 1},
                {0x1.24d037bce8b6cp+5, 0x1.8851333f8c53cp+0, 0x1.4638aca7a722bp+1, 0x1.6020b4e0da91ap-7},
                {0x1.f5a47028a379ap+5, 0x1.4b094ced91319p+3, 0x1.97e65b96537edp+8, 0x1.2050df7387627p+0},
                {0x1.5707cbf3ba2a2p-5, 0x1.369a4d44f024p-7, 0x1.193d9294d762fp-6, 0x1.fd4ed86642263p-10},
                {0x1.7a8f254b9fa5bp+1, 0x1.c6a6fd2303d4fp+1, 0x1.110535fa05b0dp+5, 0x1.47e640cce757fp+4}};
    } else {
        return {{0x1.32c914p+8f, 0x1.536974p+2f, 0x1.201be2p-3f, 0x1.c98c8p-8f},
                {0x1.806a76p+10f, 0x1.355ae8p+2f, 0x1.f479cap-7f, 0x1.31d874p-7f},
                {0x1.32c712p+1f, 0x1.2546d4p+1f, 0x1.234036p+1f, 0x1.d1cd4cp+5f},
                {0x1.ad4736p+15f, 0x1.a23aaep+9f, 0x1.89ba2ep+8f, 0x1.8cfa26p-1f},
                {0x1.9f3f24p+7f, 0x1.37faaep+2f, 0x1.4148fep+6f, 0x1.60345p-7f},
                {0x1.8bca8ap-17f, 0x1.468358p-10f, 0x1.4352a2p+5f, 0x1.bc6d0ap+5f},
                {0x1.62118p-5f, 0x1.f59a3p-7f, 0x1.634e36p-5f, 0x1.f75adep-8f},
                {0x1.1ec1d4p-16f, 0x1.9fd178p-8f, 0x1.2d7ba6p+4f, 0x1.b52bf8p+11f}};
    }
}

}

// The peak is pinned at the end of the ramp, which a symmetric rest-to-rest plan reaches at half its duration.
TEMPLATE_TEST_CASE("OnlinePlanner3rd: a rest-to-rest peak next to the cruise boundary is pinned at or below the limit",
                   "[trajectory][online_planner_3rd]", float, double)
{
    using Scalar = TestType;
    for (auto const& c : rest_to_rest_peak_commands<Scalar>()) {
        for (auto const sigma : {Scalar{1}, Scalar{-1}}) {
            CAPTURE(c.h, c.v_max, c.a_max, c.j_max, sigma);
            auto planner = make_planner<Scalar>({.v_max = c.v_max, .a_max = c.a_max, .j_max = c.j_max});
            REQUIRE(planner.update(sigma * c.h).has_value());
            auto const T = planner.diagnostics().planned_duration;
            auto t = T / Scalar{2};
            for (std::int32_t i = 0; i < 64; ++i) {
                t = std::nextafter(t, Scalar{0});
            }
            std::int32_t above = 0;
            for (std::int32_t i = 0; i <= 128; ++i, t = std::nextafter(t, T)) {
                above += (std::abs(planner.sample(t).velocity[0]) > c.v_max) ? 1 : 0;
            }
            for (std::int32_t i = 0; i <= 1024; ++i) {
                above += (std::abs(planner.sample(T * static_cast<Scalar>(i) / Scalar{1024}).velocity[0]) > c.v_max) ? 1 : 0;
            }
            CHECK(above == 0);
        }
    }
}
