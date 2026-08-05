#include "ctrlpp/trajectory/online_planner_3rd.h"

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
TEST_CASE("OnlinePlanner3rd: step response settles", "[traj][online_planner_3rd]")
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
          "[traj][online_planner_3rd][negative]")
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
          "[traj][online_planner_3rd][precision]")
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
TEST_CASE("OnlinePlanner3rd: velocity constraint", "[traj][online_planner_3rd]")
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
TEST_CASE("OnlinePlanner3rd: acceleration constraint", "[traj][online_planner_3rd]")
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
TEST_CASE("OnlinePlanner3rd: jerk constraint", "[traj][online_planner_3rd]")
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
TEST_CASE("OnlinePlanner3rd: mid-motion target change", "[traj][online_planner_3rd]")
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
TEST_CASE("OnlinePlanner3rd: negative displacement", "[traj][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    REQUIRE(planner.update(-5.0).has_value());

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(-5.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));
}

// -- Test 7: is_settled returns correct state ----------------------------------
TEST_CASE("OnlinePlanner3rd: is_settled transitions", "[traj][online_planner_3rd]")
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
TEST_CASE("OnlinePlanner3rd: reset", "[traj][online_planner_3rd]")
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
TEST_CASE("OnlinePlanner3rd: smooth acceleration profile", "[traj][online_planner_3rd]")
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
TEST_CASE("OnlinePlanner3rd: float type", "[traj][online_planner_3rd]")
{
    auto planner = make_planner<float>({.v_max = 5.0f, .a_max = 10.0f, .j_max = 50.0f});

    REQUIRE(planner.update(10.0f).has_value());
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
          "[traj][online_planner_3rd][negative]")
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

// -- Test 13: a substituted carry-velocity shape is REPORTED --------------------
//
// The planner selects the brake-then-replan fallback from two tests of its own
// (direction, and stopping distance against remaining distance), and then from a
// third condition those two do not predict: the double-S that carries the
// current velocity through the profile has a domain, and a commanded state can
// pass both tests and still fall outside it.
//
// The resulting motion is limit-respecting and reaches the target, which is
// exactly why the substitution went unreported for so long: nothing about the
// motion distinguishes it. So the assertions below are on the REPORT. The motion
// is checked afterwards, and would pass either way.
TEST_CASE("OnlinePlanner3rd: carry-velocity shape substitution is reported",
          "[traj][online_planner_3rd][diagnostics]")
{
    double constexpr v_max = 5.0;
    double constexpr a_max = 10.0;
    double constexpr j_max = 50.0;

    auto planner = make_planner<double>({.v_max = v_max, .a_max = a_max, .j_max = j_max});

    auto const cruising = cruise_up(planner, 100.0, 0.01, 200);
    REQUIRE_THAT(cruising.v, WithinAbs(v_max, 1e-9));

    double const stop_dist = jerk_limited_stop_distance(cruising.v, a_max, j_max);

    // Where the commanded target has to sit, read off the two branch conditions
    // rather than searched for:
    //
    //  * the planner calls it an overshoot when the stopping distance exceeds
    //    the remaining distance by more than its own comparison slack, and that
    //    slack is RELATIVE: the remaining distance is scaled by one plus a
    //    counted number of roundings, so a target short of the stopping point
    //    by less than the stopping distance times that slack passes the test;
    //  * the carry-velocity double-S exists only when the remaining distance is
    //    at least what the transition from the current velocity to rest already
    //    sweeps, which is the stopping distance outright, so a target short of
    //    the stopping point by ANY amount has no such shape.
    //
    // Between the two lies a band of commanded targets one slack wide that the
    // planner accepts and the shape cannot serve. Command its midpoint.
    //
    // This test is coupled to the planner's comparison slack by construction:
    // the slack is what opens the band, and a redesign of that constant changes
    // where the band is. That coupling is the point, not an accident of the
    // test. The count below is the planner's own: thirty-three roundings along
    // the stopping-distance chain that reaches the acceleration limit, plus
    // thirteen along the chain that forms the remaining distance.
    int constexpr planner_overshoot_rounding_ops = 33 + 13;
    double const band = stop_dist * static_cast<double>(planner_overshoot_rounding_ops)
                        * std::numeric_limits<double>::epsilon();
    double const shortfall = band / 2.0;
    double const target = cruising.q + stop_dist - shortfall;

    // The band is now a few tens of units in the last place of the stopping
    // distance rather than an absolute picometre, so whether the command landed
    // inside it is no longer self-evident from the arithmetic above: forming the
    // target rounds twice at the scale of the stopping POINT, which is larger
    // than the stopping distance. Assert the premise rather than assume it, so a
    // miss reports itself as a miss instead of as a surprising branch.
    double const realized_shortfall = stop_dist - (target - cruising.q);
    CAPTURE(cruising.q, cruising.v, stop_dist, band, shortfall, realized_shortfall, target);
    REQUIRE(realized_shortfall > 0.0);
    REQUIRE(realized_shortfall < band);

    REQUIRE(planner.update(target).has_value());

    // The report, first.
    auto const& diag = planner.diagnostics();
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(diag.substitution_reason
            == ctrlpp::online_planner_substitution_reason::carry_velocity_shape_unavailable);
    REQUIRE(diag.brake_duration > 0.0);
    REQUIRE(diag.commanded_target == target);
    REQUIRE(diag.initial_velocity == cruising.v);

    // Commanded against planned: the replan does not start where the command
    // did, it starts where the braking ended.
    REQUIRE(diag.replan_start_position != cruising.q);
    REQUIRE_THAT(diag.replan_start_position, WithinAbs(cruising.q + stop_dist, 1e-9));

    // Braking lands one shortfall past the commanded target. That displacement
    // is nonzero, so the bounded return profile has a positive duration even
    // though it lies below the planner's comparison slack.
    REQUIRE(diag.planned_duration > diag.brake_duration);
    REQUIRE(std::isfinite(diag.planned_duration));

    // The motion, second. It respects every limit and ends at the target, as it
    // did before the substitution was reported at all.
    double constexpr tol = 1e-6;
    double const t_end = cruising.t + diag.planned_duration;
    for (double t = cruising.t; t < t_end; t += 0.001) {
        auto const pt = planner.sample(t);
        REQUIRE(std::abs(pt.velocity[0]) <= v_max + tol);
        REQUIRE(std::abs(pt.acceleration[0]) <= a_max + tol);
        // The motion runs from the commanded start to the stopping point, never
        // turning back: the interval it stays inside is the one the report names.
        REQUIRE(pt.position[0] >= cruising.q - tol);
        REQUIRE(pt.position[0] <= diag.replan_start_position + tol);
    }

    auto const settled = planner.sample(t_end + 1.0);
    REQUIRE_THAT(settled.position[0], WithinAbs(target, 1e-6));
    REQUIRE_THAT(settled.velocity[0], WithinAbs(0.0, 1e-6));
}

// -- Test 14: the reversal substitution is reported, and is DISTINGUISHABLE -----
//
// Same disposition, different cause. A caller that only learned "something was
// substituted" could not choose between axes or explain the extra time; the
// reason enumerator is what makes the two causes separable.
TEST_CASE("OnlinePlanner3rd: reversal substitution reports a distinct reason",
          "[traj][online_planner_3rd][diagnostics]")
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
          "[traj][online_planner_3rd][diagnostics]")
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

// -- Test 16: a move commanded from inside the settle tolerance reports settled -
TEST_CASE("OnlinePlanner3rd: a settled command reports a zero-duration plan",
          "[traj][online_planner_3rd][diagnostics]")
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
          "[traj][online_planner_3rd][diagnostics]")
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
          "[traj][online_planner_3rd][settle_tolerance]")
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
          "[traj][online_planner_3rd][settle_tolerance]")
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
          "[traj][online_planner_3rd][settle_tolerance]")
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

// -- Test 21: the numerical-no-op floors are derived, and they scale ------------
//
// A separate question from the settle policy above, with a separate owner. The
// policy says when the application considers the axis arrived; this says when
// the arithmetic can no longer tell the commanded state from the one already in
// force. The floors are one unit in the last place at the scale of the limit
// that bounds each quantity, so they move with the limits rather than sitting at
// an absolute distance borrowed from an axis nobody named.
//
// The acceleration clause is the one that decides here, and that is structural
// rather than incidental: every state this planner reaches on the way to rest
// leaves the acceleration above its own floor before the velocity reaches
// velocity's. Sampling a jerk ramp a time s from rest gives an acceleration
// j*s and a speed j*s^2/2, so the two margins stand in the ratio s*a_max/(2
// v_max), which is below one for every s short of the whole acceleration
// stretch. The velocity clause is therefore exercised in the other direction,
// below, where a state carries a real speed and exactly zero acceleration.
TEST_CASE("OnlinePlanner3rd: the numerical-no-op floors scale with the limits",
          "[traj][online_planner_3rd][no_op_floor]")
{
    // Three limit sets spanning five decades of acceleration limit, with the
    // ratios held so the sweep is a sweep of scale and not of branch structure.
    constexpr std::array<std::array<double, 3>, 3> limit_sets{{
        {5.0, 10.0, 50.0},
        {0.25, 0.5, 4.0},
        {40.0, 200.0, 2000.0},
    }};

    int straddles = 0;
    for (auto const& limits : limit_sets) {
        double const v_max = limits[0];
        double const a_max = limits[1];
        double const j_max = limits[2];

        // One operation forms the compared quantity -- the comparison itself --
        // at the scale of the limit that bounds it.
        constexpr int command_acceleration_rounding_ops = 1;
        double const acceleration_floor = static_cast<double>(command_acceleration_rounding_ops)
                                          * std::numeric_limits<double>::epsilon() * a_max;

        for (double const half_or_double : {0.5, 2.0}) {
            auto planner = make_planner<double>(
                {.v_max = v_max, .a_max = a_max, .j_max = j_max});
            planner.reset(0.0);
            REQUIRE(planner.update(100.0 * v_max).has_value());

            // A time s into the opening jerk ramp leaves an acceleration j*s and
            // a speed j*s^2/2, so s is chosen to place the acceleration at half
            // the floor and at twice it.
            double const s = half_or_double * acceleration_floor / j_max;
            auto const probe = planner.sample(s);

            // Commanding the position just sampled makes the position clause an
            // exact equality, so the two magnitude clauses are what decide.
            REQUIRE(planner.update(probe.position[0]).has_value());

            CAPTURE(v_max, a_max, j_max, half_or_double,
                    std::log10(std::abs(probe.acceleration[0])),
                    std::log10(acceleration_floor), std::log10(std::abs(probe.velocity[0])));

            // The speed left behind is quadratic in s where the acceleration is
            // linear, so it sits decades below its own floor either way: the
            // acceleration clause is what flips.
            REQUIRE(std::abs(probe.velocity[0])
                    < std::numeric_limits<double>::epsilon() * v_max);

            if (half_or_double < 1.0) {
                REQUIRE(std::abs(probe.acceleration[0]) < acceleration_floor);
                REQUIRE(planner.diagnostics().disposition
                        == ctrlpp::online_planner_disposition::settled);
                REQUIRE(planner.diagnostics().planned_duration == 0.0);
            } else {
                REQUIRE(std::abs(probe.acceleration[0]) > acceleration_floor);
                REQUIRE(planner.diagnostics().disposition
                        != ctrlpp::online_planner_disposition::settled);
            }
        }
        ++straddles;
    }

    // An exact count, not a bound the loop cannot fail.
    REQUIRE(straddles == 3);

    // The velocity clause, from the other side. At cruise the acceleration is
    // EXACTLY zero -- the two jerk ramps that built the cruise added and removed
    // the same value -- so the acceleration clause cannot be what keeps the
    // short-circuit shut, and the speed is what the planner is left deciding on.
    int cruises = 0;
    for (auto const& limits : limit_sets) {
        auto planner = make_planner<double>(
            {.v_max = limits[0], .a_max = limits[1], .j_max = limits[2]});
        planner.reset(0.0);
        auto const cruising = cruise_up(planner, 1000.0 * limits[0], 0.01, 400);
        CAPTURE(limits[0], cruising.q, cruising.v);

        REQUIRE(planner.sample(cruising.t).acceleration[0] == 0.0);
        REQUIRE(std::abs(cruising.v)
                > std::numeric_limits<double>::epsilon() * limits[0]);

        REQUIRE(planner.update(cruising.q).has_value());
        REQUIRE(planner.diagnostics().disposition
                != ctrlpp::online_planner_disposition::settled);
        ++cruises;
    }
    REQUIRE(cruises == 3);
}

// -- Test 22: the acceleration-nulling phase is emitted when it moves the axis --
//
// The phase that brings a starting acceleration to zero used to be gated on the
// acceleration exceeding an absolute number, which is a comparison with no scale
// to be against. It is now gated on whether the phase changes anything: over its
// own duration |a_start| / j_max it changes the velocity by a_start^2 / (2 j_max)
// and moves the axis by at most v_max * |a_start| / j_max, and it is emitted
// unless both are below the resolution of the quantity they are measured
// against. The length side binds first, being linear where the velocity side is
// quadratic, so the boundary sits where
//
//     v_max * a_start / j_max == 5 * eps * v_max^2 / (2 a_max)
//
// -- five being two operations for the contribution and three for the scale.
// Below it the axis is left where it was; above it the phase is emitted and the
// acceleration is nulled, which is directly observable by sampling at the
// phase's own duration.
TEST_CASE("OnlinePlanner3rd: the acceleration-nulling phase is emitted when it moves the axis",
          "[traj][online_planner_3rd][no_op_floor]")
{
    constexpr std::array<std::array<double, 3>, 3> limit_sets{{
        {5.0, 10.0, 50.0},
        {0.25, 0.5, 4.0},
        {40.0, 200.0, 2000.0},
    }};

    constexpr int nulling_length_rounding_ops = 2;
    constexpr int nulling_scale_rounding_ops = 3;
    constexpr int nulling_displacement_rounding_ops = nulling_length_rounding_ops
                                                      + nulling_scale_rounding_ops;

    int straddles = 0;
    for (auto const& limits : limit_sets) {
        double const v_max = limits[0];
        double const a_max = limits[1];
        double const j_max = limits[2];

        // The starting acceleration at which the phase's displacement bound
        // equals the resolution of the stopping distance from full speed.
        double const a_boundary = static_cast<double>(nulling_displacement_rounding_ops)
                                  * std::numeric_limits<double>::epsilon() * v_max * j_max
                                  / (2.0 * a_max);

        for (double const below_or_above : {0.5, 4.0}) {
            auto planner = make_planner<double>(
                {.v_max = v_max, .a_max = a_max, .j_max = j_max});
            planner.reset(0.0);
            REQUIRE(planner.update(100.0 * v_max).has_value());

            double const s = below_or_above * a_boundary / j_max;
            auto const probe = planner.sample(s);
            double const a_start = probe.acceleration[0];
            REQUIRE(a_start > 0.0);

            // Replan to a far target from that state, then sample at exactly the
            // duration the nulling phase would occupy.
            REQUIRE(planner.update(100.0 * v_max).has_value());
            double const nulling_duration = a_start / j_max;
            auto const after = planner.sample(s + nulling_duration);

            CAPTURE(v_max, a_max, j_max, below_or_above, std::log10(a_start),
                    std::log10(a_boundary), std::log10(std::abs(after.acceleration[0]) + 1e-300));

            if (below_or_above < 1.0) {
                REQUIRE(a_start < a_boundary);
                // No phase was emitted, so the acceleration is still climbing the
                // profile's own opening ramp rather than having been nulled.
                REQUIRE(std::abs(after.acceleration[0]) >= a_start);
            } else {
                REQUIRE(a_start > a_boundary);
                // The phase was emitted and did its job: the acceleration is zero
                // at its end, to the resolution of the acceleration limit.
                REQUIRE(std::abs(after.acceleration[0])
                        < std::numeric_limits<double>::epsilon() * a_max);
            }
        }
        ++straddles;
    }
    REQUIRE(straddles == 3);
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
          "[traj][online_planner_3rd][no_op_floor]")
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
// The overshoot test compares two lengths the planner has already computed, with
// a slack that is a counted multiple of epsilon rather than a distance. Its
// verdict is therefore a property of the RATIO of the two lengths and nothing
// else, and the same relative shortfall must be called an overshoot on an axis
// whose stopping distance is metres and on one whose stopping distance is
// picometres.
//
// The lower rungs of the ladder are where an absolute slack could not follow: at
// the smallest limit set the shortfall the axis is asked to resolve is some
// three decades below the absolute slack this comparison used to carry, so that
// comparison would have reported no overshoot there.
TEST_CASE("OnlinePlanner3rd: the overshoot verdict is relative and scale-invariant",
          "[traj][online_planner_3rd][overshoot]")
{
    constexpr std::array<std::array<double, 3>, 4> limit_sets{{
        {5.0, 10.0, 50.0},
        {1e-2, 5.0, 50.0},
        {1e-4, 5.0, 50.0},
        {1e-5, 5.0, 50.0},
    }};

    // A relative shortfall many decades above the counted slack, so what is being
    // tested is the relativity and not the width.
    constexpr double relative_shortfall = 1e-6;

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
        REQUIRE(planner.sample(cruising.t).acceleration[0] == 0.0);

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
