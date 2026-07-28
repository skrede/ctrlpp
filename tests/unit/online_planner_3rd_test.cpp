#include "ctrlpp/trajectory/online_planner_3rd.h"

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
    planner.update(far_target);

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

    planner.update(-5.0);

    auto const pt = planner.sample(100.0);
    REQUIRE_THAT(pt.position[0], WithinAbs(-5.0, 1e-6));
    REQUIRE_THAT(pt.velocity[0], WithinAbs(0.0, 1e-6));
}

// -- Test 7: is_settled returns correct state ----------------------------------
TEST_CASE("OnlinePlanner3rd: is_settled transitions", "[traj][online_planner_3rd]")
{
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

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
    auto planner = make_planner<double>({.v_max = 5.0, .a_max = 10.0, .j_max = 50.0});

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
    auto planner = make_planner<float>({.v_max = 5.0f, .a_max = 10.0f, .j_max = 50.0f});

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
    auto planner = make_planner<double>({.v_max = v_max, .a_max = 10.0, .j_max = 50.0});

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
    //    the remaining distance by more than its own comparison slack, so a
    //    target short of the stopping point by LESS than one slack passes that
    //    test;
    //  * the carry-velocity double-S exists only when the remaining distance is
    //    at least what the transition from the current velocity to rest already
    //    sweeps, which is the stopping distance outright, so a target short of
    //    the stopping point by ANY amount has no such shape.
    //
    // Between the two lies a band of commanded targets one slack wide that the
    // planner accepts and the shape cannot serve. Command its midpoint: half a
    // slack is some five hundred times the rounding of a distance this size, so
    // the band is entered by construction rather than by luck.
    //
    // This test is coupled to the planner's comparison slack by construction:
    // the slack is what opens the band, and a redesign of that constant changes
    // where the band is. That coupling is the point, not an accident of the
    // test.
    double constexpr planner_comparison_slack = 1e-12;
    double const shortfall = planner_comparison_slack / 2.0;
    double const target = cruising.q + stop_dist - shortfall;

    CAPTURE(cruising.q, cruising.v, stop_dist, shortfall, target);

    planner.update(target);

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

    // The plan is at least the braking it contains. Equality is this band's own
    // signature and not a slack assertion: braking from the cruise velocity
    // lands one shortfall PAST the commanded target, so the replan that follows
    // has nothing left to travel and contributes no duration.
    REQUIRE(diag.planned_duration >= diag.brake_duration);
    REQUIRE_THAT(diag.planned_duration, WithinAbs(diag.brake_duration, 2.0 * shortfall));

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

    planner.update(target);

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

    planner.update(10.0);

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
    planner.update(200.0);
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

    planner.update(3.0);

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
    planner.update(cruising.q - 10.0);
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::braked_and_replanned);

    planner.reset(0.0);

    auto const& diag = planner.diagnostics();
    REQUIRE(diag.disposition == ctrlpp::online_planner_disposition::settled);
    REQUIRE(diag.substitution_reason == ctrlpp::online_planner_substitution_reason::none);
    REQUIRE(diag.brake_duration == 0.0);
    REQUIRE(diag.planned_duration == 0.0);
}
