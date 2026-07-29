#include "hardening_helpers.h"

#include "ctrlpp/trajectory/cubic_spline.h"
#include "ctrlpp/trajectory/synchronize.h"
#include "ctrlpp/trajectory/smoothing_spline.h"
#include "ctrlpp/trajectory/bspline_trajectory.h"
#include "ctrlpp/trajectory/online_planner_2nd.h"
#include "ctrlpp/trajectory/online_planner_3rd.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <catch2/catch_test_macros.hpp>

#include <span>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <utility>
#include <algorithm>

// ── What this file's oracles decide ───────────────────────────────────────────
//
// A rescaled profile is judged by integrating the velocity it reports against the
// displacement it was commanded, never by sampling its position near an endpoint.
// The final segment of both profiles is written as an offset backwards from the
// commanded displacement, so it returns the target position by construction: a
// profile whose velocity integrates to something else entirely still lands its
// endpoint exactly. Position is therefore asserted only through quadrature, and
// that quadrature is laid out inside each phase segment. One uniform grid over the
// whole duration is rejected: its panels straddle the phase kinks, where the
// integrand's slope jumps, and the truncation error manufactured there is orders
// of magnitude above the rounding floor these budgets describe.
//
// The out-of-domain cases assert a typed rejection naming the enumerator the
// domain reasoning predicts, and the six planner limit cases sweep the whole
// out-of-domain set -- zero, negative, NaN and infinite -- rather than one
// representative of it. Where a rejection must also leave the object alone, the
// case says so: no value, the specific enumerator, the duration unchanged, the
// phase durations unchanged.
//
// The cases that exist to reach a particular branch assert what that branch
// produces, not that a number came back. Three kinds of statement carry them,
// and which one applies is a fact about the evaluation path rather than a choice:
//
//  * Exact equality where the evaluation ASSIGNS or STORES the quantity. Past the
//    end of a planned profile the planners return the stored target with zero
//    velocity and zero acceleration; at a knot a spline's Horner form reduces to
//    its stored constant term; in cruise the acceleration is assigned zero and
//    the velocity the limit itself; a zero-displacement profile has duration zero
//    and starts where it was told. No epsilon takes part in any of those.
//  * A counted rounding budget where the evaluation COMPUTES the quantity, with
//    the operation chain enumerated in prose beside the constant and the scale
//    named. Never a round number, and never a magnitude fitted to what the code
//    happens to return today.
//  * The kinematic envelope inside every sampling loop -- velocity within the
//    velocity limit, acceleration within the acceleration limit, and for the
//    third-order planner the reported acceleration slewing no faster than the
//    jerk limit. That envelope is what a limit-respecting planner guarantees, and
//    a loop that samples hundreds of points and checks only that they are numbers
//    passes a planner that exceeds its acceleration limit throughout.
//
// Two quantities are deliberately NOT pinned, and the reason is the domain rather
// than convenience. A smoothing spline's departure from the straight line it tends
// to as its smoothing weight grows is first order in the tradeoff parameter with a
// coefficient that depends on the data through the second-difference operator's
// pseudo-inverse; the cases bound that departure against the interpolating
// spline's own departure scaled by the parameter, which is the first-order law
// itself and needs no coefficient. A rescaled profile's realized duration is the
// sum of its realized phase durations and is never assigned the request, so exact
// equality is not its contract either.
namespace
{

/// Panels laid inside each phase segment. Simpson's rule integrates every
/// polynomial up to cubic order exactly and these velocity profiles are at most
/// quadratic inside a segment, so once the panels are aligned the quadrature
/// carries no truncation error and the count buys no accuracy.
constexpr int panels_per_segment = 8;

/// create() is the only construction path on the two velocity profiles and it is
/// fallible, so every profile these cases use is built through one of these
/// helpers, which assert the command was realizable. The negative cases below do
/// not come through here: they assert the specific enumerator.
auto trapezoidal_profile(ctrlpp::trapezoidal_trajectory<double>::config const& cfg)
    -> ctrlpp::trapezoidal_trajectory<double>
{
    auto created = ctrlpp::trapezoidal_trajectory<double>::create(cfg);
    REQUIRE(created.has_value());
    return created.value();
}

auto double_s_profile(ctrlpp::double_s_trajectory<double>::config const& cfg)
    -> ctrlpp::double_s_trajectory<double>
{
    auto created = ctrlpp::double_s_trajectory<double>::create(cfg);
    REQUIRE(created.has_value());
    return created.value();
}

/// The splines and the online planners are fallible-only as well, so the same
/// build-and-assert shape covers them.
template <typename Type>
auto realizable(typename Type::config const& cfg) -> Type
{
    auto created = Type::create(cfg);
    REQUIRE(created.has_value());
    return *std::move(created);
}

/// Chained rounding operations behind one Simpson panel: each of its two fresh
/// velocity samples chains up to five multiply-adds inside evaluate(), the panel
/// itself three multiplies and three adds over its samples, and the running sums
/// one addition each. Every one is worth up to one unit in the last place at the
/// scale of the accumulated absolute area.
constexpr int rounding_ops_per_panel = 2 * 5 + 6 + 2;

/// Chained rounding operations behind a single reported quantity: a sum of phase
/// durations, or one velocity sample compared against a boundary value.
constexpr int rounding_ops_per_sample = 8;

/// Chained rounding operations behind a rescaled profile's realized duration: at
/// most six for the closed form that produces the cruise velocity (a square root
/// and five arithmetic operations), two each for the two ramp durations, three
/// each for the two ramp distances, three for the cruise duration, and two for
/// the final sum. Each is worth up to one unit in the last place at the scale of
/// the requested duration.
constexpr int rounding_ops_per_duration = 6 + 4 + 6 + 3 + 2;

struct quadrature
{
    double integral{};
    double abs_area{};
    int panels{};
    int boundaries{};
    double widest_panel{};
};

/// Kink-aligned composite Simpson integration of the reported velocity over the
/// profile's own reported duration. Segment lengths come from the profile itself,
/// are turned into boundaries by prefix sum, and each segment is integrated
/// separately.
template <typename Profile>
auto integrate_velocity(Profile const& profile) -> quadrature
{
    auto const segments = profile.phase_durations();

    quadrature out{};
    double t_start = 0.0;
    for (auto const& segment : segments) {
        double const length = static_cast<double>(segment);
        if (!(length > 0.0)) {
            continue;
        }
        ++out.boundaries;

        double const dt = length / static_cast<double>(panels_per_segment);
        out.widest_panel = std::max(out.widest_panel, dt);
        for (int p = 0; p < panels_per_segment; ++p) {
            double const a = t_start + static_cast<double>(p) * dt;
            double const b = a + dt;
            double const m = a + 0.5 * dt;
            double const panel = (dt / 6.0)
                * (profile.evaluate(a).velocity(0) + 4.0 * profile.evaluate(m).velocity(0)
                   + profile.evaluate(b).velocity(0));
            out.integral += panel;
            out.abs_area += std::abs(panel);
            ++out.panels;
        }
        t_start += length;
    }
    return out;
}

/// Assert that the profile sweeps the displacement it was commanded.
template <typename Profile>
void require_swept_displacement(Profile const& profile, double h_signed, double a_max)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();

    auto const q = integrate_velocity(profile);
    double const area_scale = std::max(q.abs_area, std::abs(h_signed));
    double const panel_tol =
        static_cast<double>(rounding_ops_per_panel * q.panels) * eps * area_scale;

    // The sample taken exactly on a segment boundary falls into the following
    // segment, where evaluate() recovers its local time by subtracting the sample
    // time from the total duration. That difference is representable only to one
    // unit in the last place of the duration, and the velocity slews at up to the
    // acceleration limit there, so the sample carries a_max * eps * T of velocity
    // uncertainty. Simpson weights an endpoint by a sixth of its panel width.
    double const boundary_tol = static_cast<double>(q.boundaries) * (q.widest_panel / 6.0)
                                * a_max * eps * static_cast<double>(profile.duration());

    CAPTURE(q.integral, h_signed, q.abs_area, q.panels, panel_tol, boundary_tol);
    REQUIRE(std::abs(q.integral - h_signed) <= panel_tol + boundary_tol);
}

/// Assert that the profile arrives at the boundary velocity it was commanded.
///
/// The bound is the acceleration limit times the step taken back from the end of
/// the move: the velocity cannot be further from its final value than the
/// acceleration is allowed to move it over that step. A profile that ended at its
/// cruise velocity instead of the commanded one misses by the full span of its
/// deceleration ramp, which is twice this bound.
template <typename Profile>
void require_terminal_velocity(Profile const& profile, double v1, double a_max)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();

    auto const segments = profile.phase_durations();
    double const delta = static_cast<double>(segments.back()) / 2.0;
    if (!(delta > 0.0)) {
        return;
    }

    double const T = static_cast<double>(profile.duration());
    double const v = profile.evaluate(T - delta).velocity(0);
    double const slew_bound = a_max * delta;
    double const sample_tol = static_cast<double>(rounding_ops_per_sample) * eps
                              * (std::max(std::abs(v1), 1.0) + a_max * T);

    CAPTURE(delta, v, v1, slew_bound, sample_tol);
    REQUIRE(std::abs(v - v1) <= slew_bound + sample_tol);
}

/// Assert that no phase duration came out negative and that they sum to the
/// reported duration.
template <typename Profile>
void require_nonnegative_phases(Profile const& profile)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();

    double sum = 0.0;
    for (auto const& segment : profile.phase_durations()) {
        double const length = static_cast<double>(segment);
        REQUIRE(std::isfinite(length));
        REQUIRE(length >= 0.0);
        sum += length;
    }
    double const T = static_cast<double>(profile.duration());
    CAPTURE(sum, T);
    REQUIRE(std::abs(sum - T) <= static_cast<double>(rounding_ops_per_sample) * eps * T);
}

/// Assert that the realized duration lands on the requested one to within the
/// rounding of the expressions that produced it. The stored duration is the sum of
/// the realized phase durations and is never assigned the request, so exact
/// equality is not the contract.
template <typename Profile>
void require_realized_duration(Profile const& profile, double T_new)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();

    double const T = static_cast<double>(profile.duration());
    double const tol = static_cast<double>(rounding_ops_per_duration) * eps * T_new;
    CAPTURE(T, T_new, tol);
    REQUIRE(std::abs(T - T_new) <= tol);
}

/// Chained rounding operations behind one velocity an online planner reports.
/// Every branch forms it the same way: a product of an acceleration and a time
/// added to a carried velocity, or a product of a deceleration and a time
/// remaining. The time itself is a difference of two sample times and, on the
/// backward-time branches, a difference of two durations; each phase duration the
/// branch selection compares against is a difference and a division. A product, a
/// sum, two differences and two divisions is six, and the two sums that formed the
/// total duration the selection ran against bring it to eight. Each is worth up to
/// one unit in the last place at the scale of the velocity limit.
constexpr int planner_velocity_rounding_ops = 8;

/// Chained rounding operations behind one acceleration the third-order planner
/// reports. Unlike the second-order planner, which assigns the limit itself, this
/// one integrates the acceleration across the constant-jerk phases it has already
/// traversed: each full phase contributes a product and a sum, and the profile
/// carries at most eleven of them, with a further product and sum for the partial
/// phase the sample lands in. Each is worth up to one unit in the last place at
/// the scale of the acceleration limit.
constexpr int planner_acceleration_rounding_ops = 2 * 11 + 2;

/// Chained rounding operations behind one interior-knot value of the smoothing
/// spline's own defining relation: three smoothed positions and three second
/// derivatives, each read back through a Horner evaluation of four operations, and
/// the six-term relation itself, which chains two differences and two divisions on
/// one side and three products, two sums and three divisions on the other. Six
/// evaluations of four is twenty-four, and the relation's fourteen bring it to
/// thirty-eight. Each is worth up to one unit in the last place at the scale of
/// the largest term the relation forms.
constexpr int moment_relation_rounding_ops = 6 * 4 + 14;

/// Chained rounding operations behind one endpoint acceleration of a cyclic cubic
/// spline solve over `spans` spans. The Sherman-Morrison reduction runs two Thomas
/// solves over the system: each forward-sweep row chains a division, a product and
/// a difference for the diagonal and the same three for the right-hand side, and
/// each back-substitution row a product, a difference and a division -- nine per
/// row across the two solves and the reduction's own recombination adds eight. The
/// endpoint acceleration is then formed from those velocities in six operations
/// for the quadratic coefficient and six for the cubic, and four more in the
/// Horner evaluation; the two endpoints being compared each carry that sixteen.
/// Each is worth up to one unit in the last place at the scale of the largest
/// acceleration the spline reports.
constexpr auto cyclic_spline_acceleration_ops(std::size_t spans) -> int
{
    return 9 * static_cast<int>(spans) + 8 + 2 * 16;
}

/// Chained rounding operations behind one spline value evaluated away from the
/// knot that stores it: three multiply-adds in the Horner form, and the quadratic
/// and cubic coefficients it runs on, which chain six operations each. Each is
/// worth up to one unit in the last place at the scale of the waypoint positions.
constexpr int spline_horner_rounding_ops = 3 * 2 + 6 + 6;

/// The kinematic envelope a limit-respecting planner guarantees at every sample.
///
/// The acceleration is asserted with NO slack. The second-order planner's
/// evaluation assigns it the limit, its negation, or zero, so a magnitude above
/// the limit is not a rounding event but a different number. The velocity is
/// accumulated across a phase and carries that accumulation's rounding.
void require_envelope_2nd(ctrlpp::trajectory_point<double, 1> const& pt,
                          double v_max,
                          double a_max)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();

    double const v_tol = static_cast<double>(planner_velocity_rounding_ops) * eps * v_max;
    CAPTURE(pt.velocity(0), pt.acceleration(0), v_max, a_max, v_tol);
    REQUIRE(std::abs(pt.velocity(0)) <= v_max + v_tol);
    REQUIRE(std::abs(pt.acceleration(0)) <= a_max);
}

/// The same envelope for the third-order planner, whose acceleration is
/// integrated across constant-jerk phases rather than assigned, so it carries a
/// counted budget where the second-order planner's carries none.
void require_envelope_3rd(ctrlpp::trajectory_point<double, 1> const& pt,
                          double v_max,
                          double a_max)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();

    double const v_tol = static_cast<double>(planner_velocity_rounding_ops) * eps * v_max;
    double const a_tol = static_cast<double>(planner_acceleration_rounding_ops) * eps * a_max;
    CAPTURE(pt.velocity(0), pt.acceleration(0), v_max, a_max, v_tol, a_tol);
    REQUIRE(std::abs(pt.velocity(0)) <= v_max + v_tol);
    REQUIRE(std::abs(pt.acceleration(0)) <= a_max + a_tol);
}

/// Assert that the reported acceleration slewed no faster than the jerk limit
/// allows over the step between two samples.
///
/// Within one constant-jerk phase the change is exactly the jerk times the step.
/// A step that crosses a phase boundary splits into two constant-jerk pieces whose
/// durations sum to the step, so the change is a convex combination of two jerks
/// and cannot exceed the limit either. This holds only inside one plan: an update
/// replaces the profile and the acceleration may step discontinuously across it.
void require_jerk_step(double a_now, double a_prev, double dt, double j_max, double a_max)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();

    double const a_tol = 2.0 * static_cast<double>(planner_acceleration_rounding_ops) * eps * a_max;
    CAPTURE(a_now, a_prev, dt, j_max, a_tol);
    REQUIRE(std::abs(a_now - a_prev) <= j_max * dt + a_tol);
}

/// Assert the relation that defines the spline the smoothing solve builds.
///
/// The solve produces smoothed positions and interior second derivatives that
/// satisfy the classical cubic-spline moment relation between them, with the
/// endpoint second derivatives held at zero. The spline stores the smoothed
/// position of span i as that span's constant term and its second derivative as
/// twice the quadratic one, so both sides of the relation are readable from
/// `evaluate` at the knots and no part of the solve is repeated here. A
/// coefficient built from the wrong smoothed positions, or second derivatives that
/// do not belong to them, breaks this relation while still returning numbers.
void require_moment_relation(ctrlpp::smoothing_spline<double> const& spline,
                             std::vector<double> const& times)
{
    constexpr double eps = std::numeric_limits<double>::epsilon();

    // The natural-like endpoint conditions are assigned, not solved: the second
    // derivative outside the interior is zero, and the first span's quadratic
    // coefficient is half of it.
    REQUIRE(spline.evaluate(times.front()).acceleration(0) == 0.0);

    for (std::size_t j = 0; j + 2 < times.size(); ++j) {
        double const h0 = times[j + 1] - times[j];
        double const h1 = times[j + 2] - times[j + 1];

        auto const p0 = spline.evaluate(times[j]);
        auto const p1 = spline.evaluate(times[j + 1]);
        auto const p2 = spline.evaluate(times[j + 2]);

        double const second_difference =
            (p0.position(0) - p1.position(0)) / h0 + (p2.position(0) - p1.position(0)) / h1;
        double const weighted_moments = h0 * p0.acceleration(0) / 6.0
                                        + (h0 + h1) * p1.acceleration(0) / 3.0
                                        + h1 * p2.acceleration(0) / 6.0;

        // Both sides of the relation cancel heavily -- the left is a difference of
        // two nearly equal slopes -- so the rounding is measured at the scale of
        // the operands that entered, not at the scale of what came out. A smoothed
        // position rounds at its own magnitude and is then divided by a span; a
        // second derivative rounds at its own magnitude and is multiplied by one.
        double const position_scale =
            std::max({std::abs(p0.position(0)), std::abs(p1.position(0)),
                      std::abs(p2.position(0))});
        double const moment_scale =
            std::max({std::abs(p0.acceleration(0)), std::abs(p1.acceleration(0)),
                      std::abs(p2.acceleration(0))});
        double const scale =
            position_scale / std::min(h0, h1) + moment_scale * (h0 + h1);
        double const tol = static_cast<double>(moment_relation_rounding_ops) * eps * scale;
        CAPTURE(j, second_difference, weighted_moments, scale, tol);
        REQUIRE(std::abs(second_difference - weighted_moments) <= tol);
    }
}

/// The largest acceleration magnitude the spline reports at its knots.
///
/// The acceleration is affine inside every span, so its extreme over the spline is
/// attained at a knot and this is the whole range rather than a sample of it. It
/// is the scale the endpoint-condition residuals are measured against: those
/// residuals are the rounding of the expressions that produced these values.
auto knot_acceleration_scale(ctrlpp::cubic_spline<double> const& spline,
                             std::vector<double> const& times) -> double
{
    double scale = 0.0;
    for (auto const time : times) {
        scale = std::max(scale, std::abs(spline.evaluate(time).acceleration(0)));
    }
    return scale;
}

/// The straight line that minimizes the squared deviation from the waypoints,
/// evaluated at `at`. It is the limit a smoothing spline tends to as its
/// smoothing weight grows without bound: the smoothed positions become the
/// projection of the waypoints onto the null space of the second-difference
/// operator, which is exactly the affine functions of time.
auto least_squares_line_at(std::vector<double> const& times,
                           std::vector<double> const& positions,
                           double at) -> double
{
    double const n = static_cast<double>(times.size());
    double t_mean = 0.0;
    double q_mean = 0.0;
    for (std::size_t i = 0; i < times.size(); ++i) {
        t_mean += times[i];
        q_mean += positions[i];
    }
    t_mean /= n;
    q_mean /= n;

    double s_tt = 0.0;
    double s_tq = 0.0;
    for (std::size_t i = 0; i < times.size(); ++i) {
        s_tt += (times[i] - t_mean) * (times[i] - t_mean);
        s_tq += (times[i] - t_mean) * (positions[i] - q_mean);
    }
    double const slope = s_tq / s_tt;
    return q_mean + slope * (at - t_mean);
}

}

// ── Cubic spline hardening ─────────────────────────────────────────────────────

TEST_CASE("Cubic spline with exactly 2 points is the straight line through them",
          "[cubic_spline][hardening][coverage]")
{
    // Two waypoints leave the natural solve a two-by-two system whose two rows are
    // the same equation, so both endpoint velocities come out equal to the slope
    // and the quadratic and cubic coefficients vanish. The spline is the straight
    // line, and its midpoint value and slope are the line's, not an approximation
    // of them.
    ctrlpp::cubic_spline<double>::config cfg{
        .times = {0.0, 1.0},
        .positions = {0.0, 1.0},
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);
    auto pt = spline.evaluate(0.5);
    REQUIRE(pt.position(0) == 0.5);
    REQUIRE(pt.velocity(0) == 1.0);
    REQUIRE(pt.acceleration(0) == 0.0);

    // The duration is the difference of the outer knots, formed by one subtraction
    // of two exactly representable operands.
    REQUIRE(spline.duration() == 1.0);
}

TEST_CASE("Cubic spline interpolation matches at knots", "[cubic_spline][hardening][precision]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> positions{0.0, 1.0, 0.5, 2.0, 1.5};

    ctrlpp::cubic_spline<double>::config cfg{
        .times = times,
        .positions = positions,
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);

    // Interpolation at a knot is not approximate. The span containing the knot
    // starts there, so the local time is exactly zero and the Horner form reduces
    // to the span's stored constant term, which is the waypoint. Nothing rounds.
    // The last knot is the exception and it is asserted separately below: the span
    // search clamps it into the final span, where the local time is that span's
    // whole width and the full Horner chain runs.
    for (std::size_t i = 0; i + 1 < times.size(); ++i) {
        auto pt = spline.evaluate(times[i]);
        CAPTURE(i);
        REQUIRE(pt.position(0) == positions[i]);
    }

    constexpr double eps = std::numeric_limits<double>::epsilon();
    double const scale = *std::max_element(positions.begin(), positions.end());
    double const tol = static_cast<double>(spline_horner_rounding_ops) * eps * scale;
    auto const last = spline.evaluate(times.back());
    CAPTURE(last.position(0), positions.back(), tol);
    REQUIRE(std::abs(last.position(0) - positions.back()) <= tol);
}

TEST_CASE("Cubic spline value does not degrade with the knot span",
          "[cubic_spline][hardening][coverage]")
{
    // The knots span fifteen decades. A spline that has lost every significant
    // digit of its coefficients still returns a finite number, so finiteness sees
    // nothing here; the conditioning claim is that the VALUE holds, and it is made
    // by asserting the analytic value with a budget that does NOT grow with the
    // span. A budget proportional to the span would be a third of the value being
    // asserted at the span this case uses, which is no assertion at all.
    //
    // Natural boundary conditions on {0, H, 2H} through {0, 1, 0} put the endpoint
    // second derivatives at zero and, by the antisymmetry of the data about the
    // middle knot, the middle velocity at zero and the outer ones at plus and
    // minus 3 / (2H). The Hermite form on the first span at its midpoint is then
    //   q(H/2) = h01(1/2) * 1 + h10(1/2) * H * (3 / (2H)) = 1/2 + 1/8 * 3/2 = 11/16
    // with no dependence on H whatever. The span cancels out of the answer, and
    // this case asserts that it cancels out of the arithmetic too.
    constexpr double analytic_midpoint = 11.0 / 16.0;
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double const tol = static_cast<double>(spline_horner_rounding_ops) * eps;

    ctrlpp::cubic_spline<double>::config cfg{
        .times = {0.0, 1e15, 2e15},
        .positions = {0.0, 1.0, 0.0},
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);
    auto pt = spline.evaluate(0.5e15);
    CAPTURE(pt.position(0), analytic_midpoint, tol);
    REQUIRE(std::abs(pt.position(0) - analytic_midpoint) <= tol);

    // The same statement swept across the spans the coefficients are formed over.
    // The half-width squared divides the cubic coefficient, so the representable
    // range of that square is the domain this holds on; the sweep stops three
    // decades inside it on either side, which is where the case's own span sits.
    for (double const span : {1e-100, 1e-15, 1.0, 1e6, 1e15, 1e100}) {
        ctrlpp::cubic_spline<double>::config swept{
            .times = {0.0, span, 2.0 * span},
            .positions = {0.0, 1.0, 0.0},
        };
        auto wide = realizable<ctrlpp::cubic_spline<double>>(swept);
        auto const mid = wide.evaluate(0.5 * span);
        CAPTURE(span, mid.position(0));
        REQUIRE(std::abs(mid.position(0) - analytic_midpoint) <= tol);
    }
}

TEST_CASE("Cubic spline declines a span whose coefficients leave the type",
          "[cubic_spline][hardening][negative]")
{
    // The construction is well posed at every span in this case: the knots are
    // strictly increasing, the waypoints are finite, and the spline exists as a
    // mathematical object. What runs out is the type. For the antisymmetric data
    // below the first span's cubic coefficient is exactly -1 / (2 H^3), so the
    // span at which it stops being representable follows from the type's own
    // limits and is not a number chosen here:
    //
    //   * below the cube root of one half the largest finite value, the
    //     coefficient overflows and every evaluation returns an infinity;
    //   * above the cube root of one half the smallest normal value, it goes
    //     subnormal -- the cubic term keeps a handful of bits and then none, and
    //     the spline answers as a quadratic with a plausible finite number that
    //     is not the curve the waypoints describe.
    //
    // Note which quantity that is. The SQUARE of the span is representable on
    // both sides of both boundaries -- at a span of 1e150 its square is 1e300,
    // comfortably finite -- so a domain test on the square of the spacing would
    // admit exactly the configurations that go on to evaluate wrongly. The
    // coefficient is the quantity that leaves the type, and it is what the
    // construction checks.
    constexpr double largest = std::numeric_limits<double>::max();
    constexpr double smallest_normal = std::numeric_limits<double>::min();
    double const overflow_span = std::cbrt(0.5 / largest);
    double const subnormal_span = std::cbrt(0.5 / smallest_normal);

    auto build = [](double span) {
        return ctrlpp::cubic_spline<double>::create({
            .times = {0.0, span, 2.0 * span},
            .positions = {0.0, 1.0, 0.0},
        });
    };

    for (double const span : {overflow_span / 10.0, subnormal_span * 10.0}) {
        auto const declined = build(span);
        CAPTURE(span);
        REQUIRE_FALSE(declined.has_value());
        REQUIRE(declined.error() == ctrlpp::spline_error::unrepresentable_spline);
    }

    // A decade inside either boundary the same construction is served, and it is
    // served correctly: the midpoint value is the span-independent 11/16 derived
    // in the conditioning case above. The rejection is narrow, not a retreat from
    // the wide range of spans the arithmetic does carry.
    constexpr double analytic_midpoint = 11.0 / 16.0;
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double const tol = static_cast<double>(spline_horner_rounding_ops) * eps;
    for (double const span : {overflow_span * 10.0, subnormal_span / 10.0}) {
        auto const served = build(span);
        CAPTURE(span);
        REQUIRE(served.has_value());
        auto const mid = served->evaluate(0.5 * span);
        CAPTURE(mid.position(0));
        REQUIRE(std::abs(mid.position(0) - analytic_midpoint) <= tol);
    }
}

TEST_CASE("Splines carry a constant-velocity segment at every span they accept",
          "[cubic_spline][smoothing_spline][hardening][coverage]")
{
    // A constant-velocity segment is the most ordinary trajectory there is, and
    // the only unusual thing here is the span. It is also the BEST conditioned
    // input either construction sees: every second derivative is zero, so both
    // curvature coefficients are zero and the spline is the straight line through
    // the waypoints.
    //
    // It is a regression case rather than a coverage case, and this is what it
    // guards. Both numerators are differences of slopes, so on collinear data
    // they are mathematically zero and arrive as the ROUNDING RESIDUAL of that
    // cancellation, not as an exact zero. A representability rule that asks
    // whether the numerator was a normal number then reads that residual as a
    // real quantity, divides it by a large span, finds the result subnormal, and
    // refuses -- refusing the flat spline while accepting a curved one at the same
    // span, which is backwards. What decides it correctly is how far the term
    // REACHES against the resolution of the position it contributes to: a
    // residual reaches a fiftieth of that resolution and a genuine coefficient
    // reaches about 1e14 times it.
    constexpr double eps = std::numeric_limits<double>::epsilon();

    SECTION("interpolating spline")
    {
        for (double const span : {1e-150, 1e-100, 1e-50, 1.0, 1e50, 1e100, 1e150, 1e300}) {
            // Waypoints on a line of slope 1/span, so the midpoint of the first
            // span is exactly one half.
            auto const built = ctrlpp::cubic_spline<double>::create({
                .times = {0.0, span, 2.0 * span},
                .positions = {0.0, 1.0, 2.0},
            });
            CAPTURE(span);
            REQUIRE(built.has_value());

            double const tol = static_cast<double>(spline_horner_rounding_ops) * eps * 2.0;
            auto const mid = built->evaluate(0.5 * span);
            CAPTURE(mid.position(0), tol);
            REQUIRE(std::abs(mid.position(0) - 0.5) <= tol);

            // The curvature of a straight line is zero at every scale, and the
            // evaluation says so to the resolution the positions are known at.
            CAPTURE(mid.acceleration(0));
            REQUIRE(std::abs(mid.acceleration(0)) * span * span <= tol);
        }
    }

    SECTION("smoothing spline")
    {
        // The smoothing spline accepts a narrower band of spans than the
        // interpolating one, and both ends of it are its own. Its system matrix
        // carries a smoothness term that grows with the spacing and a
        // regularization term that grows as the spacing shrinks, and the solve
        // squares every entry, so each end is reached when its own term can no
        // longer be squared. For the tradeoff parameter used here that band runs
        // from about 3e-77 to about 1e154 -- measured, and matching where the
        // construction stops working: a span of 1e-77 returned NaN on this
        // straight line before the domain was narrowed. Inside the band the line
        // is carried exactly as above.
        for (double const span : {1e-70, 1e-50, 1.0, 1e50, 1e100, 1e150}) {
            auto const built = ctrlpp::smoothing_spline<double>::create({
                .times = {0.0, span, 2.0 * span, 3.0 * span},
                .positions = {0.0, 1.0, 2.0, 3.0},
                .mu = 0.5,
            });
            CAPTURE(span);
            REQUIRE(built.has_value());

            double const tol = static_cast<double>(spline_horner_rounding_ops) * eps * 3.0;
            auto const mid = built->evaluate(1.5 * span);
            CAPTURE(mid.position(0), tol);
            REQUIRE(std::abs(mid.position(0) - 1.5) <= tol);
        }
    }
}

// ── Smoothing spline hardening ─────────────────────────────────────────────────

TEST_CASE("Smoothing spline at mu=1 IS the interpolating spline",
          "[smoothing_spline][hardening][coverage]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, 0.5, 2.0};

    ctrlpp::smoothing_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .mu = 1.0,
    };

    auto spline = realizable<ctrlpp::smoothing_spline<double>>(cfg);

    // At mu = 1 the regularization weight 2 (1 - mu) / (3 mu) is zero, so the
    // smoothed positions are the waypoints themselves and the spline interpolates
    // them. It does not approach them: the span containing a knot stores that
    // waypoint as its constant term, and the local time at the knot is zero, so
    // the value is the waypoint bit for bit. The last knot is clamped into the
    // final span and runs the full Horner chain, so it carries the counted budget.
    for (std::size_t i = 0; i + 1 < times.size(); ++i) {
        auto pt = spline.evaluate(times[i]);
        CAPTURE(i);
        REQUIRE(pt.position(0) == positions[i]);
    }

    constexpr double eps = std::numeric_limits<double>::epsilon();
    double const scale = *std::max_element(positions.begin(), positions.end());
    double const tol = static_cast<double>(spline_horner_rounding_ops) * eps * scale;
    auto const last = spline.evaluate(times.back());
    CAPTURE(last.position(0), positions.back(), tol);
    REQUIRE(std::abs(last.position(0) - positions.back()) <= tol);

    require_moment_relation(spline, times);
}

TEST_CASE("Smoothing spline at a large weight tends to the least-squares line",
          "[smoothing_spline][hardening][coverage]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, 0.5, 2.0};

    // A small mu makes the regularization weight large, which drives the interior
    // second derivatives toward zero and the smoothed positions toward the
    // projection of the waypoints onto the affine functions of time -- the
    // least-squares straight line. The case name has always said so and nothing
    // asserted it.
    constexpr double mu = 0.001;
    ctrlpp::smoothing_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .mu = mu,
    };

    auto spline = realizable<ctrlpp::smoothing_spline<double>>(cfg);
    require_moment_relation(spline, times);

    // The departure from that line is first order in mu, because the solve differs
    // from the limit by the regularization term divided by the weight and the
    // weight's reciprocal 3 mu / (2 (1 - mu)) is first order in mu. The
    // coefficient of that first order depends on the data through the
    // second-difference operator's pseudo-inverse and is not computed here.
    // Scaling the interpolating spline's own departure by mu bounds it without
    // that coefficient: at mu the smoothing spline must sit at least a thousandth
    // as far from the line as the mu = 1 spline does. A weight that did not grow
    // as mu shrank -- the tradeoff map inverted, or the weight applied to the wrong
    // term -- leaves the spline near the interpolating value and fails this.
    auto interpolating = realizable<ctrlpp::smoothing_spline<double>>(
        ctrlpp::smoothing_spline<double>::config{
            .times = times, .positions = positions, .mu = 1.0});

    double const probe = 1.5;
    double const line = least_squares_line_at(times, positions, probe);
    double const interpolating_departure =
        std::abs(interpolating.evaluate(probe).position(0) - line);
    double const departure = std::abs(spline.evaluate(probe).position(0) - line);

    CAPTURE(line, departure, interpolating_departure, mu);
    REQUIRE(interpolating_departure > 0.0);
    REQUIRE(departure < mu * interpolating_departure);
}

TEST_CASE("Smoothing spline declines a weight its normal equations cannot carry",
          "[smoothing_spline][hardening][negative]")
{
    // The tradeoff parameter is inside its documented domain here, so the
    // rejection is not about the domain: it is about what the regularized system
    // R + lambda Q'Q can hold. The solve forms sums of squares of that system's
    // entries, so the entries must be square-representable, which puts the bound
    // at the square root of the largest finite value. For unit knot spacing the
    // largest entry of the Gram matrix is 1 + 4 + 1, and the weight is
    // 2 (1 - mu) / (3 mu), so the parameter at which the product reaches the bound
    // is four over that square root. Every quantity in that sentence comes from
    // the type or from the spacing.
    constexpr double largest = std::numeric_limits<double>::max();
    double const parameter_bound = 4.0 / std::sqrt(largest);

    // Two waypoint counts, because the two fail differently and only one of them
    // fails loudly. Four waypoints give a two-by-two system, whose decomposition
    // turns an unsquarable entry into infinite pivots and NaN coefficients. Three
    // waypoints give a one-by-one system, which has nothing to eliminate: it
    // divides by the overflowed entry instead, returns exact zeros for the
    // interior second derivatives, and collapses the smoothed positions onto the
    // raw waypoints. That second one is the dangerous case -- a finite, plausible
    // straight-line-through-the-data answer that is not the least-squares limit
    // the weight asked for, and that nothing downstream can tell from a correct
    // one. Both are declined, and by the same enumerator, because they are the
    // same statement about the type.
    for (std::size_t const waypoints : {std::size_t{3}, std::size_t{4}}) {
        std::vector<double> times;
        std::vector<double> positions;
        for (std::size_t i = 0; i < waypoints; ++i) {
            times.push_back(static_cast<double>(i));
            positions.push_back((i % 2 == 0) ? 0.0 : 1.0);
        }
        positions.back() = 2.0;

        auto const declined = ctrlpp::smoothing_spline<double>::create({
            .times = times, .positions = positions, .mu = parameter_bound * 1e-4});
        CAPTURE(waypoints, parameter_bound);
        REQUIRE_FALSE(declined.has_value());
        REQUIRE(declined.error() == ctrlpp::spline_error::unrepresentable_spline);

        // Four decades inside the bound the same system is served, and the answer
        // it gives is the one the weight asks for: the least-squares line.
        auto const served = ctrlpp::smoothing_spline<double>::create({
            .times = times, .positions = positions, .mu = parameter_bound * 1e4});
        REQUIRE(served.has_value());

        constexpr double eps = std::numeric_limits<double>::epsilon();
        double const probe = times[waypoints / 2];
        double const line = least_squares_line_at(times, positions, probe);
        double const scale = *std::max_element(positions.begin(), positions.end());
        double const tol = static_cast<double>(spline_horner_rounding_ops) * eps * scale;
        CAPTURE(served->evaluate(probe).position(0), line, tol);
        REQUIRE(std::abs(served->evaluate(probe).position(0) - line) <= tol);
    }
}

// ── B-spline hardening ─────────────────────────────────────────────────────────

TEST_CASE("B-spline with insufficient control points", "[bspline][hardening][negative]")
{
    // Degree 3 needs at least 4 control points
    using bspline3 = ctrlpp::bspline_trajectory<double, 3>;
    bspline3::config cfg{
        .control_points = {0.0, 1.0, 2.0}, // Only 3
    };

    auto const result = bspline3::create(cfg);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::too_few_control_points);
}

TEST_CASE("B-spline with non-ascending knot vector", "[bspline][hardening][negative]")
{
    using bspline3 = ctrlpp::bspline_trajectory<double, 3>;
    bspline3::config cfg{
        .control_points = {0.0, 1.0, 2.0, 3.0, 4.0},
        .knot_vector = {0.0, 0.0, 0.0, 0.0, 0.5, 0.3, 1.0, 1.0, 1.0}, // Non-ascending
    };

    auto const result = bspline3::create(cfg);
    REQUIRE_FALSE(result.has_value());
    REQUIRE(result.error() == ctrlpp::spline_error::non_monotonic_knots);
}

// ── Trapezoidal trajectory hardening ───────────────────────────────────────────

TEST_CASE("Trapezoidal with zero distance", "[trapezoidal][hardening][coverage]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 5.0, .q1 = 5.0, .v_max = 1.0, .a_max = 1.0,
    };

    // A rest-to-rest profile over no displacement has every phase duration zero
    // and therefore a total duration of zero, and its only sample is the start
    // position it was given. Both are assigned rather than computed. No epsilon
    // takes part in either, which is the form this file uses for the same
    // quantity where a rescale is rejected on a stationary profile.
    auto const profile = trapezoidal_profile(cfg);
    REQUIRE(profile.duration() == 0.0);
    auto pt = profile.evaluate(0.0);
    REQUIRE(pt.position(0) == 5.0);
    REQUIRE(pt.velocity(0) == 0.0);
}

TEST_CASE("Trapezoidal with negative max velocity", "[trapezoidal][hardening][negative]")
{
    // A negative velocity limit puts the cruise velocity below both boundary
    // velocities, which makes the acceleration phase (v_v - v0) / a negative and
    // sends evaluate() into an undefined clamp. It is out of the domain, so it is
    // a typed rejection rather than a profile that happens to be finite.
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 1.0, .v_max = -1.0, .a_max = 1.0,
    };

    auto const created = ctrlpp::trapezoidal_trajectory<double>::create(cfg);
    REQUIRE_FALSE(created.has_value());
    REQUIRE(created.error() == ctrlpp::trajectory_error::non_positive_velocity_limit);
}

TEST_CASE("Trapezoidal triangle profile reaches correct peak velocity", "[trapezoidal][hardening][precision]")
{
    // Short distance forces triangular profile
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.5, .v_max = 10.0, .a_max = 2.0,
    };

    auto const profile = trapezoidal_profile(cfg);
    REQUIRE(profile.is_triangular());

    // A rest-to-rest triangle's two ramps meet at sqrt(a h). The profile forms it
    // by the same square root of the same product, on operands the type
    // represents exactly, so the reported peak is that value bit for bit rather
    // than a percent either side of it.
    REQUIRE(profile.peak_velocity() == std::sqrt(cfg.a_max * (cfg.q1 - cfg.q0)));

    // The traversal is asserted by quadrature. A position sample at the end
    // reproduces the commanded displacement by construction -- the deceleration
    // branch is written backwards from it -- so it cannot fail for a profile whose
    // velocity integrates to something else, and a tighter tolerance on it is a
    // tighter non-oracle.
    require_swept_displacement(profile, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(profile, cfg.v1, cfg.a_max);
    require_nonnegative_phases(profile);
}

// ── Double-S trajectory hardening ──────────────────────────────────────────────

TEST_CASE("Double-S with zero distance", "[double_s][hardening][coverage]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 3.0, .q1 = 3.0, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0,
    };

    // As for the trapezoidal profile over no displacement: every phase duration is
    // zero, the total is their sum, and the only sample is the assigned start
    // position. Exact on both counts.
    auto const profile = double_s_profile(cfg);
    REQUIRE(profile.duration() == 0.0);
    auto pt = profile.evaluate(0.0);
    REQUIRE(pt.position(0) == 3.0);
    REQUIRE(pt.velocity(0) == 0.0);
}

TEST_CASE("Double-S with negative jerk limit", "[double_s][hardening][negative]")
{
    // The jerk limit divides every jerk-phase duration, so its domain is finite
    // and strictly positive and a negative one is a typed rejection.
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 1.0, .v_max = 1.0, .a_max = 1.0, .j_max = -1.0,
    };

    auto const created = ctrlpp::double_s_trajectory<double>::create(cfg);
    REQUIRE_FALSE(created.has_value());
    REQUIRE(created.error() == ctrlpp::trajectory_error::non_positive_jerk_limit);
}

// ── Online planner 2nd hardening ───────────────────────────────────────────────

TEST_CASE("Online planner 2nd rejects out-of-domain velocity limit",
          "[online_planner_2nd][hardening][negative]")
{
    // v_max divides in the planner math (cruise duration h / v_v), so the
    // domain is finite and strictly positive.
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    for (double const v_max : {0.0, -1.0, nan, inf}) {
        auto const result =
            ctrlpp::online_planner_2nd<double>::create({.v_max = v_max, .a_max = 1.0});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_velocity_limit);
    }
}

TEST_CASE("Online planner 2nd rejects out-of-domain acceleration limit",
          "[online_planner_2nd][hardening][negative]")
{
    // a_max divides in the planner math (stopping distance v^2 / (2 a_max),
    // ramp durations v_v / a_max), so the domain is finite and strictly positive.
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    for (double const a_max : {0.0, -1.0, nan, inf}) {
        auto const result =
            ctrlpp::online_planner_2nd<double>::create({.v_max = 1.0, .a_max = a_max});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_acceleration_limit);
    }
}

TEST_CASE("Online planner 2nd with instant target flip", "[online_planner_2nd][hardening][coverage]")
{
    constexpr double target = -5.0;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    REQUIRE(planner.update(5.0).has_value());
    planner.sample(0.1);
    planner.sample(0.2);

    // Flip target mid-motion. The velocity carried into the flip points away from
    // the new target, so the commanded shape does not exist and the planner
    // substitutes a brake-then-replan. It reports which, so the branch this case
    // reaches is asserted rather than assumed.
    REQUIRE(planner.update(target).has_value());
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(planner.diagnostics().substitution_reason
            == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
    REQUIRE(planner.diagnostics().commanded_target == target);
    REQUIRE(planner.diagnostics().brake_duration > 0.0);

    // The substitution costs time, not correctness: the recovery respects both
    // limits at every sample and ends at the commanded target. That is the whole
    // contract of a limit-respecting replan, and a single finite sample sees none
    // of it.
    double t = 0.2;
    for (int i = 0; i < 2000; ++i) {
        t += 0.01;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }

    auto const settled = planner.sample(t + 0.01);
    REQUIRE(planner.is_settled());
    REQUIRE(settled.position(0) == target);
    REQUIRE(settled.velocity(0) == 0.0);
    REQUIRE(settled.acceleration(0) == 0.0);
}

TEST_CASE("Online planner 2nd reaches target", "[online_planner_2nd][hardening][convergence]")
{
    constexpr double target = 3.0;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 1.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    REQUIRE(planner.update(target).has_value());
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::commanded_profile);

    double t = 0.0;
    for (int i = 0; i < 1000; ++i) {
        t += 0.01;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }

    // Past the end of the profile the evaluation returns the stored target with
    // zero velocity: it does not evaluate a polynomial there. Exact equality is
    // what says so, and it is what would catch an evaluation that did.
    auto pt = planner.sample(t + 0.01);
    REQUIRE(planner.is_settled());
    REQUIRE(pt.position(0) == target);
    REQUIRE(pt.velocity(0) == 0.0);
}

// ── Online planner 3rd hardening ───────────────────────────────────────────────

TEST_CASE("Online planner 3rd rejects out-of-domain velocity limit",
          "[online_planner_3rd][hardening][negative]")
{
    // v_max divides in the planner math (cruise duration h / v_max), so the
    // domain is finite and strictly positive.
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    for (double const v_max : {0.0, -1.0, nan, inf}) {
        auto const result = ctrlpp::online_planner_3rd<double>::create(
            {.v_max = v_max, .a_max = 1.0, .j_max = 1.0});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_velocity_limit);
    }
}

TEST_CASE("Online planner 3rd rejects out-of-domain acceleration limit",
          "[online_planner_3rd][hardening][negative]")
{
    // a_max divides in the planner math (constant-deceleration duration
    // |v| / a_max), so the domain is finite and strictly positive.
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    for (double const a_max : {0.0, -1.0, nan, inf}) {
        auto const result = ctrlpp::online_planner_3rd<double>::create(
            {.v_max = 1.0, .a_max = a_max, .j_max = 1.0});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_acceleration_limit);
    }
}

TEST_CASE("Online planner 3rd rejects out-of-domain jerk limit",
          "[online_planner_3rd][hardening][negative]")
{
    // j_max divides in the planner math (jerk-phase durations a_max / j_max
    // and |a| / j_max), so the domain is finite and strictly positive.
    auto constexpr nan = std::numeric_limits<double>::quiet_NaN();
    auto constexpr inf = std::numeric_limits<double>::infinity();

    for (double const j_max : {0.0, -1.0, nan, inf}) {
        auto const result = ctrlpp::online_planner_3rd<double>::create(
            {.v_max = 1.0, .a_max = 1.0, .j_max = j_max});
        REQUIRE_FALSE(result.has_value());
        REQUIRE(result.error() == ctrlpp::trajectory_error::non_positive_jerk_limit);
    }
}

TEST_CASE("Online planner 3rd with instant target reversal",
          "[online_planner_3rd][hardening][coverage]")
{
    constexpr double target = -5.0;
    constexpr double step = 0.01;
    ctrlpp::online_planner_3rd<double>::config cfg{.v_max = 1.0, .a_max = 2.0, .j_max = 5.0};
    auto planner = realizable<ctrlpp::online_planner_3rd<double>>(cfg);

    REQUIRE(planner.update(5.0).has_value());
    planner.sample(0.1);

    // Reversing the commanded direction puts the carried velocity outside the
    // domain of the shape that would carry it through, so the planner brakes to
    // rest and replans from the stopping point. It reports that it did.
    REQUIRE(planner.update(target).has_value());
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(planner.diagnostics().substitution_reason
            == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
    REQUIRE(planner.diagnostics().brake_duration > 0.0);

    // The third-order planner bounds a jerk as well, so the envelope has a third
    // member: the reported acceleration may not slew faster than the jerk limit
    // over a sampling step. All three are checked at every sample of one plan.
    double t = 0.1;
    double previous_acceleration = planner.sample(t).acceleration(0);
    for (int i = 0; i < 2000; ++i) {
        t += step;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_3rd(pt, cfg.v_max, cfg.a_max);
        require_jerk_step(pt.acceleration(0), previous_acceleration, step, cfg.j_max, cfg.a_max);
        previous_acceleration = pt.acceleration(0);
    }

    auto const settled = planner.sample(t + step);
    REQUIRE(planner.is_settled());
    REQUIRE(settled.position(0) == target);
    REQUIRE(settled.velocity(0) == 0.0);
    REQUIRE(settled.acceleration(0) == 0.0);
}

TEST_CASE("Online planner 3rd reaches target", "[online_planner_3rd][hardening][convergence]")
{
    constexpr double target = 3.0;
    constexpr double step = 0.01;
    ctrlpp::online_planner_3rd<double>::config cfg{.v_max = 2.0, .a_max = 1.0, .j_max = 5.0};
    auto planner = realizable<ctrlpp::online_planner_3rd<double>>(cfg);

    REQUIRE(planner.update(target).has_value());
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::commanded_profile);

    double t = 0.0;
    double previous_acceleration = 0.0;
    for (int i = 0; i < 2000; ++i) {
        t += step;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_3rd(pt, cfg.v_max, cfg.a_max);
        require_jerk_step(pt.acceleration(0), previous_acceleration, step, cfg.j_max, cfg.a_max);
        previous_acceleration = pt.acceleration(0);
    }

    auto pt = planner.sample(t + step);
    REQUIRE(planner.is_settled());
    REQUIRE(pt.position(0) == target);
    REQUIRE(pt.velocity(0) == 0.0);
    REQUIRE(pt.acceleration(0) == 0.0);
}

// ── Coverage gap-filling tests ────────────────────────────────────────────────

TEST_CASE("Smoothing spline with 2 points degenerates to linear",
          "[smoothing_spline][hardening][coverage]")
{
    ctrlpp::smoothing_spline<double>::config cfg{
        .times = {0.0, 1.0},
        .positions = {0.0, 5.0},
        .mu = 0.5,
    };
    // Two waypoints leave no interior knot for the smoothing solve to act on, so
    // the construction takes the linear path: the span's constant term is the
    // first position, its linear term the slope, and its quadratic and cubic terms
    // are assigned zero. The tradeoff parameter never enters. The midpoint value,
    // the slope and the vanishing curvature are therefore exact, not approximate.
    auto spline = realizable<ctrlpp::smoothing_spline<double>>(cfg);

    auto pt = spline.evaluate(0.5);
    REQUIRE(pt.position(0) == 2.5);
    REQUIRE(pt.velocity(0) == 5.0);
    REQUIRE(pt.acceleration(0) == 0.0);
    REQUIRE(spline.duration() == 1.0);
}

TEST_CASE("Smoothing spline rejects non-finite two-point geometry",
          "[smoothing_spline][hardening][negative]")
{
    auto const max = std::numeric_limits<double>::max();

    SECTION("finite endpoints with an overflowing span")
    {
        auto const rejected = ctrlpp::smoothing_spline<double>::create({
            .times = {-max, max},
            .positions = {0.0, 1.0},
            .mu = 0.5,
        });
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::spline_error::unrepresentable_spline);
    }

    SECTION("an infinite endpoint")
    {
        auto const rejected = ctrlpp::smoothing_spline<double>::create({
            .times = {0.0, std::numeric_limits<double>::infinity()},
            .positions = {0.0, 1.0},
            .mu = 0.5,
        });
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::spline_error::non_finite_input);
    }

    SECTION("a non-finite waypoint")
    {
        auto const rejected = ctrlpp::smoothing_spline<double>::create({
            .times = {0.0, 1.0},
            .positions = {0.0, std::numeric_limits<double>::quiet_NaN()},
            .mu = 0.5,
        });
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::spline_error::non_finite_input);
    }
}

TEST_CASE("Trapezoidal trajectory rescale_to extends motion",
          "[trapezoidal][hardening][coverage]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 2.0,
    };
    auto profile = trapezoidal_profile(cfg);
    auto const original_T = profile.duration();

    auto const rescaled = profile.rescale_to(original_T * 2.0);
    REQUIRE(rescaled.has_value());
    REQUIRE(profile.duration() > original_T);
    require_realized_duration(profile, original_T * 2.0);

    // The traversed displacement is asserted by integrating the reported velocity
    // over the profile's own phase segments. A position sample at or near the end
    // is NOT a valid check: the final segment is written as an offset backwards
    // from the commanded displacement, so it returns the target position even on a
    // profile whose velocity integrates to something else entirely.
    require_swept_displacement(profile, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(profile, cfg.v1, cfg.a_max);
    require_nonnegative_phases(profile);
    REQUIRE(std::abs(profile.peak_velocity()) <= cfg.v_max);
}

TEST_CASE("Trapezoidal trajectory preserves representable large-velocity motion",
          "[trapezoidal][hardening][coverage]")
{
    SECTION("equal boundary velocities produce a constant-velocity segment")
    {
        ctrlpp::trapezoidal_trajectory<double>::config cfg{
            .q0 = 0.0,
            .q1 = 1e200,
            .v_max = 2e200,
            .a_max = 1.0,
            .v0 = 1e200,
            .v1 = 1e200,
        };

        auto const built = ctrlpp::trapezoidal_trajectory<double>::create(cfg);
        REQUIRE(built.has_value());
        CHECK(built->duration() == 1.0);
        CHECK(built->peak_velocity() == cfg.v0);
        CHECK(built->phase_durations()[0] == 0.0);
        CHECK(built->phase_durations()[1] == 1.0);
        CHECK(built->phase_durations()[2] == 0.0);
        require_swept_displacement(*built, cfg.q1 - cfg.q0, cfg.a_max);
    }

    SECTION("nearby boundary velocities do not overflow the feasibility test")
    {
        for(int ulps = 1; ulps <= 8; ++ulps)
        {
            double terminal_velocity = 1e150;
            for(int step = 0; step < ulps; ++step)
            {
                terminal_velocity =
                    std::nextafter(terminal_velocity,
                                   std::numeric_limits<double>::infinity());
            }

            ctrlpp::trapezoidal_trajectory<double>::config cfg{
                .q0 = 0.0,
                .q1 = 1e150,
                .v_max = 2e150,
                .a_max = 1e150,
                .v0 = 1e150,
                .v1 = terminal_velocity,
            };
            auto const built =
                ctrlpp::trapezoidal_trajectory<double>::create(cfg);
            CAPTURE(ulps, terminal_velocity);
            REQUIRE(built.has_value());
            REQUIRE(std::isfinite(built->duration()));
            REQUIRE(std::isfinite(built->peak_velocity()));
        }
    }
}

TEST_CASE("Trapezoidal trajectory rescale_to shorter than current is rejected",
          "[trapezoidal][hardening][negative]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 2.0,
    };
    auto profile = trapezoidal_profile(cfg);
    auto const original_T = profile.duration();
    auto const original_phases = profile.phase_durations();

    auto const rescaled = profile.rescale_to(original_T * 0.5);
    REQUIRE(!rescaled.has_value());
    REQUIRE(rescaled.error() == ctrlpp::trajectory_error::duration_shorter_than_current);

    // Nothing may have moved: a rejected request leaves the profile untouched.
    REQUIRE(profile.duration() == original_T);
    REQUIRE(profile.phase_durations() == original_phases);
}

TEST_CASE("Trapezoidal trajectory rescale_to emits the valley shape",
          "[trapezoidal][hardening][coverage]")
{
    // Both boundary velocities sit well above the cruise velocity a long duration
    // needs, so the profile decelerates away from the initial velocity, holds a
    // low cruise velocity, and accelerates back up to the final one. Solving the
    // plateau branch here would hand back a negative acceleration phase.
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 1.0, .v_max = 2.0, .a_max = 1.0, .v0 = 0.9, .v1 = 0.9,
    };
    auto profile = trapezoidal_profile(cfg);

    auto const rescaled = profile.rescale_to(10.0);
    REQUIRE(rescaled.has_value());

    auto const phases = profile.phase_durations();
    REQUIRE(phases[0] > 0.0);
    REQUIRE(phases[1] > 0.0);
    REQUIRE(phases[2] > 0.0);
    REQUIRE(std::abs(profile.peak_velocity()) < cfg.v0);

    require_realized_duration(profile, 10.0);
    require_swept_displacement(profile, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(profile, cfg.v1, cfg.a_max);
}

TEST_CASE("Double-S trajectory rescale_to rebuilds a rest-to-rest profile",
          "[double_s][hardening][coverage]")
{
    // The family the old cruise-padding rescale corrupted in every measured case,
    // and the one where the scale follows in closed form because the duration is
    // exactly proportional to the reciprocal of the scale at rest.
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0,
    };
    auto profile = double_s_profile(cfg);
    auto const original_T = profile.duration();

    auto const rescaled = profile.rescale_to(original_T * 3.0);
    REQUIRE(rescaled.has_value());
    require_realized_duration(profile, original_T * 3.0);
    require_swept_displacement(profile, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(profile, cfg.v1, cfg.a_max);
    require_nonnegative_phases(profile);

    // The rebuilt profile respects the scaled limits it was built under.
    REQUIRE(std::abs(profile.peak_velocity()) <= cfg.v_max);
}

TEST_CASE("Double-S trajectory rescale_to rebuilds with nonzero boundary velocities",
          "[double_s][hardening][coverage]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0,
        .v0 = 1.0, .v1 = 0.5,
    };
    auto profile = double_s_profile(cfg);
    auto const original_T = profile.duration();

    auto const rescaled = profile.rescale_to(original_T * 1.5);
    REQUIRE(rescaled.has_value());
    require_realized_duration(profile, original_T * 1.5);
    require_swept_displacement(profile, cfg.q1 - cfg.q0, cfg.a_max);
    require_nonnegative_phases(profile);

    // The boundary velocities are the command and are left unscaled, so the
    // profile still arrives at the one it was given.
    require_terminal_velocity(profile, cfg.v1, cfg.a_max);
}

TEST_CASE("Double-S trajectory rescale_to rejects shortening and unreachable requests",
          "[double_s][hardening][negative]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0,
        .v0 = 1.0, .v1 = 0.5,
    };
    auto profile = double_s_profile(cfg);
    auto const original_T = profile.duration();
    auto const original_phases = profile.phase_durations();

    auto const shorter = profile.rescale_to(original_T * 0.5);
    REQUIRE(!shorter.has_value());
    REQUIRE(shorter.error() == ctrlpp::trajectory_error::duration_shorter_than_current);
    REQUIRE(profile.duration() == original_T);
    REQUIRE(profile.phase_durations() == original_phases);

    // Slowing the profile down means lowering the velocity limit, and the limit
    // cannot fall below the boundary velocities the caller commanded. That pins a
    // finite reachable maximum well short of ten times the current duration.
    auto const unreachable = profile.rescale_to(original_T * 10.0);
    REQUIRE(!unreachable.has_value());
    REQUIRE(unreachable.error() == ctrlpp::trajectory_error::unreachable_duration);
    REQUIRE(profile.duration() == original_T);
    REQUIRE(profile.phase_durations() == original_phases);

    // Its own duration is a success no-op on both profiles.
    REQUIRE(profile.rescale_to(original_T).has_value());
    REQUIRE(profile.duration() == original_T);
}

TEST_CASE("Trapezoidal trajectory rescale_to rejects a duration past the reachable maximum",
          "[trapezoidal][hardening][negative]")
{
    // The commanded displacement is below the boundary-kinetic term
    // (v0^2 + v1^2) / (2 a), so the cruise duration reaches zero at a strictly
    // positive cruise velocity and the reachable durations stop at
    // (v0 + v1 - 2 sqrt((v0^2 + v1^2) / 2 - a h)) / a rather than growing without
    // bound.
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.1, .v_max = 5.0, .a_max = 1.0, .v0 = 1.0, .v1 = 1.0,
    };
    auto profile = trapezoidal_profile(cfg);
    auto const original_T = profile.duration();

    auto const h = std::abs(cfg.q1 - cfg.q0);
    auto const v_min = std::sqrt((cfg.v0 * cfg.v0 + cfg.v1 * cfg.v1) / 2.0 - cfg.a_max * h);
    auto const T_max = (cfg.v0 + cfg.v1 - 2.0 * v_min) / cfg.a_max;
    REQUIRE(T_max > original_T);

    // Inside the reachable set the request is served.
    auto const inside = profile.rescale_to(0.5 * (original_T + T_max));
    REQUIRE(inside.has_value());

    // Beyond it the request is a typed rejection, not a clamped success.
    auto other = trapezoidal_profile(cfg);
    auto const outside = other.rescale_to(T_max * 2.0);
    REQUIRE(!outside.has_value());
    REQUIRE(outside.error() == ctrlpp::trajectory_error::unreachable_duration);
    REQUIRE(other.duration() == original_T);
}

TEST_CASE("Online planner 2nd retargets while moving triggers braking",
          "[online_planner_2nd][hardening][coverage]")
{
    constexpr double target = -5.0;
    constexpr double step = 0.05;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 4.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Start moving to 10
    REQUIRE(planner.update(10.0).has_value());
    // Sample partway through to build up velocity
    for (int i = 0; i < 20; ++i) {
        planner.sample(step * static_cast<double>(i + 1));
    }
    auto const carried = planner.sample(step * 20.0);
    REQUIRE(carried.velocity(0) > 0.0);

    // Retarget behind the direction of travel. The braking the case name claims is
    // a branch of the planner and the planner names it: nothing about the motion
    // alone distinguishes a brake-then-replan from a commanded profile that
    // happens to decelerate first.
    REQUIRE(planner.update(target).has_value());
    auto const& report = planner.diagnostics();
    REQUIRE(report.disposition == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(report.substitution_reason
            == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
    REQUIRE(report.initial_velocity == carried.velocity(0));

    // The brake runs at the acceleration limit from the carried velocity to rest,
    // so its duration is that velocity over that limit and the ground it covers is
    // the mean of the two velocities times it. Both are the closed forms the
    // branch is built from, evaluated on the same operands in the same order.
    REQUIRE(report.brake_duration == std::abs(carried.velocity(0)) / cfg.a_max);
    REQUIRE(report.replan_start_position
            == carried.position(0) + carried.velocity(0) * report.brake_duration / 2.0);

    double t = step * 20.0;
    for (int i = 0; i < 200; ++i) {
        t += step;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }

    auto const settled = planner.sample(t + step);
    REQUIRE(planner.is_settled());
    REQUIRE(settled.position(0) == target);
    REQUIRE(settled.velocity(0) == 0.0);
}

TEST_CASE("Online planner 2nd with same position target is a no-op",
          "[online_planner_2nd][hardening][coverage]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // A target at the current position with no velocity to shed is a zero-duration
    // profile, and every sample of it returns the stored target with zero velocity
    // and zero acceleration. Nothing is computed, so nothing rounds.
    REQUIRE(planner.update(0.0).has_value());
    REQUIRE(planner.diagnostics().disposition == ctrlpp::online_planner_disposition::settled);
    REQUIRE(planner.diagnostics().planned_duration == 0.0);

    auto pt = planner.sample(0.01);
    REQUIRE(pt.position(0) == 0.0);
    REQUIRE(pt.velocity(0) == 0.0);
    REQUIRE(pt.acceleration(0) == 0.0);
}

TEST_CASE("Trapezoidal negative cruise duration clamped to zero",
          "[trapezoidal][hardening][coverage]")
{
    // The velocity limit is far above what the displacement can reach, so the
    // cruise duration the three-phase solve produces is negative and is clamped to
    // zero. The shape predicate is the observable that says the clamp fired: a
    // profile with a cruise phase does not report itself triangular.
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.1, .v_max = 100.0, .a_max = 1.0,
    };
    auto profile = trapezoidal_profile(cfg);

    REQUIRE(profile.is_triangular());
    REQUIRE(profile.phase_durations()[1] == 0.0);

    // With the cruise phase gone the two ramps meet at sqrt(a h), exactly as they
    // do for a triangle the limits never constrained.
    REQUIRE(profile.peak_velocity() == std::sqrt(cfg.a_max * (cfg.q1 - cfg.q0)));

    // The traversal is asserted by quadrature rather than by the endpoint sample
    // the deceleration branch reproduces by construction.
    require_swept_displacement(profile, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(profile, cfg.v1, cfg.a_max);
    require_nonnegative_phases(profile);
}

// ── Online planner 2nd: overshoot and braking coverage ────────────────────────

TEST_CASE("Online planner 2nd overshoot recovery brakes and reverses",
          "[online_planner_2nd][hardening][coverage]")
{
    // Start moving AWAY from target: positive velocity, negative target
    constexpr double target = -3.0;
    constexpr double step = 0.01;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 4.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // First, build up positive velocity toward +10
    REQUIRE(planner.update(10.0).has_value());
    double t = 0.0;
    for (int i = 0; i < 30; ++i) {
        t += step;
        planner.sample(t);
    }
    auto const carried = planner.sample(t);
    REQUIRE(carried.velocity(0) > 0.0);

    // Now retarget behind us: the carried velocity points away from the target, so
    // the commanded shape does not exist and the planner brakes first.
    REQUIRE(planner.update(target).has_value());
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(planner.diagnostics().substitution_reason
            == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
    REQUIRE(planner.diagnostics().brake_duration
            == std::abs(carried.velocity(0)) / cfg.a_max);

    // The recovery is where the kinematic envelope is the contract: a braking
    // reversal that exceeded either limit would be a motion the machine cannot
    // execute, and every sample of it is a number.
    for (int i = 0; i < 200; ++i) {
        t += step;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }

    // Eventually should reach the target
    for (int i = 0; i < 800; ++i) {
        t += step;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }
    auto final_pt = planner.sample(t + step);
    REQUIRE(planner.is_settled());
    REQUIRE(final_pt.position(0) == target);
    REQUIRE(final_pt.velocity(0) == 0.0);
}

TEST_CASE("Online planner 2nd wrong-direction: positive velocity, target behind",
          "[online_planner_2nd][hardening][coverage]")
{
    constexpr double step = 0.01;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 3.0, .a_max = 5.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build up positive velocity by targeting +5
    REQUIRE(planner.update(5.0).has_value());
    double t = 0.0;
    for (int i = 0; i < 20; ++i) {
        t += step;
        planner.sample(t);
    }

    // The precondition this case needs is that the velocity points forward, which
    // is what makes the retarget below a reversal. That is the statement; a
    // threshold fitted to the velocity the planner happens to have built is not.
    auto mid = planner.sample(t);
    REQUIRE(mid.velocity(0) > 0.0);

    // Now target is behind current position (same sign but smaller). The stopping
    // distance exceeds the remaining displacement, so overshoot is detected.
    double const target = mid.position(0) - 0.001;
    REQUIRE(planner.update(target).has_value());
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(planner.diagnostics().substitution_reason
            == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);

    // Sample through -- must brake, stop, and reverse slightly
    for (int i = 0; i < 500; ++i) {
        t += step;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }
    auto settled = planner.sample(t + step);
    REQUIRE(planner.is_settled());
    REQUIRE(settled.velocity(0) == 0.0);
    REQUIRE(settled.position(0) == target);
}

TEST_CASE("Online planner 2nd near-zero displacement with velocity triggers braking",
          "[online_planner_2nd][hardening][coverage]")
{
    constexpr double step = 0.01;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build up velocity
    REQUIRE(planner.update(5.0).has_value());
    double t = 0.0;
    for (int i = 0; i < 50; ++i) {
        t += step;
        planner.sample(t);
    }

    auto current = planner.sample(t);
    REQUIRE(current.velocity(0) > 0.0);

    // Retarget to exactly the current position: the displacement vanishes while
    // the velocity does not, which is the third of the three conditions that
    // select the brake-then-replan branch.
    REQUIRE(planner.update(current.position(0)).has_value());
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(planner.diagnostics().substitution_reason
            == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
    REQUIRE(planner.diagnostics().brake_duration
            == std::abs(current.velocity(0)) / cfg.a_max);

    for (int i = 0; i < 500; ++i) {
        t += step;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }
    REQUIRE(planner.is_settled());

    auto const settled = planner.sample(t + step);
    REQUIRE(settled.position(0) == current.position(0));
    REQUIRE(settled.velocity(0) == 0.0);
}

TEST_CASE("Online planner 2nd evaluate_profile at and past T boundary",
          "[online_planner_2nd][hardening][coverage]")
{
    constexpr double target = 1.0;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    REQUIRE(planner.update(target).has_value());

    // Sample well past when we should have arrived
    double t = 0.0;
    for (int i = 0; i < 2000; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    // Past the end of the profile the evaluation returns the stored target and
    // zeros; it does not run a polynomial out beyond the interval it was fitted
    // on. All three quantities are assigned, so all three are exact -- and an
    // evaluation that ran the polynomial instead would drift and be caught here
    // where a tolerance of a hundred millionth would not notice.
    auto pt = planner.sample(t + 100.0);
    REQUIRE(pt.position(0) == target);
    REQUIRE(pt.velocity(0) == 0.0);
    REQUIRE(pt.acceleration(0) == 0.0);
}

TEST_CASE("Online planner 2nd braking phase is entered and evaluated as a brake",
          "[online_planner_2nd][hardening][coverage]")
{
    constexpr double target = 5.0;
    constexpr double fine_step = 0.005;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 3.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build velocity in negative direction
    REQUIRE(planner.update(-8.0).has_value());
    double t = 0.0;
    for (int i = 0; i < 40; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    auto state = planner.sample(t);
    REQUIRE(state.velocity(0) < 0.0);

    // Now target in the positive direction. Which branch this takes is not a
    // property of the motion -- both branches reach the target under both limits
    // -- so the planner reports it, and that report is what turns a claim about
    // branch coverage into something that can fail.
    REQUIRE(planner.update(target).has_value());
    auto const& report = planner.diagnostics();
    REQUIRE(report.disposition == ctrlpp::online_planner_disposition::braked_and_replanned);
    REQUIRE(report.substitution_reason
            == ctrlpp::online_planner_substitution_reason::reversal_or_overshoot);
    REQUIRE(report.initial_velocity == state.velocity(0));
    REQUIRE(report.brake_duration == std::abs(state.velocity(0)) / cfg.a_max);
    REQUIRE(report.replan_start_position
            == state.position(0) + state.velocity(0) * report.brake_duration / 2.0);

    // Sample inside the braking window, which the fine step stays well within. The
    // brake is a constant deceleration opposing the carried velocity, so its
    // signature is exact: the acceleration is the limit itself, signed against the
    // velocity, and the velocity is the carried one walked back along it. Neither
    // of the other two branches produces that.
    double const brake_start = t;
    REQUIRE(10.0 * fine_step < report.brake_duration);
    for (int i = 0; i < 10; ++i) {
        t += fine_step;
        auto const pt = planner.sample(t);
        double const elapsed = t - brake_start;
        double const expected_velocity = state.velocity(0) + cfg.a_max * elapsed;
        constexpr double eps = std::numeric_limits<double>::epsilon();
        double const tol =
            static_cast<double>(planner_velocity_rounding_ops) * eps * cfg.v_max;
        CAPTURE(i, t, elapsed, pt.velocity(0), expected_velocity, tol);
        REQUIRE(pt.acceleration(0) == cfg.a_max);
        REQUIRE(std::abs(pt.velocity(0) - expected_velocity) <= tol);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }

    // Continue through the replanned rest-to-rest move that follows the brake.
    for (int i = 0; i < 1000; ++i) {
        t += 0.01;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }
    auto final_pt = planner.sample(t + 0.01);
    REQUIRE(planner.is_settled());
    REQUIRE(final_pt.position(0) == target);
    REQUIRE(final_pt.velocity(0) == 0.0);
}

TEST_CASE("Online planner 2nd reset clears state",
          "[online_planner_2nd][hardening][coverage]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 3.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    REQUIRE(planner.update(10.0).has_value());
    planner.sample(0.5);

    // Reset stores the position and zeros everything else, so the sample that
    // follows returns what was stored rather than anything computed from it.
    planner.reset(7.0);
    REQUIRE(planner.is_settled());
    REQUIRE(planner.diagnostics().disposition == ctrlpp::online_planner_disposition::settled);
    REQUIRE(planner.diagnostics().commanded_target == 7.0);
    REQUIRE(planner.diagnostics().planned_duration == 0.0);

    auto pt = planner.sample(0.0);
    REQUIRE(pt.position(0) == 7.0);
    REQUIRE(pt.velocity(0) == 0.0);
    REQUIRE(pt.acceleration(0) == 0.0);
}

// ── Cubic spline: periodic and clamped boundary conditions ────────────────────

TEST_CASE("Cubic spline periodic BC wraps velocity and acceleration",
          "[cubic_spline][hardening][coverage]")
{
    // Periodic BC requires q_0 == q_n
    std::vector<double> times{0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> positions{1.0, 2.0, 0.5, 2.5, 1.0};

    ctrlpp::cubic_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .bc = ctrlpp::boundary_condition::periodic,
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);

    // Matching endpoint derivatives are not an outcome of the solve, they are the
    // constraint the cyclic system imposes: the last velocity is assigned the
    // first, and the closing row of the system equates the accelerations. What is
    // left over is the rounding of the expressions that carried those velocities
    // into the two end spans' coefficients, so one counted budget covers all three
    // quantities rather than three unexplained magnitudes.
    constexpr double eps = std::numeric_limits<double>::epsilon();
    auto const spans = times.size() - 1;
    double const acceleration_scale = knot_acceleration_scale(spline, times);
    double const position_scale = *std::max_element(positions.begin(), positions.end());
    double const derivative_tol = static_cast<double>(cyclic_spline_acceleration_ops(spans))
                                  * eps * acceleration_scale;

    auto start = spline.evaluate(times.front());
    auto end = spline.evaluate(times.back());
    CAPTURE(acceleration_scale, derivative_tol);
    REQUIRE(std::abs(start.velocity(0) - end.velocity(0)) <= derivative_tol);
    REQUIRE(std::abs(start.acceleration(0) - end.acceleration(0)) <= derivative_tol);

    // The first knot's value is the span's stored constant term, so the start
    // position is the waypoint exactly; the last is clamped into the final span
    // and runs the Horner chain, so it carries the smaller budget for that chain.
    REQUIRE(start.position(0) == positions.front());
    double const horner_tol =
        static_cast<double>(spline_horner_rounding_ops) * eps * position_scale;
    CAPTURE(end.position(0), horner_tol);
    REQUIRE(std::abs(end.position(0) - positions.back()) <= horner_tol);
}

TEST_CASE("Cubic spline periodic endpoint comparison is overflow safe",
          "[cubic_spline][hardening][negative]")
{
    SECTION("a large finite mismatch is rejected")
    {
        auto const rejected = ctrlpp::cubic_spline<double>::create({
            .times = {0.0, 1.0, 2.0},
            .positions = {1e308, 9.5e307, 9e307},
            .bc = ctrlpp::boundary_condition::periodic,
        });
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::spline_error::periodic_endpoint_mismatch);
    }

    SECTION("non-finite endpoints are rejected before comparison")
    {
        auto const rejected = ctrlpp::cubic_spline<double>::create({
            .times = {0.0, 1.0, 2.0},
            .positions = {std::numeric_limits<double>::infinity(), 0.0,
                          std::numeric_limits<double>::infinity()},
            .bc = ctrlpp::boundary_condition::periodic,
        });
        REQUIRE_FALSE(rejected.has_value());
        CHECK(rejected.error() == ctrlpp::spline_error::non_finite_input);
    }

    SECTION("large matching endpoints remain closed")
    {
        auto const built = ctrlpp::cubic_spline<double>::create({
            .times = {0.0, 1.0, 2.0},
            .positions = {1e150, 0.0, 1e150},
            .bc = ctrlpp::boundary_condition::periodic,
        });
        REQUIRE(built.has_value());
        auto const start = built->evaluate(0.0);
        auto const end = built->evaluate(2.0);
        double const tolerance =
            static_cast<double>(spline_horner_rounding_ops)
            * std::numeric_limits<double>::epsilon() * 1e150;
        CHECK(start.position(0) == 1e150);
        CHECK(std::abs(end.position(0) - 1e150) <= tolerance);
    }
}

TEST_CASE("Cubic spline clamped BC with exactly 2 waypoints",
          "[cubic_spline][hardening][coverage]")
{
    ctrlpp::cubic_spline<double>::config cfg{
        .times = {0.0, 1.0},
        .positions = {0.0, 1.0},
        .bc = ctrlpp::boundary_condition::clamped,
        .v0 = 2.0,
        .vn = -1.0,
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);

    // With two waypoints the clamped solve has no interior unknowns at all: both
    // endpoint velocities are the supplied ones, assigned rather than solved, and
    // the first span's constant term is the first waypoint. Three of the four
    // quantities below are therefore exact. Only the last knot's position runs the
    // Horner chain, because the span search clamps it into the single span.
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double const tol = static_cast<double>(spline_horner_rounding_ops) * eps
                       * std::max(std::abs(cfg.positions.back()), std::abs(cfg.vn));

    auto start = spline.evaluate(0.0);
    auto end = spline.evaluate(1.0);
    REQUIRE(start.position(0) == 0.0);
    REQUIRE(start.velocity(0) == cfg.v0);
    CAPTURE(end.position(0), end.velocity(0), tol);
    REQUIRE(std::abs(end.position(0) - 1.0) <= tol);
    REQUIRE(std::abs(end.velocity(0) - cfg.vn) <= tol);
}

TEST_CASE("Cubic spline clamped BC with interior knots",
          "[cubic_spline][hardening][coverage]")
{
    ctrlpp::cubic_spline<double>::config cfg{
        .times = {0.0, 1.0, 2.0, 3.0},
        .positions = {0.0, 1.0, 0.5, 2.0},
        .bc = ctrlpp::boundary_condition::clamped,
        .v0 = 0.5,
        .vn = 1.0,
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);

    // The endpoint velocities are imposed, not solved: the clamped path assigns
    // them before reducing the interior system, so the first one survives into the
    // first span's linear coefficient untouched. The last one reaches the sample
    // through the final span's Horner chain, which is what its budget covers, and
    // the same budget covers the last knot's position for the same reason.
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double const scale = std::max(*std::max_element(cfg.positions.begin(), cfg.positions.end()),
                                  std::abs(cfg.vn));
    double const tol = static_cast<double>(spline_horner_rounding_ops) * eps * scale;

    auto start = spline.evaluate(0.0);
    auto end = spline.evaluate(3.0);
    REQUIRE(start.velocity(0) == cfg.v0);
    CAPTURE(end.velocity(0), tol);
    REQUIRE(std::abs(end.velocity(0) - cfg.vn) <= tol);

    // Knot interpolation, exact at every knot the span search does not clamp.
    for (std::size_t i = 0; i + 1 < cfg.times.size(); ++i) {
        auto pt = spline.evaluate(cfg.times[i]);
        CAPTURE(i);
        REQUIRE(pt.position(0) == cfg.positions[i]);
    }
    CAPTURE(end.position(0));
    REQUIRE(std::abs(end.position(0) - cfg.positions.back()) <= tol);
}

TEST_CASE("Cubic spline find_span at exact knot time returns correct span",
          "[cubic_spline][hardening][coverage]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, 0.0, 1.0};

    ctrlpp::cubic_spline<double>::config cfg{
        .times = times,
        .positions = positions,
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);

    // The span search has no accessor, so the span it chose is asserted through
    // what that choice produces. At an interior knot the correct span is the one
    // that STARTS there: its local time is zero and the Horner form collapses to
    // its stored constant term, so the value is the waypoint bit for bit. The
    // preceding span reaches the same waypoint only through its full cubic chain,
    // which reproduces it to within rounding and not generally to the bit. Exact
    // equality is therefore the sharpest statement available about which span ran,
    // and it is why the exact form is used here rather than a tolerance that both
    // spans would satisfy.
    for (std::size_t i = 0; i + 1 < times.size(); ++i) {
        auto pt = spline.evaluate(times[i]);
        CAPTURE(i);
        REQUIRE(pt.position(0) == positions[i]);
    }

    // The last knot has no span starting at it, so the search clamps it into the
    // final span. Past the last knot the time clamp does the rest, and a sample
    // there must be the SAME evaluation -- identical bits, not merely a nearby
    // value, because it is the identical span at the identical local time.
    auto pt = spline.evaluate(times.back());
    auto past = spline.evaluate(times.back() + 1.0);
    REQUIRE(past.position(0) == pt.position(0));
    REQUIRE(past.velocity(0) == pt.velocity(0));
    REQUIRE(past.acceleration(0) == pt.acceleration(0));

    constexpr double eps = std::numeric_limits<double>::epsilon();
    double const scale = *std::max_element(positions.begin(), positions.end());
    double const tol = static_cast<double>(spline_horner_rounding_ops) * eps * scale;
    CAPTURE(pt.position(0), positions.back(), tol);
    REQUIRE(std::abs(pt.position(0) - positions.back()) <= tol);
}

TEST_CASE("Cubic spline periodic BC on the smallest cyclic system this case builds",
          "[cubic_spline][hardening][coverage]")
{
    // Four waypoints, three spans, and a cyclic system of three unknowns. The
    // Sherman-Morrison reduction the periodic solve uses is exact down to two, so
    // this is not the smallest system it admits; it is the smallest one this file
    // exercises, and the case name used to claim three waypoints while the
    // configuration below has always had four.
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, -1.0, 0.0};

    ctrlpp::cubic_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .bc = ctrlpp::boundary_condition::periodic,
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);

    // As above: the wrap is imposed by the cyclic system, so the residual is the
    // rounding of the coefficients that carried it, at the scale of the largest
    // acceleration the spline reports. One budget, both quantities.
    constexpr double eps = std::numeric_limits<double>::epsilon();
    auto const spans = times.size() - 1;
    double const acceleration_scale = knot_acceleration_scale(spline, times);
    double const tol = static_cast<double>(cyclic_spline_acceleration_ops(spans))
                       * eps * acceleration_scale;

    auto start = spline.evaluate(0.0);
    auto end = spline.evaluate(3.0);
    CAPTURE(acceleration_scale, tol, start.velocity(0), end.velocity(0));
    REQUIRE(std::abs(start.velocity(0) - end.velocity(0)) <= tol);
    REQUIRE(std::abs(start.acceleration(0) - end.acceleration(0)) <= tol);
}

// ── Synchronize: empty vector and additional edge cases ───────────────────────

TEST_CASE("Synchronize empty vector is safe no-op",
          "[synchronize][hardening][coverage]")
{
    std::vector<ctrlpp::trapezoidal_trajectory<double>> empty;
    REQUIRE(ctrlpp::synchronize(std::span{empty}).has_value());
    REQUIRE(empty.empty());
}

TEST_CASE("Synchronize vector with single element is no-op",
          "[synchronize][hardening][coverage]")
{
    std::vector<ctrlpp::trapezoidal_trajectory<double>> axes;
    axes.push_back(trapezoidal_profile({.q0 = 0.0, .q1 = 5.0, .v_max = 2.0, .a_max = 1.0}));
    auto const dur_before = axes[0].duration();

    REQUIRE(ctrlpp::synchronize(std::span{axes}).has_value());

    // A single axis is passed a bit-exact copy of its own duration, so it takes
    // the success no-op path and nothing changes at all.
    REQUIRE(axes[0].duration() == dur_before);
}

// ── Trapezoidal: non-zero BCs and rescale edge cases ──────────────────────────

TEST_CASE("Trapezoidal reports the acceleration it raised the command to",
          "[trapezoidal][hardening][coverage]")
{
    // v0 and v1 are high relative to the displacement, so the two boundary
    // velocities cannot be reconciled over it at the commanded limit and the
    // construction raises the limit rather than refusing the command.
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.1, .v_max = 5.0, .a_max = 1.0,
        .v0 = 3.0, .v1 = 2.0,
    };

    auto const profile = trapezoidal_profile(cfg);
    auto const& disp = profile.disposition();

    REQUIRE(disp.commanded_acceleration == cfg.a_max);
    REQUIRE(disp.realized_acceleration > disp.commanded_acceleration);

    // The remedy's own closed form (B&M eq. (3.15)): the smallest acceleration
    // at which the two ramps cover the commanded displacement exactly is the
    // half-difference of the squared boundary velocities divided by that
    // displacement, plus the unit in the last place that keeps the raised value
    // on the feasible side of the test it was derived from. Evaluated here in
    // the same operation order the library uses, on the same operands, so exact
    // equality is the contract rather than a tolerance.
    double const abs_h = std::abs(cfg.q1 - cfg.q0);
    double const v_diff_sq = std::abs(cfg.v0 * cfg.v0 - cfg.v1 * cfg.v1) / 2.0;
    double const expected_a = v_diff_sq / abs_h + std::numeric_limits<double>::epsilon();
    CAPTURE(disp.realized_acceleration, expected_a);
    REQUIRE(disp.realized_acceleration == expected_a);

    // The profile is correct and limit-respecting under the REALIZED limit,
    // which is the whole reason the raise is a disposition and not a failure.
    // The commanded limit is not the scale these contracts are measured at --
    // it is not the limit the ramps run at.
    require_nonnegative_phases(profile);
    require_swept_displacement(profile, cfg.q1 - cfg.q0, disp.realized_acceleration);
    require_terminal_velocity(profile, cfg.v1, disp.realized_acceleration);
}

TEST_CASE("Trapezoidal reports an unraised acceleration as equal to the commanded one",
          "[trapezoidal][hardening][coverage]")
{
    // The same shape with room to spare: the displacement is large enough that
    // the boundary velocities are feasible at the commanded limit, so nothing is
    // raised. Without this case the assertion above cannot tell a disposition
    // that is always set from one that is set correctly.
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 20.0, .v_max = 5.0, .a_max = 1.0,
        .v0 = 3.0, .v1 = 2.0,
    };

    auto const profile = trapezoidal_profile(cfg);
    auto const& disp = profile.disposition();

    REQUIRE(disp.commanded_acceleration == cfg.a_max);
    REQUIRE(disp.realized_acceleration == cfg.a_max);

    require_nonnegative_phases(profile);
    require_swept_displacement(profile, cfg.q1 - cfg.q0, disp.realized_acceleration);
    require_terminal_velocity(profile, cfg.v1, disp.realized_acceleration);
}

TEST_CASE("Trapezoidal rescale_to very long duration",
          "[trapezoidal][hardening][coverage]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 5.0, .v_max = 2.0, .a_max = 1.0,
    };
    auto profile = trapezoidal_profile(cfg);

    // Rescale to a very long duration (100x original)
    auto const original_T = profile.duration();
    auto const rescaled = profile.rescale_to(original_T * 100.0);
    REQUIRE(rescaled.has_value());

    REQUIRE(profile.duration() > original_T * 10.0);
    require_realized_duration(profile, original_T * 100.0);

    // The start position is not a by-construction value, so it is still worth
    // asserting; the traversal is asserted by quadrature rather than by the end
    // position, which the final segment reproduces by construction. At time zero
    // the acceleration branch adds nothing to the commanded start, so the sample
    // is that start exactly.
    auto start = profile.evaluate(0.0);
    REQUIRE(start.position(0) == cfg.q0);
    REQUIRE(start.velocity(0) == cfg.v0);
    require_swept_displacement(profile, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(profile, cfg.v1, cfg.a_max);
    require_nonnegative_phases(profile);
}

TEST_CASE("Trapezoidal rescale_to rejects a zero-distance trajectory",
          "[trapezoidal][hardening][negative]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 3.0, .q1 = 3.0, .v_max = 2.0, .a_max = 1.0,
    };
    auto profile = trapezoidal_profile(cfg);

    // A stationary profile reaches its own duration and nothing longer: the
    // reachable maximum derived from the vanishing-cruise limit is zero at rest,
    // so every longer request is a typed rejection. No epsilon takes part in that.
    REQUIRE(profile.duration() == 0.0);

    auto const rescaled = profile.rescale_to(1.0);
    REQUIRE(!rescaled.has_value());
    REQUIRE(rescaled.error() == ctrlpp::trajectory_error::unreachable_duration);
    REQUIRE(profile.duration() == 0.0);
}

TEST_CASE("Trapezoidal negative direction with non-zero BCs",
          "[trapezoidal][hardening][coverage]")
{
    // Negative displacement with initial/final velocities
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 10.0, .q1 = 2.0, .v_max = 3.0, .a_max = 2.0,
        .v0 = -1.0, .v1 = -0.5,
    };

    auto profile = trapezoidal_profile(cfg);

    // Nothing was raised here, so the contracts are measured at the commanded
    // limit; the disposition says which.
    REQUIRE(profile.disposition().realized_acceleration == cfg.a_max);
    require_nonnegative_phases(profile);

    // A negative traversal with nonzero boundary velocities is what the swept
    // displacement and terminal velocity contracts were written for: the signed
    // displacement goes in as it stands and the commanded final velocity with it.
    // The end sample they replace returns the commanded value by construction and
    // was carrying half a percent of slack to say nothing.
    require_swept_displacement(profile, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(profile, cfg.v1, cfg.a_max);

    auto start = profile.evaluate(0.0);
    REQUIRE(start.position(0) == cfg.q0);
    REQUIRE(start.velocity(0) == cfg.v0);
}

// ── Smoothing spline: 2-point linear degeneration ─────────────────────────────

TEST_CASE("Smoothing spline with 2 points and mu near zero is still linear",
          "[smoothing_spline][hardening][coverage]")
{
    ctrlpp::smoothing_spline<double>::config cfg{
        .times = {0.0, 2.0},
        .positions = {1.0, 5.0},
        .mu = 0.01,
    };

    auto spline = realizable<ctrlpp::smoothing_spline<double>>(cfg);

    // With two waypoints the construction never reaches the smoothing solve: it
    // takes the linear path, where the tradeoff parameter plays no part and the
    // quadratic and cubic coefficients are assigned zero. All three quantities are
    // exact, and the same three the case's tolerances were describing to a percent.
    auto mid = spline.evaluate(1.0);
    REQUIRE(mid.position(0) == 3.0);
    REQUIRE(mid.velocity(0) == 2.0);
    REQUIRE(mid.acceleration(0) == 0.0);
}

TEST_CASE("Smoothing spline at the bottom of the tradeoff domain saturates on the line",
          "[smoothing_spline][hardening][coverage]")
{
    // There is no clamp here to exercise. The construction rejects a tradeoff
    // parameter outside (0, 1] and passes everything inside it straight into the
    // weight 2 (1 - mu) / (3 mu); nothing rounds the parameter up to anything.
    // What a parameter this small does is saturate: the weight is large enough
    // that the interior second derivatives have been driven to the bottom of the
    // significand and the smoothed positions ARE the least-squares line, so the
    // spline evaluates to that line rather than merely near it.
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, 0.5, 2.0};

    ctrlpp::smoothing_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .mu = 1e-15,
    };

    auto spline = realizable<ctrlpp::smoothing_spline<double>>(cfg);
    require_moment_relation(spline, times);

    constexpr double eps = std::numeric_limits<double>::epsilon();
    double const probe = 1.5;
    double const line = least_squares_line_at(times, positions, probe);
    double const scale = *std::max_element(positions.begin(), positions.end());
    double const tol = static_cast<double>(spline_horner_rounding_ops) * eps * scale;

    auto pt = spline.evaluate(probe);
    CAPTURE(pt.position(0), line, tol);
    REQUIRE(std::abs(pt.position(0) - line) <= tol);

    // Saturated means further descent changes nothing: a parameter an order
    // smaller lands on the same line to the same budget. That is the statement a
    // clamp would have made, made against the behavior that is actually there.
    auto deeper = realizable<ctrlpp::smoothing_spline<double>>(
        ctrlpp::smoothing_spline<double>::config{
            .times = times, .positions = positions, .mu = 1e-16});
    CAPTURE(deeper.evaluate(probe).position(0));
    REQUIRE(std::abs(deeper.evaluate(probe).position(0) - line) <= tol);
}

// ── Online planner 2nd: rest-to-rest subroutine ───────────────────────────────

TEST_CASE("Online planner 2nd rest-to-rest with near-zero displacement after brake",
          "[online_planner_2nd][hardening][coverage]")
{
    constexpr double step = 0.01;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build velocity, then retarget to a position very close to where we will
    // stop after braking -- exercises rest-to-rest with near-zero abs_h
    REQUIRE(planner.update(5.0).has_value());
    double t = 0.0;
    for (int i = 0; i < 20; ++i) {
        t += step;
        planner.sample(t);
    }

    auto state = planner.sample(t);
    // Estimate where we would stop after braking: q + v^2/(2*a)
    double stop_pos = state.position(0)
                      + state.velocity(0) * std::abs(state.velocity(0))
                            / (2.0 * cfg.a_max);
    // Target exactly at the stop position. The remaining displacement is the
    // rounding of that estimate against the planner's own stopping distance, which
    // is what puts the rest-to-rest solve on its degenerate branch. The planner
    // can still carry the velocity through, so which branch it selects is its
    // report to make rather than this case's to assume; what the case pins is
    // that the motion respects both limits and ends where it was sent.
    REQUIRE(planner.update(stop_pos).has_value());
    REQUIRE(planner.diagnostics().commanded_target == stop_pos);
    REQUIRE(planner.diagnostics().initial_velocity == state.velocity(0));

    for (int i = 0; i < 500; ++i) {
        t += step;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
    }
    REQUIRE(planner.is_settled());

    auto const settled = planner.sample(t + step);
    REQUIRE(settled.position(0) == stop_pos);
    REQUIRE(settled.velocity(0) == 0.0);
}

TEST_CASE("Online planner 2nd cruise phase with initial velocity",
          "[online_planner_2nd][hardening][coverage]")
{
    // Large displacement so the planner enters cruise phase even with initial velocity
    constexpr double target = 100.0;
    constexpr double step = 0.01;
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 4.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build some velocity first
    REQUIRE(planner.update(target).has_value());
    double t = 0.0;
    for (int i = 0; i < 10; ++i) {
        t += step;
        planner.sample(t);
    }

    // Retarget with the same large displacement, now carrying a velocity: the
    // commanded shape exists, so the planner carries it through rather than
    // braking, and the profile has room for a cruise phase.
    REQUIRE(planner.update(target).has_value());
    REQUIRE(planner.diagnostics().disposition
            == ctrlpp::online_planner_disposition::commanded_profile);
    double const planned_duration = planner.diagnostics().planned_duration;
    REQUIRE(planned_duration > 0.0);

    // In cruise the evaluation assigns the acceleration zero and the velocity the
    // cruise velocity, which on this profile IS the velocity limit, signed. Both
    // are assignments, so the detection is exact. A band fitted around them can be
    // set by a sample that is merely near cruise and missed on one that is exactly
    // in it, which is the opposite of what a coverage flag is for.
    int cruise_samples = 0;
    int const steps = static_cast<int>(planned_duration / step) + 200;
    for (int i = 0; i < steps; ++i) {
        t += step;
        auto const pt = planner.sample(t);
        CAPTURE(i, t);
        require_envelope_2nd(pt, cfg.v_max, cfg.a_max);
        if (pt.acceleration(0) == 0.0 && std::abs(pt.velocity(0)) == cfg.v_max) {
            ++cruise_samples;
        }
    }
    REQUIRE(cruise_samples > 0);

    // The loop now runs past the end of the profile, which the old one did not:
    // its last sample sat mid-deceleration and a tolerance of one percent of the
    // target was what let that pass for a settled value.
    auto final_pt = planner.sample(t + step);
    REQUIRE(planner.is_settled());
    REQUIRE(final_pt.position(0) == target);
    REQUIRE(final_pt.velocity(0) == 0.0);
}
