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
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <span>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
#include <utility>
#include <algorithm>

using Catch::Matchers::WithinAbs;

// ── Oracle for time-rescaled velocity profiles ────────────────────────────────
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

}

// ── Cubic spline hardening ─────────────────────────────────────────────────────

TEST_CASE("Cubic spline with exactly 2 points", "[cubic_spline][hardening][negative]")
{
    ctrlpp::cubic_spline<double>::config cfg{
        .times = {0.0, 1.0},
        .positions = {0.0, 1.0},
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);
    auto pt = spline.evaluate(0.5);
    REQUIRE(std::isfinite(pt.position(0)));
    REQUIRE(std::isfinite(pt.velocity(0)));
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

    for (std::size_t i = 0; i < times.size(); ++i) {
        auto pt = spline.evaluate(times[i]);
        REQUIRE_THAT(pt.position(0), WithinAbs(positions[i], 1e-10));
    }
}

TEST_CASE("Cubic spline with huge span", "[cubic_spline][hardening][negative]")
{
    ctrlpp::cubic_spline<double>::config cfg{
        .times = {0.0, 1e15, 2e15},
        .positions = {0.0, 1.0, 0.0},
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);
    auto pt = spline.evaluate(0.5e15);
    REQUIRE(std::isfinite(pt.position(0)));
}

// ── Smoothing spline hardening ─────────────────────────────────────────────────

TEST_CASE("Smoothing spline with mu=1 approaches interpolation", "[smoothing_spline][hardening][negative]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, 0.5, 2.0};

    ctrlpp::smoothing_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .mu = 1.0,
    };

    auto spline = realizable<ctrlpp::smoothing_spline<double>>(cfg);

    // At mu=1, should pass through all waypoints
    for (std::size_t i = 0; i < times.size(); ++i) {
        auto pt = spline.evaluate(times[i]);
        REQUIRE_THAT(pt.position(0), WithinAbs(positions[i], 0.01));
    }
}

TEST_CASE("Smoothing spline with very large lambda", "[smoothing_spline][hardening][negative]")
{
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, 0.5, 2.0};

    // mu near zero -> lambda very large -> maximum smoothness (near straight line)
    ctrlpp::smoothing_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .mu = 0.001,
    };

    auto spline = realizable<ctrlpp::smoothing_spline<double>>(cfg);
    auto pt = spline.evaluate(1.5);
    REQUIRE(std::isfinite(pt.position(0)));
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

TEST_CASE("Trapezoidal with zero distance", "[trapezoidal][hardening][negative]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 5.0, .q1 = 5.0, .v_max = 1.0, .a_max = 1.0,
    };

    auto const traj = trapezoidal_profile(cfg);
    REQUIRE_THAT(traj.duration(), WithinAbs(0.0, 1e-12));
    auto pt = traj.evaluate(0.0);
    REQUIRE_THAT(pt.position(0), WithinAbs(5.0, 1e-12));
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

    auto const traj = trapezoidal_profile(cfg);
    REQUIRE(traj.is_triangular());

    // Peak velocity for triangle: sqrt(a * h) = sqrt(2 * 0.5) = 1.0
    REQUIRE_THAT(traj.peak_velocity(), WithinAbs(1.0, 0.01));

    // Final position should be q1
    auto pt = traj.evaluate(traj.duration());
    REQUIRE_THAT(pt.position(0), WithinAbs(0.5, 1e-6));
}

// ── Double-S trajectory hardening ──────────────────────────────────────────────

TEST_CASE("Double-S with zero distance", "[double_s][hardening][negative]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 3.0, .q1 = 3.0, .v_max = 1.0, .a_max = 1.0, .j_max = 1.0,
    };

    auto const traj = double_s_profile(cfg);
    REQUIRE_THAT(traj.duration(), WithinAbs(0.0, 1e-12));
    auto pt = traj.evaluate(0.0);
    REQUIRE_THAT(pt.position(0), WithinAbs(3.0, 1e-12));
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

TEST_CASE("Online planner 2nd with instant target flip", "[online_planner_2nd][hardening][negative]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    planner.update(5.0);
    planner.sample(0.1);
    planner.sample(0.2);

    // Flip target mid-motion
    planner.update(-5.0);
    auto pt = planner.sample(0.3);
    REQUIRE(std::isfinite(pt.position(0)));
    REQUIRE(std::isfinite(pt.velocity(0)));
}

TEST_CASE("Online planner 2nd reaches target", "[online_planner_2nd][hardening][convergence]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 1.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    planner.update(3.0);

    double t = 0.0;
    for (int i = 0; i < 1000; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    auto pt = planner.sample(t + 0.01);
    REQUIRE_THAT(pt.position(0), WithinAbs(3.0, 1e-6));
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

TEST_CASE("Online planner 3rd with instant target reversal", "[online_planner_3rd][hardening][negative]")
{
    ctrlpp::online_planner_3rd<double>::config cfg{.v_max = 1.0, .a_max = 2.0, .j_max = 5.0};
    auto planner = realizable<ctrlpp::online_planner_3rd<double>>(cfg);

    planner.update(5.0);
    planner.sample(0.1);

    // Reverse direction
    planner.update(-5.0);
    auto pt = planner.sample(0.2);
    REQUIRE(std::isfinite(pt.position(0)));
}

TEST_CASE("Online planner 3rd reaches target", "[online_planner_3rd][hardening][convergence]")
{
    ctrlpp::online_planner_3rd<double>::config cfg{.v_max = 2.0, .a_max = 1.0, .j_max = 5.0};
    auto planner = realizable<ctrlpp::online_planner_3rd<double>>(cfg);

    planner.update(3.0);

    double t = 0.0;
    for (int i = 0; i < 2000; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    auto pt = planner.sample(t + 0.01);
    REQUIRE_THAT(pt.position(0), WithinAbs(3.0, 1e-4));
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
    auto spline = realizable<ctrlpp::smoothing_spline<double>>(cfg);

    auto pt = spline.evaluate(0.5);
    REQUIRE_THAT(pt.position(0), WithinAbs(2.5, 0.01));
    REQUIRE_THAT(spline.duration(), WithinAbs(1.0, 1e-12));
}

TEST_CASE("Trapezoidal trajectory rescale_to extends motion",
          "[trapezoidal][hardening][coverage]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 2.0,
    };
    auto traj = trapezoidal_profile(cfg);
    auto const original_T = traj.duration();

    auto const rescaled = traj.rescale_to(original_T * 2.0);
    REQUIRE(rescaled.has_value());
    REQUIRE(traj.duration() > original_T);
    require_realized_duration(traj, original_T * 2.0);

    // The traversed displacement is asserted by integrating the reported velocity
    // over the profile's own phase segments. A position sample at or near the end
    // is NOT a valid check: the final segment is written as an offset backwards
    // from the commanded displacement, so it returns the target position even on a
    // profile whose velocity integrates to something else entirely.
    require_swept_displacement(traj, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(traj, cfg.v1, cfg.a_max);
    require_nonnegative_phases(traj);
    REQUIRE(std::abs(traj.peak_velocity()) <= cfg.v_max);
}

TEST_CASE("Trapezoidal trajectory rescale_to shorter than current is rejected",
          "[trapezoidal][hardening][negative]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 2.0,
    };
    auto traj = trapezoidal_profile(cfg);
    auto const original_T = traj.duration();
    auto const original_phases = traj.phase_durations();

    auto const rescaled = traj.rescale_to(original_T * 0.5);
    REQUIRE(!rescaled.has_value());
    REQUIRE(rescaled.error() == ctrlpp::trajectory_error::duration_shorter_than_current);

    // Nothing may have moved: a rejected request leaves the profile untouched.
    REQUIRE(traj.duration() == original_T);
    REQUIRE(traj.phase_durations() == original_phases);
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
    auto traj = trapezoidal_profile(cfg);

    auto const rescaled = traj.rescale_to(10.0);
    REQUIRE(rescaled.has_value());

    auto const phases = traj.phase_durations();
    REQUIRE(phases[0] > 0.0);
    REQUIRE(phases[1] > 0.0);
    REQUIRE(phases[2] > 0.0);
    REQUIRE(std::abs(traj.peak_velocity()) < cfg.v0);

    require_realized_duration(traj, 10.0);
    require_swept_displacement(traj, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(traj, cfg.v1, cfg.a_max);
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
    auto traj = double_s_profile(cfg);
    auto const original_T = traj.duration();

    auto const rescaled = traj.rescale_to(original_T * 3.0);
    REQUIRE(rescaled.has_value());
    require_realized_duration(traj, original_T * 3.0);
    require_swept_displacement(traj, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(traj, cfg.v1, cfg.a_max);
    require_nonnegative_phases(traj);

    // The rebuilt profile respects the scaled limits it was built under.
    REQUIRE(std::abs(traj.peak_velocity()) <= cfg.v_max);
}

TEST_CASE("Double-S trajectory rescale_to rebuilds with nonzero boundary velocities",
          "[double_s][hardening][coverage]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0,
        .v0 = 1.0, .v1 = 0.5,
    };
    auto traj = double_s_profile(cfg);
    auto const original_T = traj.duration();

    auto const rescaled = traj.rescale_to(original_T * 1.5);
    REQUIRE(rescaled.has_value());
    require_realized_duration(traj, original_T * 1.5);
    require_swept_displacement(traj, cfg.q1 - cfg.q0, cfg.a_max);
    require_nonnegative_phases(traj);

    // The boundary velocities are the command and are left unscaled, so the
    // profile still arrives at the one it was given.
    require_terminal_velocity(traj, cfg.v1, cfg.a_max);
}

TEST_CASE("Double-S trajectory rescale_to rejects shortening and unreachable requests",
          "[double_s][hardening][negative]")
{
    ctrlpp::double_s_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0,
        .v0 = 1.0, .v1 = 0.5,
    };
    auto traj = double_s_profile(cfg);
    auto const original_T = traj.duration();
    auto const original_phases = traj.phase_durations();

    auto const shorter = traj.rescale_to(original_T * 0.5);
    REQUIRE(!shorter.has_value());
    REQUIRE(shorter.error() == ctrlpp::trajectory_error::duration_shorter_than_current);
    REQUIRE(traj.duration() == original_T);
    REQUIRE(traj.phase_durations() == original_phases);

    // Slowing the profile down means lowering the velocity limit, and the limit
    // cannot fall below the boundary velocities the caller commanded. That pins a
    // finite reachable maximum well short of ten times the current duration.
    auto const unreachable = traj.rescale_to(original_T * 10.0);
    REQUIRE(!unreachable.has_value());
    REQUIRE(unreachable.error() == ctrlpp::trajectory_error::unreachable_duration);
    REQUIRE(traj.duration() == original_T);
    REQUIRE(traj.phase_durations() == original_phases);

    // Its own duration is a success no-op on both profiles.
    REQUIRE(traj.rescale_to(original_T).has_value());
    REQUIRE(traj.duration() == original_T);
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
    auto traj = trapezoidal_profile(cfg);
    auto const original_T = traj.duration();

    auto const h = std::abs(cfg.q1 - cfg.q0);
    auto const v_min = std::sqrt((cfg.v0 * cfg.v0 + cfg.v1 * cfg.v1) / 2.0 - cfg.a_max * h);
    auto const T_max = (cfg.v0 + cfg.v1 - 2.0 * v_min) / cfg.a_max;
    REQUIRE(T_max > original_T);

    // Inside the reachable set the request is served.
    auto const inside = traj.rescale_to(0.5 * (original_T + T_max));
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
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 4.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Start moving to 10
    planner.update(10.0);
    // Sample partway through to build up velocity
    for (int i = 0; i < 20; ++i) {
        planner.sample(0.05 * static_cast<double>(i + 1));
    }
    // Retarget to opposite direction -- triggers braking
    planner.update(-5.0);
    auto pt = planner.sample(0.05 * 21);
    REQUIRE(std::isfinite(pt.position(0)));
    REQUIRE(std::isfinite(pt.velocity(0)));
}

TEST_CASE("Online planner 2nd with same position target is near-zero motion",
          "[online_planner_2nd][hardening][coverage]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Target at current position (0)
    planner.update(0.0);
    auto pt = planner.sample(0.01);
    REQUIRE_THAT(pt.position(0), WithinAbs(0.0, 1e-10));
}

TEST_CASE("Trapezoidal negative cruise duration clamped to zero",
          "[trapezoidal][hardening][coverage]")
{
    // Very short distance relative to max velocity -> T_v < 0 path
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 0.1, .v_max = 100.0, .a_max = 1.0,
    };
    auto traj = trapezoidal_profile(cfg);

    // Should be triangular (no cruise phase)
    REQUIRE(traj.is_triangular());
    auto pt = traj.evaluate(traj.duration());
    REQUIRE_THAT(pt.position(0), WithinAbs(0.1, 0.01));
}

// ── Online planner 2nd: overshoot and braking coverage ────────────────────────

TEST_CASE("Online planner 2nd overshoot recovery brakes and reverses",
          "[online_planner_2nd][hardening][coverage]")
{
    // Start moving AWAY from target: positive velocity, negative target
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 4.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // First, build up positive velocity toward +10
    planner.update(10.0);
    double t = 0.0;
    for (int i = 0; i < 30; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    // Now retarget behind us -- current velocity is positive, target is negative
    // This triggers wrong_direction detection and braking
    planner.update(-3.0);

    // Sample through the braking phase
    for (int i = 0; i < 200; ++i) {
        t += 0.01;
        auto pt = planner.sample(t);
        REQUIRE(std::isfinite(pt.position(0)));
        REQUIRE(std::isfinite(pt.velocity(0)));
        REQUIRE(std::isfinite(pt.acceleration(0)));
    }

    // Eventually should reach the target
    for (int i = 0; i < 800; ++i) {
        t += 0.01;
        planner.sample(t);
    }
    auto final_pt = planner.sample(t + 0.01);
    REQUIRE_THAT(final_pt.position(0), WithinAbs(-3.0, 1e-4));
    REQUIRE(planner.is_settled());
}

TEST_CASE("Online planner 2nd wrong-direction: positive velocity, target behind",
          "[online_planner_2nd][hardening][coverage]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 3.0, .a_max = 5.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build up positive velocity by targeting +5
    planner.update(5.0);
    double t = 0.0;
    for (int i = 0; i < 20; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    // Verify we have positive velocity
    auto mid = planner.sample(t);
    REQUIRE(mid.velocity(0) > 0.1);

    // Now target is behind current position (same sign but smaller)
    // This will trigger overshoot detection since stopping distance > displacement
    planner.update(mid.position(0) - 0.001);

    // Sample through -- must brake, stop, and reverse slightly
    for (int i = 0; i < 500; ++i) {
        t += 0.01;
        auto pt = planner.sample(t);
        REQUIRE(std::isfinite(pt.position(0)));
    }
    auto settled = planner.sample(t + 0.01);
    REQUIRE_THAT(settled.velocity(0), WithinAbs(0.0, 1e-6));
}

TEST_CASE("Online planner 2nd near-zero displacement with velocity triggers braking",
          "[online_planner_2nd][hardening][coverage]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build up velocity
    planner.update(5.0);
    double t = 0.0;
    for (int i = 0; i < 50; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    auto current = planner.sample(t);
    // Retarget to exactly the current position -- displacement is near zero but
    // velocity is nonzero, so the wrong_direction/overshoot check fires
    planner.update(current.position(0));

    for (int i = 0; i < 500; ++i) {
        t += 0.01;
        auto pt = planner.sample(t);
        REQUIRE(std::isfinite(pt.position(0)));
    }
    REQUIRE(planner.is_settled());
}

TEST_CASE("Online planner 2nd evaluate_profile at and past T boundary",
          "[online_planner_2nd][hardening][coverage]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    planner.update(1.0);

    // Sample well past when we should have arrived
    double t = 0.0;
    for (int i = 0; i < 2000; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    // At this point the planner is settled. Further samples past T should
    // return target position with zero velocity and acceleration.
    auto pt = planner.sample(t + 100.0);
    REQUIRE_THAT(pt.position(0), WithinAbs(1.0, 1e-8));
    REQUIRE_THAT(pt.velocity(0), WithinAbs(0.0, 1e-8));
    REQUIRE_THAT(pt.acceleration(0), WithinAbs(0.0, 1e-8));
}

TEST_CASE("Online planner 2nd braking phase evaluation covers all branches",
          "[online_planner_2nd][hardening][coverage]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 3.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build velocity in negative direction
    planner.update(-8.0);
    double t = 0.0;
    for (int i = 0; i < 40; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    auto state = planner.sample(t);
    REQUIRE(state.velocity(0) < -0.1);

    // Now target in positive direction -- triggers wrong_direction braking
    planner.update(5.0);

    // Sample finely through the braking phase to cover dt < T_brake_ branch
    for (int i = 0; i < 10; ++i) {
        t += 0.005;
        auto pt = planner.sample(t);
        REQUIRE(std::isfinite(pt.position(0)));
        REQUIRE(std::isfinite(pt.velocity(0)));
        REQUIRE(std::isfinite(pt.acceleration(0)));
    }

    // Continue through rest-to-rest phase after braking
    for (int i = 0; i < 1000; ++i) {
        t += 0.01;
        planner.sample(t);
    }
    auto final_pt = planner.sample(t + 0.01);
    REQUIRE_THAT(final_pt.position(0), WithinAbs(5.0, 1e-4));
}

TEST_CASE("Online planner 2nd reset clears state",
          "[online_planner_2nd][hardening][coverage]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 3.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    planner.update(10.0);
    planner.sample(0.5);

    // Reset to a new position
    planner.reset(7.0);
    REQUIRE(planner.is_settled());

    auto pt = planner.sample(0.0);
    REQUIRE_THAT(pt.position(0), WithinAbs(7.0, 1e-12));
    REQUIRE_THAT(pt.velocity(0), WithinAbs(0.0, 1e-12));
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

    // Velocity at t=0 should match velocity at t=4 (periodic wrap)
    auto start = spline.evaluate(times.front());
    auto end = spline.evaluate(times.back());
    REQUIRE_THAT(start.velocity(0), WithinAbs(end.velocity(0), 1e-8));

    // Acceleration at t=0 should match acceleration at t=4
    REQUIRE_THAT(start.acceleration(0), WithinAbs(end.acceleration(0), 1e-8));

    // Position at endpoints must match
    REQUIRE_THAT(start.position(0), WithinAbs(end.position(0), 1e-10));
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

    // Endpoints should match
    auto start = spline.evaluate(0.0);
    auto end = spline.evaluate(1.0);
    REQUIRE_THAT(start.position(0), WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(end.position(0), WithinAbs(1.0, 1e-10));

    // Velocities at endpoints should match the clamped values
    REQUIRE_THAT(start.velocity(0), WithinAbs(2.0, 1e-8));
    REQUIRE_THAT(end.velocity(0), WithinAbs(-1.0, 1e-8));
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

    // Endpoint velocities must match clamped values
    auto start = spline.evaluate(0.0);
    auto end = spline.evaluate(3.0);
    REQUIRE_THAT(start.velocity(0), WithinAbs(0.5, 1e-8));
    REQUIRE_THAT(end.velocity(0), WithinAbs(1.0, 1e-8));

    // Knot interpolation
    for (std::size_t i = 0; i < cfg.times.size(); ++i) {
        auto pt = spline.evaluate(cfg.times[i]);
        REQUIRE_THAT(pt.position(0), WithinAbs(cfg.positions[i], 1e-10));
    }
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

    // Evaluate exactly at each knot -- should not crash and return matching position
    for (std::size_t i = 0; i < times.size(); ++i) {
        auto pt = spline.evaluate(times[i]);
        REQUIRE_THAT(pt.position(0), WithinAbs(positions[i], 1e-10));
    }

    // Evaluate at the last knot time (edge case for find_span clamp)
    auto pt = spline.evaluate(times.back());
    REQUIRE_THAT(pt.position(0), WithinAbs(positions.back(), 1e-10));

    // Evaluate past the last knot (should be clamped)
    auto past = spline.evaluate(times.back() + 1.0);
    REQUIRE_THAT(past.position(0), WithinAbs(positions.back(), 1e-10));
}

TEST_CASE("Cubic spline periodic BC with 3 points (minimum for cyclic Thomas)",
          "[cubic_spline][hardening][coverage]")
{
    // Minimum periodic: 3 points, 2 spans, cyclic system size = 2
    // But cyclic_thomas_solve requires n >= 3. With n_pts=3, n=2 spans,
    // the periodic solver uses n=2 unknowns. Let us use 4 points instead.
    std::vector<double> times{0.0, 1.0, 2.0, 3.0};
    std::vector<double> positions{0.0, 1.0, -1.0, 0.0};

    ctrlpp::cubic_spline<double>::config cfg{
        .times = times,
        .positions = positions,
        .bc = ctrlpp::boundary_condition::periodic,
    };

    auto spline = realizable<ctrlpp::cubic_spline<double>>(cfg);

    auto start = spline.evaluate(0.0);
    auto end = spline.evaluate(3.0);
    REQUIRE_THAT(start.velocity(0), WithinAbs(end.velocity(0), 1e-8));
    REQUIRE_THAT(start.acceleration(0), WithinAbs(end.acceleration(0), 1e-8));
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

    auto const traj = trapezoidal_profile(cfg);
    auto const& disp = traj.disposition();

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
    require_nonnegative_phases(traj);
    require_swept_displacement(traj, cfg.q1 - cfg.q0, disp.realized_acceleration);
    require_terminal_velocity(traj, cfg.v1, disp.realized_acceleration);
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

    auto const traj = trapezoidal_profile(cfg);
    auto const& disp = traj.disposition();

    REQUIRE(disp.commanded_acceleration == cfg.a_max);
    REQUIRE(disp.realized_acceleration == cfg.a_max);

    require_nonnegative_phases(traj);
    require_swept_displacement(traj, cfg.q1 - cfg.q0, disp.realized_acceleration);
    require_terminal_velocity(traj, cfg.v1, disp.realized_acceleration);
}

TEST_CASE("Trapezoidal rescale_to very long duration",
          "[trapezoidal][hardening][coverage]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 0.0, .q1 = 5.0, .v_max = 2.0, .a_max = 1.0,
    };
    auto traj = trapezoidal_profile(cfg);

    // Rescale to a very long duration (100x original)
    auto const original_T = traj.duration();
    auto const rescaled = traj.rescale_to(original_T * 100.0);
    REQUIRE(rescaled.has_value());

    REQUIRE(traj.duration() > original_T * 10.0);
    require_realized_duration(traj, original_T * 100.0);

    // The start position is not a by-construction value, so it is still worth
    // asserting; the traversal is asserted by quadrature rather than by the end
    // position, which the final segment reproduces by construction.
    auto start = traj.evaluate(0.0);
    REQUIRE_THAT(start.position(0), WithinAbs(0.0, 1e-10));
    require_swept_displacement(traj, cfg.q1 - cfg.q0, cfg.a_max);
    require_terminal_velocity(traj, cfg.v1, cfg.a_max);
    require_nonnegative_phases(traj);
}

TEST_CASE("Trapezoidal rescale_to rejects a zero-distance trajectory",
          "[trapezoidal][hardening][negative]")
{
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 3.0, .q1 = 3.0, .v_max = 2.0, .a_max = 1.0,
    };
    auto traj = trapezoidal_profile(cfg);

    // A stationary profile reaches its own duration and nothing longer: the
    // reachable maximum derived from the vanishing-cruise limit is zero at rest,
    // so every longer request is a typed rejection. No epsilon takes part in that.
    REQUIRE(traj.duration() == 0.0);

    auto const rescaled = traj.rescale_to(1.0);
    REQUIRE(!rescaled.has_value());
    REQUIRE(rescaled.error() == ctrlpp::trajectory_error::unreachable_duration);
    REQUIRE(traj.duration() == 0.0);
}

TEST_CASE("Trapezoidal negative direction with non-zero BCs",
          "[trapezoidal][hardening][coverage]")
{
    // Negative displacement with initial/final velocities
    ctrlpp::trapezoidal_trajectory<double>::config cfg{
        .q0 = 10.0, .q1 = 2.0, .v_max = 3.0, .a_max = 2.0,
        .v0 = -1.0, .v1 = -0.5,
    };

    auto traj = trapezoidal_profile(cfg);
    REQUIRE(std::isfinite(traj.duration()));

    auto start = traj.evaluate(0.0);
    auto end = traj.evaluate(traj.duration());
    REQUIRE_THAT(start.position(0), WithinAbs(10.0, 1e-10));
    REQUIRE_THAT(end.position(0), WithinAbs(2.0, 0.01));
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

    // 2-point case always degenerates to linear regardless of mu
    auto mid = spline.evaluate(1.0);
    REQUIRE_THAT(mid.position(0), WithinAbs(3.0, 0.01));

    // Velocity should be constant (slope = 2.0)
    REQUIRE_THAT(mid.velocity(0), WithinAbs(2.0, 0.01));

    // Acceleration should be zero for linear
    REQUIRE_THAT(mid.acceleration(0), WithinAbs(0.0, 1e-10));
}

TEST_CASE("Smoothing spline with mu at machine epsilon clamp",
          "[smoothing_spline][hardening][coverage]")
{
    // mu very close to 0 triggers the clamp to eps
    ctrlpp::smoothing_spline<double>::config cfg{
        .times = {0.0, 1.0, 2.0, 3.0},
        .positions = {0.0, 1.0, 0.5, 2.0},
        .mu = 1e-15,
    };

    auto spline = realizable<ctrlpp::smoothing_spline<double>>(cfg);
    auto pt = spline.evaluate(1.5);
    REQUIRE(std::isfinite(pt.position(0)));
    REQUIRE(std::isfinite(pt.velocity(0)));
}

// ── Online planner 2nd: rest-to-rest subroutine ───────────────────────────────

TEST_CASE("Online planner 2nd rest-to-rest with near-zero displacement after brake",
          "[online_planner_2nd][hardening][coverage]")
{
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 1.0, .a_max = 2.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build velocity, then retarget to a position very close to where we will
    // stop after braking -- exercises rest-to-rest with near-zero abs_h
    planner.update(5.0);
    double t = 0.0;
    for (int i = 0; i < 20; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    auto state = planner.sample(t);
    // Estimate where we would stop after braking: q + v^2/(2*a)
    double stop_pos = state.position(0)
                      + state.velocity(0) * std::abs(state.velocity(0))
                            / (2.0 * cfg.a_max);
    // Target exactly at the stop position
    planner.update(stop_pos);

    for (int i = 0; i < 500; ++i) {
        t += 0.01;
        auto pt = planner.sample(t);
        REQUIRE(std::isfinite(pt.position(0)));
    }
    REQUIRE(planner.is_settled());
}

TEST_CASE("Online planner 2nd cruise phase with initial velocity",
          "[online_planner_2nd][hardening][coverage]")
{
    // Large displacement so the planner enters cruise phase even with initial velocity
    ctrlpp::online_planner_2nd<double>::config cfg{.v_max = 2.0, .a_max = 4.0};
    auto planner = realizable<ctrlpp::online_planner_2nd<double>>(cfg);

    // Build some velocity first
    planner.update(100.0);
    double t = 0.0;
    for (int i = 0; i < 10; ++i) {
        t += 0.01;
        planner.sample(t);
    }

    // Retarget with large displacement -- will use compute_with_initial_velocity
    // and should enter cruise phase (v_tri >= v_max)
    planner.update(100.0);

    bool saw_cruise = false;
    for (int i = 0; i < 5000; ++i) {
        t += 0.01;
        auto pt = planner.sample(t);
        REQUIRE(std::isfinite(pt.position(0)));
        // Cruise phase: velocity near v_max, acceleration near zero
        if (std::abs(pt.acceleration(0)) < 0.01
            && std::abs(std::abs(pt.velocity(0)) - cfg.v_max) < 0.1) {
            saw_cruise = true;
        }
    }
    REQUIRE(saw_cruise);

    auto final_pt = planner.sample(t + 0.01);
    REQUIRE_THAT(final_pt.position(0), WithinAbs(100.0, 1.0));
}
