// This anchor holds the double-S constructor to the contract that matters for
// a velocity profile: over its own reported duration it must sweep exactly the
// displacement it was commanded, and it must leave and arrive at the boundary
// velocities it was given.
//
// The oracle is deliberately NOT a position sample taken at or near an
// endpoint. The profile's final segment is written as an offset backwards from
// the commanded displacement, so it returns the target position by
// construction; a profile whose velocity integrates to a thousand times the
// commanded distance still lands its endpoint exactly. Position is therefore
// asserted only through numerical integration of the reported velocity, and
// that quadrature is aligned to the profile's own segment boundaries -- see the
// rounding-op budget below for why a single uniform grid over the whole
// duration is rejected. Two further checks stand beside it: the velocity a step
// inside each end against the commanded boundary velocity, and an independent
// two-sided bound on the reported duration that does not borrow its scale from
// the duration being judged.
//
// The cruise-free family is served exhaustively: either the construction finds
// the peak velocity that sweeps the commanded displacement, or the displacement
// is smaller than the fastest admissible transition between the two boundary
// velocities and the configuration is rejected as
// trajectory_error::unreachable_boundary_velocity. There is no third outcome in
// which a finite profile is returned that does not traverse its own command.

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <algorithm>

using namespace ctrlpp;

namespace
{

/// Seed for every randomized family below. Named and captured by each case so a
/// failure prints its own reproducer instead of requiring a hunt through the
/// generator.
constexpr std::uint32_t sweep_seed = 20260727;

// Rounding-op budget for the kink-aligned quadrature.
//
// Simpson's rule integrates every polynomial up to cubic order exactly, and the
// double-S velocity is at most quadratic inside each of its seven segments, so
// a quadrature whose panels are aligned to the segment boundaries carries no
// truncation error at all -- the entire residual is floating-point rounding.
// One uniform grid spanning the whole duration does NOT have that property: its
// panels straddle the segment kinks, where the integrand's derivative jumps,
// and the truncation error manufactured there is orders of magnitude above the
// rounding floor this budget describes. That is why the panels below are laid
// out inside each segment and never as a single grid over the whole move; do
// not "simplify" this back to a uniform grid.
//
// Three chains contribute, each worth up to one unit in the last place of the
// accumulated absolute area:
//   * every velocity sample: evaluate() forms a cubic in a local time and, in
//     the deceleration segments, a subtraction of that time from the end of the
//     move -- at most 8 chained rounding operations;
//   * every Simpson panel: three multiplies and three adds over its samples;
//   * the two running sums: one addition each per panel.
// Two fresh samples are taken per panel, so the budget per panel is
// 2 * 8 + 6 + 2 rounding operations.
constexpr int rounding_ops_per_panel = 2 * 8 + 6 + 2;

/// Panels laid inside each of the seven segments. Simpson is already exact on
/// the integrand, so the count buys no accuracy; it is kept small so the budget
/// above stays a tight description of the residual.
constexpr int panels_per_segment = 4;

/// Chained rounding operations behind a single reported quantity (a phase
/// duration sum, or one velocity sample compared against a boundary value).
constexpr int rounding_ops_per_sample = 8;

struct sweep_config
{
    double q0, q1, v_max, a_max, j_max, v0, v1;
};

/// Duration of the fastest admissible jerk-limited transition spanning a
/// velocity change. The constant-acceleration segment vanishes when the change
/// is too small to build up to the acceleration limit, and the ramp is then
/// triangular in acceleration.
auto ramp_duration(double dv, double a_max, double j_max) -> double
{
    double const d = std::abs(dv);
    if(d * j_max < a_max * a_max)
        return 2.0 * std::sqrt(d / j_max);
    return d / a_max + a_max / j_max;
}

/// Distance swept by the fastest admissible transition between the two
/// positive-frame boundary velocities. Each ramp's acceleration profile is
/// point-symmetric about its own midpoint, so the mean velocity it sweeps is
/// exactly the average of its endpoints. No profile in this family covers less
/// ground than that, which makes it the feasibility floor on the commanded
/// displacement.
auto minimum_displacement(double pv0, double pv1, double a_max, double j_max) -> double
{
    double const hi = std::max(pv0, pv1);
    double const lo = std::min(pv0, pv1);
    return 0.5 * (hi + lo) * ramp_duration(hi - lo, a_max, j_max);
}

struct quadrature
{
    double integral{};
    double abs_area{};
    int panels{};
};

template <typename Trajectory>
auto sample_velocity(Trajectory const& profile, double t) -> double
{
    using Scalar = typename Trajectory::scalar_type;
    return static_cast<double>(profile.evaluate(static_cast<Scalar>(t)).velocity(0));
}

/// Kink-aligned composite Simpson integration of the reported velocity over the
/// profile's own reported duration. The segment lengths come from the profile
/// itself, are turned into boundaries by prefix sum, and each segment is
/// integrated separately.
template <typename Trajectory>
auto integrate_velocity(Trajectory const& profile) -> quadrature
{
    auto const segments = profile.phase_durations();

    quadrature out{};
    double t_start = 0.0;
    for(auto const& segment : segments)
    {
        double const length = static_cast<double>(segment);
        if(!(length > 0.0))
            continue;

        double const dt = length / static_cast<double>(panels_per_segment);
        for(int p = 0; p < panels_per_segment; ++p)
        {
            double const a = t_start + static_cast<double>(p) * dt;
            double const b = a + dt;
            double const m = a + 0.5 * dt;
            double const panel = (dt / 6.0)
                * (sample_velocity(profile, a) + 4.0 * sample_velocity(profile, m) + sample_velocity(profile, b));
            out.integral += panel;
            out.abs_area += std::abs(panel);
            ++out.panels;
        }
        t_start += length;
    }
    return out;
}

/// Assert the constructor's fidelity contract for one configuration on one
/// scalar tier: swept displacement, boundary velocities, segment bookkeeping,
/// and an independent duration bound.
template <typename Scalar>
void check_constructor_fidelity(sweep_config const& cfg)
{
    constexpr double eps = static_cast<double>(std::numeric_limits<Scalar>::epsilon());

    double_s_trajectory<Scalar> profile({.q0 = static_cast<Scalar>(cfg.q0),
                                         .q1 = static_cast<Scalar>(cfg.q1),
                                         .v_max = static_cast<Scalar>(cfg.v_max),
                                         .a_max = static_cast<Scalar>(cfg.a_max),
                                         .j_max = static_cast<Scalar>(cfg.j_max),
                                         .v0 = static_cast<Scalar>(cfg.v0),
                                         .v1 = static_cast<Scalar>(cfg.v1)});

    double const T = static_cast<double>(profile.duration());
    double const h_signed = static_cast<double>(static_cast<Scalar>(cfg.q1)) - static_cast<double>(static_cast<Scalar>(cfg.q0));
    double const h = std::abs(h_signed);

    CAPTURE(T, h_signed);
    REQUIRE(std::isfinite(T));
    REQUIRE(T > 0.0);

    // Duration, bounded independently of the reported duration itself. A bound
    // scaled by T is vacuous exactly when T is the thing that is wrong.
    //
    // Lower: no admissible velocity exceeds the larger of the velocity limit and
    // the two boundary velocities, so the swept distance cannot exceed that
    // ceiling times the duration.
    //
    // Upper: the slowest admissible traversal cruises the whole displacement at
    // the velocity limit and spends one full ramp at each end. A single ramp
    // spans at most twice the velocity ceiling and lasts no longer than the
    // slower of its two admissible shapes, and a ramp whose mean velocity is
    // negative can only add its own swept distance back to the cruise segment.
    double const v_ceiling = std::max({cfg.v_max, std::abs(cfg.v0), std::abs(cfg.v1)});
    double const ramp_bound = std::max(2.0 * v_ceiling / cfg.a_max + cfg.a_max / cfg.j_max,
                                       2.0 * std::sqrt(2.0 * v_ceiling / cfg.j_max));
    double const T_min = h / v_ceiling;
    double const T_max = (h + 2.0 * v_ceiling * ramp_bound) / cfg.v_max + 2.0 * ramp_bound;
    double const T_bound_tol = static_cast<double>(rounding_ops_per_sample) * eps * T_max;
    CAPTURE(T_min, T_max);
    REQUIRE(T >= T_min - T_bound_tol);
    REQUIRE(T <= T_max + T_bound_tol);

    // Segment bookkeeping.
    auto const segments = profile.phase_durations();
    double segment_sum = 0.0;
    for(auto const& segment : segments)
    {
        double const length = static_cast<double>(segment);
        REQUIRE(std::isfinite(length));
        REQUIRE(length >= 0.0);
        segment_sum += length;
    }
    CAPTURE(segment_sum);
    REQUIRE(std::abs(segment_sum - T) <= static_cast<double>(rounding_ops_per_sample) * eps * T);

    // Swept displacement, by kink-aligned quadrature of the reported velocity.
    auto const q = integrate_velocity(profile);
    double const area_scale = std::max(q.abs_area, h);
    double const displacement_tol =
        static_cast<double>(rounding_ops_per_panel * q.panels) * eps * area_scale;
    CAPTURE(q.integral, q.abs_area, q.panels, displacement_tol);
    REQUIRE(std::abs(q.integral - h_signed) <= displacement_tol);

    // Boundary velocities. The acceleration is zero at both ends and the jerk
    // limit bounds how fast it can leave zero, so the velocity a step delta
    // inside either end differs from the boundary value by at most
    // j_max * delta^2 / 2. The step is the profile's own first (respectively
    // last) segment length, which is where that bound is attained exactly.
    //
    // Two chains set the margin. The sample itself carries the same 8 chained
    // rounding operations as any other evaluation, at the velocity scale. The
    // probe time is the more expensive one: a time measured backwards from the
    // end of the move is representable only to one unit in the last place of the
    // duration, and the profile recovers its own local time by subtracting that
    // probe time from the duration again, so the step it actually sees differs
    // from the step asked for at that scale. The velocity slews at up to the
    // acceleration limit there, which converts that time uncertainty into
    // a_max * eps * T of velocity. The term is not a duration-scaled tolerance
    // standing in for a duration check -- the duration is bounded independently
    // above, and this term stays far below the defect signal even when it is not.
    double const v_scale = std::max(v_ceiling, 1.0);
    double const v_tol = static_cast<double>(rounding_ops_per_sample) * eps * (v_scale + cfg.a_max * T);

    double const delta_start = static_cast<double>(segments.front());
    if(delta_start > 0.0)
    {
        double const v = sample_velocity(profile, delta_start);
        double const jerk_bound = 0.5 * cfg.j_max * delta_start * delta_start;
        CAPTURE(delta_start, v, jerk_bound);
        REQUIRE(std::abs(v - cfg.v0) <= jerk_bound + v_tol);
    }

    double const delta_end = static_cast<double>(segments.back());
    if(delta_end > 0.0)
    {
        double const v = sample_velocity(profile, T - delta_end);
        double const jerk_bound = 0.5 * cfg.j_max * delta_end * delta_end;
        CAPTURE(delta_end, v, jerk_bound);
        REQUIRE(std::abs(v - cfg.v1) <= jerk_bound + v_tol);
    }
}

enum class boundary_family
{
    zero,       ///< both boundary velocities zero: the family the constructor already served
    equal,      ///< equal nonzero boundary velocities
    differing,  ///< differing nonzero boundary velocities
    reversing,  ///< initial motion away from the target
    no_cruise,  ///< displacement just above the feasibility floor, so no cruise segment exists
    cruising,   ///< displacement far above the floor, so a cruise segment exists
    below_floor, ///< displacement below the feasibility floor: not realizable at all
};

constexpr int configs_per_family = 60;

/// Seeded configurations for one family, generated in the positive-displacement
/// frame and then folded onto a random direction of travel. The commanded
/// displacement is always drawn above the family's own feasibility floor, so
/// every configuration produced here is realizable; the unrealizable region has
/// its own cases below.
auto make_family_configs(boundary_family family) -> std::array<sweep_config, configs_per_family>
{
    std::mt19937 gen(sweep_seed + static_cast<std::uint32_t>(family));
    std::uniform_real_distribution<double> q0_dist(-20.0, 20.0);
    std::uniform_real_distribution<double> v_max_dist(1.0, 6.0);
    std::uniform_real_distribution<double> a_max_dist(0.5, 20.0);
    std::uniform_real_distribution<double> j_max_dist(1.0, 80.0);
    std::uniform_real_distribution<double> frac_dist(0.05, 0.9);
    std::uniform_real_distribution<double> reverse_frac_dist(-0.8, -0.05);
    std::uniform_real_distribution<double> decade_dist(-2.0, 2.0);
    std::uniform_real_distribution<double> slack_dist(1.3, 4.0);
    std::uniform_real_distribution<double> shortfall_dist(0.05, 0.8);
    std::uniform_real_distribution<double> cruise_dist(20.0, 400.0);
    std::bernoulli_distribution flip(0.5);

    std::array<sweep_config, configs_per_family> configs{};
    for(auto& cfg : configs)
    {
        double const v_max = v_max_dist(gen);
        double const a_max = a_max_dist(gen);
        double const j_max = j_max_dist(gen);

        double pv0 = 0.0;
        double pv1 = 0.0;
        switch(family)
        {
        case boundary_family::zero:
            break;
        case boundary_family::equal:
            pv0 = frac_dist(gen) * v_max;
            pv1 = pv0;
            break;
        case boundary_family::reversing:
            pv0 = reverse_frac_dist(gen) * v_max;
            pv1 = frac_dist(gen) * v_max;
            if(flip(gen))
                std::swap(pv0, pv1);
            break;
        default:
            pv0 = frac_dist(gen) * v_max;
            pv1 = frac_dist(gen) * v_max;
            break;
        }

        double const floor_h = minimum_displacement(pv0, pv1, a_max, j_max);
        double h = 0.0;
        switch(family)
        {
        case boundary_family::no_cruise:
            h = slack_dist(gen) * std::max(floor_h, v_max * v_max / (2.0 * a_max) * 0.01);
            break;
        case boundary_family::cruising:
            h = cruise_dist(gen) * v_max * v_max / a_max;
            break;
        case boundary_family::below_floor:
            h = shortfall_dist(gen) * floor_h;
            break;
        default:
            h = std::pow(10.0, decade_dist(gen));
            break;
        }
        if(family != boundary_family::below_floor)
            h = std::max(h, slack_dist(gen) * floor_h);

        double const direction = flip(gen) ? 1.0 : -1.0;
        double const q0 = q0_dist(gen);
        cfg = {q0, q0 + direction * h, v_max, a_max, j_max, direction * pv0, direction * pv1};
    }
    return configs;
}

/// Recorded configurations whose constructed profile did not traverse its
/// commanded displacement: the reported duration exceeded the admissible bound
/// by three orders of magnitude and the integrated velocity missed the command
/// by a comparable factor.
constexpr sweep_config recorded_forward_violation{
    .q0 = 5.911630,
    .q1 = 6.893902,
    .v_max = 5.905937,
    .a_max = 18.316431,
    .j_max = 66.404132,
    .v0 = 2.447245,
    .v1 = 0.193441,
};

constexpr sweep_config recorded_reversed_violation{
    .q0 = 1.744413,
    .q1 = -2.633139,
    .v_max = 5.808603,
    .a_max = 13.340846,
    .j_max = 3.971480,
    .v0 = -1.630792,
    .v1 = -3.772260,
};

}

TEST_CASE("double-S profiles sweep their commanded displacement at zero boundary velocity", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);
    for(auto const& cfg : make_family_configs(boundary_family::zero))
    {
        CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.j_max, cfg.v0, cfg.v1);
        check_constructor_fidelity<double>(cfg);
        check_constructor_fidelity<float>(cfg);
    }
}

TEST_CASE("double-S profiles sweep their commanded displacement at equal nonzero boundary velocities", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);
    for(auto const& cfg : make_family_configs(boundary_family::equal))
    {
        CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.j_max, cfg.v0, cfg.v1);
        check_constructor_fidelity<double>(cfg);
        check_constructor_fidelity<float>(cfg);
    }
}

TEST_CASE("double-S profiles sweep their commanded displacement at differing nonzero boundary velocities", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);
    for(auto const& cfg : make_family_configs(boundary_family::differing))
    {
        CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.j_max, cfg.v0, cfg.v1);
        check_constructor_fidelity<double>(cfg);
        check_constructor_fidelity<float>(cfg);
    }
}

TEST_CASE("double-S profiles sweep their commanded displacement when the motion starts away from the target", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);
    for(auto const& cfg : make_family_configs(boundary_family::reversing))
    {
        CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.j_max, cfg.v0, cfg.v1);
        check_constructor_fidelity<double>(cfg);
        check_constructor_fidelity<float>(cfg);
    }
}

TEST_CASE("double-S profiles sweep their commanded displacement with no cruise segment", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);
    for(auto const& cfg : make_family_configs(boundary_family::no_cruise))
    {
        CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.j_max, cfg.v0, cfg.v1);
        check_constructor_fidelity<double>(cfg);
        check_constructor_fidelity<float>(cfg);
    }
}

TEST_CASE("double-S profiles sweep their commanded displacement with a cruise segment", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);
    for(auto const& cfg : make_family_configs(boundary_family::cruising))
    {
        CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.j_max, cfg.v0, cfg.v1);
        check_constructor_fidelity<double>(cfg);
        check_constructor_fidelity<float>(cfg);
    }
}

TEST_CASE("double-S construction rejects displacements below the transition it already sweeps", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);
    for(auto const& cfg : make_family_configs(boundary_family::below_floor))
    {
        CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.j_max, cfg.v0, cfg.v1);

        auto const result = double_s_trajectory<double>::try_create({.q0 = cfg.q0,
                                                                    .q1 = cfg.q1,
                                                                    .v_max = cfg.v_max,
                                                                    .a_max = cfg.a_max,
                                                                    .j_max = cfg.j_max,
                                                                    .v0 = cfg.v0,
                                                                    .v1 = cfg.v1});
        REQUIRE(!result.has_value());
        REQUIRE(result.error() == trajectory_error::unreachable_boundary_velocity);

        // The non-fallible constructor holds the start position for a zero
        // duration instead of reporting a traversal it never performs.
        double_s_trajectory<double> profile({.q0 = cfg.q0,
                                             .q1 = cfg.q1,
                                             .v_max = cfg.v_max,
                                             .a_max = cfg.a_max,
                                             .j_max = cfg.j_max,
                                             .v0 = cfg.v0,
                                             .v1 = cfg.v1});
        CAPTURE(profile.duration());
        REQUIRE(profile.duration() == 0.0);
        REQUIRE(profile.is_degenerate());
    }
}

TEST_CASE("double-S construction rejects a displacement shorter than the boundary velocities allow", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);

    // Slowing from 3.0 to 0.5 already sweeps roughly 0.55 within these limits,
    // so a command of 0.1 cannot be realized without overshooting and returning.
    sweep_config const cfg{
        .q0 = 0.0, .q1 = 0.1, .v_max = 5.0, .a_max = 20.0, .j_max = 100.0, .v0 = 3.0, .v1 = 0.5};
    CAPTURE(minimum_displacement(cfg.v0, cfg.v1, cfg.a_max, cfg.j_max));

    auto const rejected = double_s_trajectory<double>::try_create({.q0 = cfg.q0,
                                                                   .q1 = cfg.q1,
                                                                   .v_max = cfg.v_max,
                                                                   .a_max = cfg.a_max,
                                                                   .j_max = cfg.j_max,
                                                                   .v0 = cfg.v0,
                                                                   .v1 = cfg.v1});
    REQUIRE(!rejected.has_value());
    REQUIRE(rejected.error() == trajectory_error::unreachable_boundary_velocity);

    // Widening the command past that floor is accepted and traverses it.
    sweep_config const widened{
        .q0 = 0.0, .q1 = 2.0, .v_max = 5.0, .a_max = 20.0, .j_max = 100.0, .v0 = 3.0, .v1 = 0.5};
    auto const accepted = double_s_trajectory<double>::try_create({.q0 = widened.q0,
                                                                   .q1 = widened.q1,
                                                                   .v_max = widened.v_max,
                                                                   .a_max = widened.a_max,
                                                                   .j_max = widened.j_max,
                                                                   .v0 = widened.v0,
                                                                   .v1 = widened.v1});
    REQUIRE(accepted.has_value());
    check_constructor_fidelity<double>(widened);
    check_constructor_fidelity<float>(widened);
}

TEST_CASE("double-S profile sweeps its commanded displacement at a recorded forward violation", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);
    check_constructor_fidelity<double>(recorded_forward_violation);
    check_constructor_fidelity<float>(recorded_forward_violation);
}

TEST_CASE("double-S profile sweeps its commanded displacement at a recorded reversed violation", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);
    check_constructor_fidelity<double>(recorded_reversed_violation);
    check_constructor_fidelity<float>(recorded_reversed_violation);
}
