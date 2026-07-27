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
//
// The second half of the file holds both profiles to the same contract after a
// time rescaling: the retimed profile must still sweep its commanded
// displacement, still leave and arrive at its commanded boundary velocities,
// still respect its kinematic limits, stay continuous through the derivative
// order it guarantees, and land on the requested duration within a budget
// derived from the conditioning of its own solve. A request that cannot be
// realized is a specific typed rejection, never a silently wrong profile.

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

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
    double back_time_exposure{};
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
    double widest_panel = 0.0;
    double last_length = 0.0;
    for(auto const& segment : segments)
    {
        double const length = static_cast<double>(segment);
        if(!(length > 0.0))
            continue;

        double const dt = length / static_cast<double>(panels_per_segment);
        widest_panel = std::max(widest_panel, dt);
        last_length = length;
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
    out.back_time_exposure = last_length + widest_panel;
    return out;
}

/// Assert the constructor's fidelity contract for one configuration on one
/// scalar tier: swept displacement, boundary velocities, segment bookkeeping,
/// and an independent duration bound.
template <typename Scalar>
void check_constructor_fidelity(sweep_config const& cfg)
{
    constexpr double eps = static_cast<double>(std::numeric_limits<Scalar>::epsilon());

    auto const created = double_s_trajectory<Scalar>::create({.q0 = static_cast<Scalar>(cfg.q0),
                                                              .q1 = static_cast<Scalar>(cfg.q1),
                                                              .v_max = static_cast<Scalar>(cfg.v_max),
                                                              .a_max = static_cast<Scalar>(cfg.a_max),
                                                              .j_max = static_cast<Scalar>(cfg.j_max),
                                                              .v0 = static_cast<Scalar>(cfg.v0),
                                                              .v1 = static_cast<Scalar>(cfg.v1)});
    REQUIRE(created.has_value());
    auto const& profile = created.value();

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

// ---------------------------------------------------------------------------
// Time-scaling contract
// ---------------------------------------------------------------------------

/// Multiples of a profile's own duration requested from it. A small multiple
/// keeps the retimed profile close to the shape it started with, a moderate one
/// moves it a shape or two along, and a large one drives the cruise velocity
/// down towards the vanishing-cruise limit where the solve is worst conditioned.
constexpr double duration_multiples[] = {1.05, 3.0, 40.0};

/// Sample count for the continuity and limit scan. The grid is laid over the
/// COMMANDED duration, never the realized one: a step derived from a duration
/// that is itself wrong widens the bound by exactly the error it is meant to
/// catch.
constexpr int continuity_samples = 64;

/// Chained rounding operations behind the trapezoidal cruise velocity and the
/// duration recomputed from it. Twelve form the root -- two for the linear
/// coefficient, three for the constant, three for the discriminant, one for its
/// square root, and three for the cancellation-free root selection -- and
/// fifteen form the duration -- two for each ramp duration, three for each ramp
/// distance, three for the cruise duration, and two for the three-term sum.
constexpr int duration_ops_trapezoidal = 27;

/// Chained rounding operations behind the double-S duration on the rest-to-rest
/// path, where the scale is a single quotient of the two durations: one for that
/// quotient, six for the three scaled limits (one, two, and three multiplies),
/// five for the phase durations the rebuild derives from them, and four for the
/// sum.
constexpr int duration_ops_double_s_closed_form = 16;

/// The same chain on the solved path, plus one: the bracket is halved until its
/// midpoint lands on an endpoint, so the accepted scale sits within one unit in
/// the last place of the crossing, and within a fixed segment shape the duration
/// is proportional to the reciprocal of that scale, so it inherits that unit at
/// its own magnitude.
constexpr int duration_ops_double_s_solved = duration_ops_double_s_closed_form + 1;

struct scaling_limits
{
    double v_max{};
    double a_max{};
    double j_max{}; ///< zero for a profile that does not bound jerk
};

/// Conditioning of the trapezoidal cruise-velocity solve.
///
/// The realized duration is recomputed from the cruise velocity the solve
/// returns, and the duration's sensitivity to that velocity is what turns the
/// solve's roundings into a duration error. Within a shape,
/// T(v) = (+-)(v - v0 - v1)/a + D/v with D the shape's residual displacement, so
/// |dT/dv| * |dv| is bounded by the relative error of v times
/// (v/a + |D|/v), and the second term is the duration itself less its ramp part.
/// The relative error of v is NOT a handful of epsilons: each shape's root
/// carries a difference of two nearly equal quantities, and the ratio of the
/// operands entering that difference to the difference itself is the
/// amplification.
///
///  * plateau and valley: the constant term is a h +- (v0^2 + v1^2)/2 and the
///    discriminant is b^2 - 4c; both can cancel, and the larger amplification of
///    the two governs.
///  * ramp-through: the residual is h minus the distance the two ramps already
///    sweep, taken against the scale of the operands that formed it.
///
/// Returning the amplification rather than folding a fixed factor in keeps the
/// bound tight where the problem is well conditioned, which is where a defect
/// would actually have to hide.
auto trapezoidal_solve_conditioning(double v_cruise, double a, double h, double v0, double v1,
                                    double T_target) -> double
{
    double const v_lo = std::min(v0, v1);
    double const v_hi = std::max(v0, v1);
    double const v_sum_sq = v0 * v0 + v1 * v1;

    auto const ratio = [](double operands, double residual) -> double {
        double const scale = std::abs(residual);
        return (scale > 0.0) ? std::abs(operands) / scale : 1.0;
    };

    if(v_cruise > v_lo && v_cruise < v_hi)
    {
        double const ramp_distance = (v_hi - v_lo) * (v_hi + v_lo) / (2.0 * a);
        return std::max(1.0, ratio(std::max(h, v_sum_sq / (2.0 * a)), h - ramp_distance));
    }

    bool const plateau = (v_cruise >= v_hi);
    double const b = plateau ? ((v0 + v1) + a * T_target) : (a * T_target - (v0 + v1));
    double const c = plateau ? (a * h + v_sum_sq / 2.0) : (v_sum_sq / 2.0 - a * h);
    double const coefficient = ratio(a * h + v_sum_sq / 2.0, c);
    double const discriminant = ratio(b * b + 4.0 * std::abs(c), b * b - 4.0 * c);
    return std::max({1.0, coefficient, discriminant});
}

/// Assert the whole time-scaling contract on a profile that has just accepted a
/// requested duration.
///
/// Position is asserted ONLY by kink-aligned quadrature of the reported velocity
/// against the commanded displacement. It is deliberately not sampled at either
/// end of the interval, and not near either end either: both profiles write
/// their final segment as an offset backwards from the commanded displacement,
/// so the end position is returned correctly by construction even by a profile
/// whose traversal is wrong by orders of magnitude.
/// @param duration_ops chained roundings behind the realized duration on the
///        path the profile actually took
/// @param conditioning amplification the solve applies to those roundings; one
///        for a solve whose result is not a canceling difference
template <typename Trajectory>
void check_time_scaling_contract(Trajectory const& profile, double h_signed, scaling_limits const& lim,
                                 double v0, double v1, double T_target, int duration_ops,
                                 double conditioning)
{
    using Scalar = typename Trajectory::scalar_type;
    constexpr double eps = static_cast<double>(std::numeric_limits<Scalar>::epsilon());

    double const T = static_cast<double>(profile.duration());
    CAPTURE(T, T_target, h_signed);
    REQUIRE(std::isfinite(T));
    REQUIRE(T > 0.0);

    // Segment bookkeeping: no phase may be negative, and the reported duration
    // must be what the phases add up to.
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

    // Swept displacement, by kink-aligned quadrature.
    //
    // Two terms. The first is the quadrature's own rounding floor: Simpson is
    // exact on the integrand once the panels lie inside a single segment, so the
    // whole residual there is rounding, at the scale of the absolute area the
    // panels accumulated. The second is the deceleration branch's backward-time
    // recovery: it forms its local time by subtracting the sample time from the
    // duration, which is representable only to one unit in the last place of the
    // duration, and the velocity slews at up to the acceleration limit there, so
    // a sample taken in that branch carries up to a_max * eps * T of velocity
    // error. The quadrature weight exposed to it is the final segment plus one
    // panel of its predecessor, whose right endpoint is the shared boundary.
    auto const q = integrate_velocity(profile);
    double const area_scale = std::max(q.abs_area, std::abs(h_signed));
    double const displacement_tol =
        static_cast<double>(rounding_ops_per_panel * q.panels) * eps * area_scale
        + static_cast<double>(rounding_ops_per_sample) * eps * lim.a_max * T * q.back_time_exposure;
    CAPTURE(q.integral, q.abs_area, q.panels, q.back_time_exposure, displacement_tol);
    REQUIRE(std::abs(q.integral - h_signed) <= displacement_tol);

    // Realized duration, against a budget derived from the solve's own
    // conditioning. Never an exact-equality assertion: neither profile snaps its
    // stored duration to the request, and on the trapezoidal path the cruise
    // velocity is a genuinely ill-conditioned function of the request whenever
    // the shape's residual displacement is a canceling difference.
    double const duration_tol = static_cast<double>(duration_ops) * eps * T_target * conditioning;
    CAPTURE(duration_tol, conditioning);
    REQUIRE(std::abs(T - T_target) <= duration_tol);

    // Boundary velocities. Nothing anchors the velocity the way the final
    // position segment is anchored, so a step taken inside either end is a real
    // check. The step is the profile's own first (respectively last) nonempty
    // segment, and the excursion across it cannot exceed the acceleration limit
    // times its length whatever segment shape it has.
    //
    // Two chains set the margin: the sample itself, at the velocity scale, and
    // the probe time, which for a time measured backwards from the end of the
    // move is representable only to one unit in the last place of the duration
    // and converts into a_max * eps * T of velocity through the slew there.
    double const v_ceiling = std::max({lim.v_max, std::abs(v0), std::abs(v1)});
    double const v_tol =
        static_cast<double>(rounding_ops_per_sample) * eps * (v_ceiling + lim.a_max * T);

    double first_length = 0.0;
    double last_length = 0.0;
    for(auto const& segment : segments)
    {
        double const length = static_cast<double>(segment);
        if(!(length > 0.0))
            continue;
        if(!(first_length > 0.0))
            first_length = length;
        last_length = length;
    }

    if(first_length > 0.0)
    {
        double const v = sample_velocity(profile, first_length);
        CAPTURE(first_length, v);
        REQUIRE(std::abs(v - v0) <= lim.a_max * first_length + v_tol);
    }
    if(last_length > 0.0)
    {
        double const v = sample_velocity(profile, T - last_length);
        CAPTURE(last_length, v);
        REQUIRE(std::abs(v - v1) <= lim.a_max * last_length + v_tol);
    }

    // Kinematic limits and continuity, on a grid laid over the COMMANDED
    // duration. The limit margins are one sample's worth of rounding at each
    // limit's own scale; the continuity bounds are the kinematic bounds on how
    // far each quantity can move in one step of that grid, and the velocity
    // margin carries the same backward-time term as the boundary check.
    //
    // The position bound also carries the conditioning of a phase duration,
    // which the seam between two phases inherits: the cruise velocity is a
    // difference of two nearly equal quantities against the boundary velocities,
    // so its absolute error lives at the velocity scale; dividing it by the
    // acceleration to form a phase duration multiplies that error by the
    // reciprocal of the acceleration, and the position the two phases have to
    // agree on multiplies it by the velocity once more. That is a conditioning
    // limit of the parametrization, not a discontinuity in the profile.
    double const dt = T_target / static_cast<double>(continuity_samples);
    double const seam_conditioning = static_cast<double>(rounding_ops_per_sample) * eps * v_ceiling
                                     * v_ceiling / lim.a_max;
    double const position_bound = v_ceiling * dt + 0.5 * lim.a_max * dt * dt;
    double const velocity_bound = lim.a_max * dt + 0.5 * lim.j_max * dt * dt;
    double const acceleration_bound = lim.j_max * dt;

    double previous_q = 0.0;
    double previous_v = 0.0;
    double previous_a = 0.0;
    for(int i = 0; i <= continuity_samples; ++i)
    {
        double const t = static_cast<double>(i) * dt;
        auto const point = profile.evaluate(static_cast<Scalar>(t));
        double const position = static_cast<double>(point.position(0));
        double const velocity = static_cast<double>(point.velocity(0));
        double const acceleration = static_cast<double>(point.acceleration(0));

        CAPTURE(i, t, position, velocity, acceleration);
        REQUIRE(std::isfinite(position));
        // The velocity margin is the same one the boundary check uses, and for
        // the same reason: a sample taken through the deceleration branch is
        // resolved only to one unit in the last place of the duration, which the
        // slew there turns into a_max * eps * T of velocity.
        REQUIRE(std::abs(velocity) <= v_ceiling + v_tol);
        REQUIRE(std::abs(acceleration)
                <= lim.a_max + static_cast<double>(rounding_ops_per_sample) * eps * lim.a_max);

        if(i > 0)
        {
            double const position_scale =
                std::max({std::abs(position), std::abs(previous_q), std::abs(h_signed)});
            double const position_tol =
                static_cast<double>(rounding_ops_per_sample) * eps * (position_scale + lim.a_max * T * dt)
                + seam_conditioning;
            REQUIRE(std::abs(position - previous_q) <= position_bound + position_tol);
            REQUIRE(std::abs(velocity - previous_v) <= velocity_bound + v_tol);
            if(lim.j_max > 0.0)
            {
                double const acceleration_tol =
                    static_cast<double>(rounding_ops_per_sample) * eps * (lim.a_max + lim.j_max * T);
                REQUIRE(std::abs(acceleration - previous_a) <= acceleration_bound + acceleration_tol);
            }
        }

        previous_q = position;
        previous_v = velocity;
        previous_a = acceleration;
    }
}

/// The three shapes a rescaled trapezoidal profile can take, named by where the
/// cruise velocity lands relative to the two boundary velocities. Read off the
/// profile's own reported peak, folded back into the positive-displacement
/// frame; no internal state is inspected.
enum class cruise_shape
{
    plateau,      ///< cruise at or above both boundary velocities
    ramp_through, ///< cruise strictly between them: both ramps run one way
    valley,       ///< cruise at or below both: the axis dips and climbs back
};

struct shape_census
{
    int plateau{};
    int ramp_through{};
    int valley{};
};

auto classify_cruise_shape(double v_cruise, double v0, double v1) -> cruise_shape
{
    if(v_cruise >= std::max(v0, v1))
        return cruise_shape::plateau;
    if(v_cruise > std::min(v0, v1))
        return cruise_shape::ramp_through;
    return cruise_shape::valley;
}

/// Configurations for the time-scaling families.
///
/// The commanded displacement is drawn with headroom over two separate floors:
/// the trapezoidal constructor's boundary-velocity feasibility adjustment, which
/// would otherwise raise the acceleration above the commanded limit and make the
/// shape residual cancel to the last bit, and the double-S realizability floor.
/// Both floors have their own dedicated coverage elsewhere in this file; here
/// the point is to exercise the retiming, not the constructor's corners.
auto make_scaling_configs(boundary_family family) -> std::array<sweep_config, configs_per_family>
{
    std::mt19937 gen(sweep_seed + 100u + static_cast<std::uint32_t>(family));
    std::uniform_real_distribution<double> q0_dist(-20.0, 20.0);
    std::uniform_real_distribution<double> v_max_dist(1.0, 6.0);
    std::uniform_real_distribution<double> a_max_dist(0.5, 20.0);
    std::uniform_real_distribution<double> j_max_dist(1.0, 80.0);
    std::uniform_real_distribution<double> frac_dist(0.05, 0.9);
    std::uniform_real_distribution<double> reverse_frac_dist(-0.8, -0.05);
    std::uniform_real_distribution<double> decade_dist(-2.0, 2.0);
    std::uniform_real_distribution<double> slack_dist(2.0, 5.0);
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

        double const slack = slack_dist(gen);
        double h = std::pow(10.0, decade_dist(gen));
        h = std::max(h, slack * std::abs(pv0 * pv0 - pv1 * pv1) / (2.0 * a_max));
        h = std::max(h, slack * std::abs(minimum_displacement(pv0, pv1, a_max, j_max)));

        double const direction = flip(gen) ? 1.0 : -1.0;
        double const q0 = q0_dist(gen);
        cfg = {q0, q0 + direction * h, v_max, a_max, j_max, direction * pv0, direction * pv1};
    }
    return configs;
}

/// Request each duration multiple of a trapezoidal profile's own duration and
/// hold every accepted result to the full contract. A rejected request is
/// asserted by its specific enumerator, never merely observed.
template <typename Scalar>
void sweep_trapezoidal_scaling(sweep_config const& cfg, shape_census& census)
{
    for(double multiple : duration_multiples)
    {
        typename trapezoidal_trajectory<Scalar>::config const tcfg{
            .q0 = static_cast<Scalar>(cfg.q0),
            .q1 = static_cast<Scalar>(cfg.q1),
            .v_max = static_cast<Scalar>(cfg.v_max),
            .a_max = static_cast<Scalar>(cfg.a_max),
            .v0 = static_cast<Scalar>(cfg.v0),
            .v1 = static_cast<Scalar>(cfg.v1)};
        auto built = trapezoidal_trajectory<Scalar>::create(tcfg);
        if(!built.has_value())
            continue;

        auto profile = built.value();
        double const T_current = static_cast<double>(profile.duration());
        if(!(T_current > 0.0))
            continue;

        auto const target = static_cast<Scalar>(T_current * multiple);
        CAPTURE(multiple, T_current, static_cast<double>(target));

        auto const rescaled = profile.rescale_to(target);
        if(!rescaled.has_value())
        {
            // The request is longer than the current duration and finite and
            // positive, so the only rejection the contract admits here is the
            // unreachable one.
            REQUIRE(rescaled.error() == trajectory_error::unreachable_duration);
            continue;
        }

        double const h_signed =
            static_cast<double>(tcfg.q1) - static_cast<double>(tcfg.q0);
        double const sigma = (h_signed >= 0.0) ? 1.0 : -1.0;
        double const v_cruise = sigma * static_cast<double>(profile.peak_velocity());
        double const pv0 = sigma * static_cast<double>(tcfg.v0);
        double const pv1 = sigma * static_cast<double>(tcfg.v1);

        switch(classify_cruise_shape(v_cruise, pv0, pv1))
        {
        case cruise_shape::plateau: ++census.plateau; break;
        case cruise_shape::ramp_through: ++census.ramp_through; break;
        case cruise_shape::valley: ++census.valley; break;
        }

        double const conditioning = trapezoidal_solve_conditioning(
            v_cruise, cfg.a_max, std::abs(h_signed), pv0, pv1, static_cast<double>(target));
        CAPTURE(v_cruise, conditioning);

        check_time_scaling_contract(profile, h_signed,
                                    {.v_max = cfg.v_max, .a_max = cfg.a_max, .j_max = 0.0},
                                    static_cast<double>(tcfg.v0), static_cast<double>(tcfg.v1),
                                    static_cast<double>(target), duration_ops_trapezoidal,
                                    conditioning);
    }
}

/// The same for the double-S profile. The rest-to-rest family takes the
/// closed-form path and every other family takes the bracketed solve, so the
/// operation count is selected by the commanded boundary velocities.
template <typename Scalar>
void sweep_double_s_scaling(sweep_config const& cfg)
{
    for(double multiple : duration_multiples)
    {
        typename double_s_trajectory<Scalar>::config const dcfg{
            .q0 = static_cast<Scalar>(cfg.q0),
            .q1 = static_cast<Scalar>(cfg.q1),
            .v_max = static_cast<Scalar>(cfg.v_max),
            .a_max = static_cast<Scalar>(cfg.a_max),
            .j_max = static_cast<Scalar>(cfg.j_max),
            .v0 = static_cast<Scalar>(cfg.v0),
            .v1 = static_cast<Scalar>(cfg.v1)};

        auto created = double_s_trajectory<Scalar>::create(dcfg);
        if(!created.has_value())
            continue;

        auto profile = created.value();
        double const T_current = static_cast<double>(profile.duration());
        if(!(T_current > 0.0))
            continue;

        auto const target = static_cast<Scalar>(T_current * multiple);
        CAPTURE(multiple, T_current, static_cast<double>(target));

        auto const rescaled = profile.rescale_to(target);
        if(!rescaled.has_value())
        {
            REQUIRE(rescaled.error() == trajectory_error::unreachable_duration);
            continue;
        }

        bool const rest_to_rest = (dcfg.v0 == Scalar{0} && dcfg.v1 == Scalar{0});
        int const duration_ops =
            rest_to_rest ? duration_ops_double_s_closed_form : duration_ops_double_s_solved;

        double const h_signed = static_cast<double>(dcfg.q1) - static_cast<double>(dcfg.q0);
        check_time_scaling_contract(
            profile, h_signed, {.v_max = cfg.v_max, .a_max = cfg.a_max, .j_max = cfg.j_max},
            static_cast<double>(dcfg.v0), static_cast<double>(dcfg.v1),
            static_cast<double>(target), duration_ops, 1.0);
    }
}

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

        auto const result = double_s_trajectory<double>::create({.q0 = cfg.q0,
                                                                    .q1 = cfg.q1,
                                                                    .v_max = cfg.v_max,
                                                                    .a_max = cfg.a_max,
                                                                    .j_max = cfg.j_max,
                                                                    .v0 = cfg.v0,
                                                                    .v1 = cfg.v1});
        REQUIRE(!result.has_value());
        REQUIRE(result.error() == trajectory_error::unreachable_boundary_velocity);
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

    auto const rejected = double_s_trajectory<double>::create({.q0 = cfg.q0,
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
    auto const accepted = double_s_trajectory<double>::create({.q0 = widened.q0,
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

// ---------------------------------------------------------------------------
// Time scaling
// ---------------------------------------------------------------------------

TEST_CASE("trapezoidal profiles hold their traversal contract when retimed", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);

    shape_census census{};
    for(auto family : {boundary_family::zero, boundary_family::equal, boundary_family::differing,
                       boundary_family::reversing})
    {
        for(auto const& cfg : make_scaling_configs(family))
        {
            CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.v0, cfg.v1);
            sweep_trapezoidal_scaling<double>(cfg, census);
            sweep_trapezoidal_scaling<float>(cfg, census);
        }
    }

    // All three shapes the retimed profile can take are exercised by the sweep,
    // identified by where the reported cruise velocity landed rather than by
    // reaching into the profile.
    CAPTURE(census.plateau, census.ramp_through, census.valley);
    REQUIRE(census.plateau > 0);
    REQUIRE(census.ramp_through > 0);
    REQUIRE(census.valley > 0);
}

TEST_CASE("double-S profiles hold their traversal contract when retimed", "[trajectory][anchor]")
{
    CAPTURE(sweep_seed);

    for(auto family : {boundary_family::zero, boundary_family::equal, boundary_family::differing,
                       boundary_family::reversing})
    {
        for(auto const& cfg : make_scaling_configs(family))
        {
            CAPTURE(cfg.q0, cfg.q1, cfg.v_max, cfg.a_max, cfg.j_max, cfg.v0, cfg.v1);
            sweep_double_s_scaling<double>(cfg);
            sweep_double_s_scaling<float>(cfg);
        }
    }
}

TEST_CASE("a retimed trapezoidal profile takes each of its three shapes", "[trajectory][anchor]")
{
    // Slowing from 1.5 to 0.5 over ten units at unit acceleration. The duration
    // falls as the cruise velocity rises, so the three shapes occupy three
    // consecutive ranges of the requested duration, and the boundaries between
    // them are the durations realized when the cruise velocity sits exactly on a
    // boundary velocity: about seven at the larger, about nineteen at the
    // smaller.
    sweep_config const cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 3.0, .a_max = 1.0, .j_max = 0.0, .v0 = 1.5, .v1 = 0.5};

    struct expectation
    {
        double target;
        cruise_shape shape;
    };

    for(auto const& expected :
        {expectation{5.0, cruise_shape::plateau}, expectation{10.0, cruise_shape::ramp_through},
         expectation{25.0, cruise_shape::valley}})
    {
        auto built = trapezoidal_trajectory<double>::create({.q0 = cfg.q0,
                                                             .q1 = cfg.q1,
                                                             .v_max = cfg.v_max,
                                                             .a_max = cfg.a_max,
                                                             .v0 = cfg.v0,
                                                             .v1 = cfg.v1});
        REQUIRE(built.has_value());
        auto profile = built.value();
        CAPTURE(expected.target, profile.duration());
        REQUIRE(profile.rescale_to(expected.target).has_value());

        double const v_cruise = profile.peak_velocity();
        CAPTURE(v_cruise);
        REQUIRE(classify_cruise_shape(v_cruise, cfg.v0, cfg.v1) == expected.shape);

        // Each shape also shows itself in the phase structure: a plateau
        // accelerates then decelerates, a valley does the opposite, and a
        // ramp-through runs both ramps the same way.
        auto const segments = profile.phase_durations();
        double const a_start = profile.evaluate(0.5 * segments[0]).acceleration(0);
        double const a_end = profile.evaluate(profile.duration() - 0.5 * segments[2]).acceleration(0);
        CAPTURE(a_start, a_end);
        if(expected.shape == cruise_shape::plateau)
        {
            REQUIRE(a_start > 0.0);
            REQUIRE(a_end < 0.0);
        }
        else if(expected.shape == cruise_shape::valley)
        {
            REQUIRE(a_start < 0.0);
            REQUIRE(a_end > 0.0);
        }
        else
        {
            REQUIRE(a_start < 0.0);
            REQUIRE(a_end < 0.0);
        }

        double const conditioning = trapezoidal_solve_conditioning(
            v_cruise, cfg.a_max, cfg.q1 - cfg.q0, cfg.v0, cfg.v1, expected.target);
        check_time_scaling_contract(profile, cfg.q1 - cfg.q0,
                                    {.v_max = cfg.v_max, .a_max = cfg.a_max, .j_max = 0.0}, cfg.v0,
                                    cfg.v1, expected.target, duration_ops_trapezoidal, conditioning);
    }
}

TEST_CASE("a retimed profile that starts away from its target holds its contract",
          "[trajectory][anchor]")
{
    // The axis begins moving in the direction opposite to its target, so its
    // transformed boundary velocity is negative. The branch algebra covers that
    // sign; this case is what tells a correct evaluator apart from one that is
    // wrong there, rather than leaving the corner invisible.
    sweep_config const cfg{
        .q0 = 0.0, .q1 = 10.0, .v_max = 3.0, .a_max = 1.0, .j_max = 40.0, .v0 = -1.0, .v1 = 0.5};
    REQUIRE(cfg.v0 * (cfg.q1 - cfg.q0) < 0.0);

    shape_census census{};
    sweep_trapezoidal_scaling<double>(cfg, census);
    sweep_trapezoidal_scaling<float>(cfg, census);
    CAPTURE(census.plateau, census.ramp_through, census.valley);
    REQUIRE(census.plateau + census.ramp_through + census.valley > 0);

    // The valley shape needs a cruise velocity at or below both boundary
    // velocities, and no positive cruise velocity can sit below a negative one,
    // so this family reaches only the other two shapes.
    REQUIRE(census.valley == 0);

    sweep_double_s_scaling<double>(cfg);
    sweep_double_s_scaling<float>(cfg);
}

TEST_CASE("retiming rejects a duration shorter than the profile already takes", "[trajectory][anchor]")
{
    auto trapezoidal_built = trapezoidal_trajectory<double>::create(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .v0 = 0.5, .v1 = 0.25});
    REQUIRE(trapezoidal_built.has_value());
    auto trapezoidal = trapezoidal_built.value();
    auto const trapezoidal_before = trapezoidal.duration();
    auto const trapezoidal_result = trapezoidal.rescale_to(0.5 * trapezoidal_before);
    REQUIRE(!trapezoidal_result.has_value());
    REQUIRE(trapezoidal_result.error() == trajectory_error::duration_shorter_than_current);
    REQUIRE(trapezoidal.duration() == trapezoidal_before);

    auto created = double_s_trajectory<double>::create(
        {.q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0, .v0 = 0.5, .v1 = 0.25});
    REQUIRE(created.has_value());
    auto double_s = created.value();
    auto const double_s_before = double_s.duration();
    auto const double_s_result = double_s.rescale_to(0.5 * double_s_before);
    REQUIRE(!double_s_result.has_value());
    REQUIRE(double_s_result.error() == trajectory_error::duration_shorter_than_current);
    REQUIRE(double_s.duration() == double_s_before);
}

TEST_CASE("trapezoidal retiming rejects a duration past its own reachable maximum",
          "[trajectory][anchor]")
{
    // Both boundary velocities are positive and the commanded displacement sits
    // below the distance the two ramps sweep on their own, so the shape that
    // reaches the vanishing-cruise limit is the valley and its supremum is
    // finite: the cruise duration runs out before the duration can diverge.
    sweep_config const cfg{
        .q0 = 0.0, .q1 = 1.0, .v_max = 3.0, .a_max = 1.0, .j_max = 0.0, .v0 = 2.0, .v1 = 2.0};

    double const h = cfg.q1 - cfg.q0;
    double const v_sum_sq = cfg.v0 * cfg.v0 + cfg.v1 * cfg.v1;
    double const residual = h - v_sum_sq / (2.0 * cfg.a_max);
    REQUIRE(residual < 0.0);

    // Vanishing-cruise duration: the cruise velocity reaches
    // sqrt((v0^2 + v1^2)/2 - a h) and both ramps run from their boundary
    // velocity down to it and back up.
    double const v_min = std::sqrt(-cfg.a_max * residual);
    double const T_sup = (cfg.v0 + cfg.v1 - 2.0 * v_min) / cfg.a_max;
    CAPTURE(v_min, T_sup);

    auto make = [&] {
        auto built = trapezoidal_trajectory<double>::create({.q0 = cfg.q0,
                                                            .q1 = cfg.q1,
                                                            .v_max = cfg.v_max,
                                                            .a_max = cfg.a_max,
                                                            .v0 = cfg.v0,
                                                            .v1 = cfg.v1});
        REQUIRE(built.has_value());
        return built.value();
    };

    auto rejected = make();
    REQUIRE(T_sup > rejected.duration());
    auto const beyond = rejected.rescale_to(T_sup * 1.1);
    REQUIRE(!beyond.has_value());
    REQUIRE(beyond.error() == trajectory_error::unreachable_duration);

    // Just inside the supremum the request is served, which is what makes the
    // rejection above a boundary rather than a blanket refusal.
    auto accepted = make();
    double const inside = 0.5 * (accepted.duration() + T_sup);
    REQUIRE(accepted.rescale_to(inside).has_value());
    double const conditioning = trapezoidal_solve_conditioning(accepted.peak_velocity(), cfg.a_max, h,
                                                               cfg.v0, cfg.v1, inside);
    check_time_scaling_contract(accepted, h, {.v_max = cfg.v_max, .a_max = cfg.a_max, .j_max = 0.0},
                                cfg.v0, cfg.v1, inside, duration_ops_trapezoidal, conditioning);
}

TEST_CASE("double-S retiming rejects a duration past its velocity-scale lower bound",
          "[trajectory][anchor]")
{
    // The retiming slows the profile by dividing its velocity limit by the scale
    // it solves for, while the commanded boundary velocities stay fixed. The
    // scale therefore cannot fall below the point where those boundary
    // velocities themselves reach the scaled limit, and here that point is close
    // to one, so only a narrow band of longer durations is reachable.
    sweep_config const cfg{
        .q0 = 0.0, .q1 = 5.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0, .v0 = 4.5, .v1 = 0.0};
    double const scale_floor = std::max(std::abs(cfg.v0), std::abs(cfg.v1)) / cfg.v_max;
    CAPTURE(scale_floor);
    REQUIRE(scale_floor < 1.0);

    auto created = double_s_trajectory<double>::create({.q0 = cfg.q0,
                                                            .q1 = cfg.q1,
                                                            .v_max = cfg.v_max,
                                                            .a_max = cfg.a_max,
                                                            .j_max = cfg.j_max,
                                                            .v0 = cfg.v0,
                                                            .v1 = cfg.v1});
    REQUIRE(created.has_value());

    auto profile = created.value();
    double const T_current = profile.duration();
    CAPTURE(T_current);

    auto const beyond = profile.rescale_to(T_current / (scale_floor * scale_floor));
    REQUIRE(!beyond.has_value());
    REQUIRE(beyond.error() == trajectory_error::unreachable_duration);
    REQUIRE(profile.duration() == T_current);
}

TEST_CASE("retiming rejects a profile with nothing to traverse", "[trajectory][anchor]")
{
    auto trapezoidal_built =
        trapezoidal_trajectory<double>::create({.q0 = 2.0, .q1 = 2.0, .v_max = 5.0, .a_max = 10.0});
    REQUIRE(trapezoidal_built.has_value());
    auto trapezoidal = trapezoidal_built.value();
    REQUIRE(trapezoidal.duration() == 0.0);
    auto const trapezoidal_result = trapezoidal.rescale_to(1.0);
    REQUIRE(!trapezoidal_result.has_value());
    REQUIRE(trapezoidal_result.error() == trajectory_error::unreachable_duration);
    REQUIRE(trapezoidal.duration() == 0.0);

    auto double_s_built = double_s_trajectory<double>::create(
        {.q0 = 2.0, .q1 = 2.0, .v_max = 5.0, .a_max = 10.0, .j_max = 100.0});
    REQUIRE(double_s_built.has_value());
    auto double_s = double_s_built.value();
    REQUIRE(double_s.duration() == 0.0);
    auto const double_s_result = double_s.rescale_to(1.0);
    REQUIRE(!double_s_result.has_value());
    REQUIRE(double_s_result.error() == trajectory_error::unreachable_duration);
    REQUIRE(double_s.duration() == 0.0);
}

TEST_CASE("trapezoidal retiming rejects a valley root that lands on its own shape boundary",
          "[trajectory][anchor]")
{
    // Both boundary velocities sit within a few units in the last place of the
    // velocity limit while the acceleration is ten orders of magnitude smaller,
    // so the valley quadratic's discriminant b^2 - 4c is a difference of two
    // quantities that agree to the full width of the significand. It cancels to
    // exactly zero and the closed form returns the smaller boundary velocity
    // itself, a root with no significant digits in it.
    //
    // That root is not a solution. The valley shape is selected only when the
    // request exceeds the duration the smaller boundary velocity already
    // realizes, and the total duration is strictly decreasing in the cruise
    // velocity, so the root answering such a request lies STRICTLY below that
    // boundary. Accepting the boundary itself would hand back a profile that
    // realizes the same duration for every request past it while reporting
    // success, which is a silently wrong retiming. It is a typed rejection, and
    // the profile is left bitwise untouched.
    auto built = trapezoidal_trajectory<double>::create({.q0 = 0.0,
                                                         .q1 = 0.0004425048828125,
                                                         .v_max = 0.9999999999999996,
                                                         .a_max = 1e-6,
                                                         .v0 = 0.9999999999999994,
                                                         .v1 = 0.9999999050046706});
    REQUIRE(built.has_value());
    auto profile = built.value();

    auto const T_current = profile.duration();
    auto const phases_before = profile.phase_durations();
    CAPTURE(T_current);

    auto const rejected = profile.rescale_to(T_current * 1.0000000002328306);
    REQUIRE(!rejected.has_value());
    REQUIRE(rejected.error() == trajectory_error::unreachable_duration);
    REQUIRE(profile.duration() == T_current);
    REQUIRE(profile.phase_durations() == phases_before);
}
