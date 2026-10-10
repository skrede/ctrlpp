#ifndef HPP_GUARD_CTRLPP_TESTS_FUZZ_ONLINE_PLANNER_3RD_CASE_H
#define HPP_GUARD_CTRLPP_TESTS_FUZZ_ONLINE_PLANNER_3RD_CASE_H

#include "ctrlpp/trajectory/online_planner_3rd.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace ctrlpp::fuzz
{

// Dense scan: many fixed-size steps exercise every jerk/accel/cruise/decel
// segment (and the settled tail) regardless of the fuzzed target distance.
inline constexpr std::int32_t online_planner_3rd_samples = 2000;

// Rounding-op margin: sample() integrates through up to 11 constant-jerk
// phases, each a chained multiply-add (cubic position, quadratic
// velocity, linear acceleration) plus phase-boundary subtractions, so a
// single sample can accumulate several ULPs of rounding at its own
// operand scale, not just the single operation a bare epsilon assumes;
// 32 is a generous round count of those chained operations.
inline constexpr double online_planner_3rd_rounding_op_margin = 32.0;

inline constexpr double online_planner_3rd_largest_limit = 1e6;
inline constexpr double online_planner_3rd_largest_dt = 1.0;

template <typename Scalar>
struct online_planner_3rd_case
{
    Scalar pos;
    Scalar target;
    Scalar v_max;
    Scalar a_max;
    Scalar j_max;
    Scalar dt;
    Scalar t_start;
    Scalar retarget_target;
    std::int32_t retarget_index;
};

// The largest admitted a_max over the spacing of representable times at the largest sample time the scan forms,
// 2 * samples * dt_max, so that a_max / j_max reaches below the sample-time resolution.
template <typename Scalar>
auto online_planner_3rd_largest_jerk() -> double
{
    auto const t_max = static_cast<Scalar>(2.0 * online_planner_3rd_samples * online_planner_3rd_largest_dt);
    auto const spacing = std::nextafter(t_max, std::numeric_limits<Scalar>::infinity()) - t_max;
    return static_cast<double>(static_cast<Scalar>(online_planner_3rd_largest_limit) / spacing);
}

template <typename Scalar>
auto narrow_field(double value, Scalar& out) -> bool
{
    out = static_cast<Scalar>(value);
    return std::isfinite(out);
}

// 9 doubles: pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target = 72 bytes.
// The clamps run in double and each field is then narrowed to the scalar under test.
template <typename Scalar>
auto decode_online_planner_3rd_case(const std::uint8_t* data, std::size_t size, online_planner_3rd_case<Scalar>& out) -> bool
{
    if (size < 72)
        return false;
    double buf[9];
    std::memcpy(buf, data, 72);
    for (int i = 0; i < 9; ++i) {
        if (!std::isfinite(buf[i]))
            return false;
    }
    constexpr double largest = online_planner_3rd_largest_limit;
    double const dt = std::clamp(std::abs(buf[5]), 1e-6, online_planner_3rd_largest_dt);
    // Clamp the scan's absolute time origin to within one scan window's worth
    // of the reference time (t_ref_ = 0): sample() computes elapsed time as
    // t - t_ref_, and an absolute t far larger than the dt-sized increments
    // being added to it would carry that increment below t's own floating-
    // point resolution, corrupting the finite differences below with pure
    // representation noise rather than trajectory behavior.
    double const t_start_bound = static_cast<double>(online_planner_3rd_samples) * dt;
    out.retarget_index = static_cast<std::int32_t>(std::clamp(std::abs(buf[7]), 0.0, 1.0) * (online_planner_3rd_samples - 1));
    return narrow_field(buf[0], out.pos) && narrow_field(buf[1], out.target)
           && narrow_field(std::clamp(std::abs(buf[2]), 1e-3, largest), out.v_max)
           && narrow_field(std::clamp(std::abs(buf[3]), 1e-3, largest), out.a_max)
           && narrow_field(std::clamp(std::abs(buf[4]), 1e-3, online_planner_3rd_largest_jerk<Scalar>()), out.j_max)
           && narrow_field(dt, out.dt) && narrow_field(std::clamp(buf[6], -t_start_bound, t_start_bound), out.t_start)
           && narrow_field(buf[8], out.retarget_target);
}

template <typename Scalar>
struct online_planner_3rd_state
{
    Scalar q;
    Scalar v;
    Scalar a;
};

// Envelope: the analytically-returned velocity/acceleration must not
// exceed the configured limits beyond a rounding-scaled tolerance.
// The velocity check adds an a_max*dt physical term and the
// acceleration check adds a j_max*dt physical term: when a phase
// transition falls within a single dt-sized sampling step (coarse
// sampling relative to that phase's own, possibly tiny, duration),
// the returned value can differ from the clean limit/0 value by up
// to one step's worth of the next-higher derivative's limit.
// After the retarget, |v| and |a| are held to the limits with no tolerance,
// because the planner evaluates each phase-group end from its analytic state
// instead of accumulating rounding.
template <typename Scalar>
void check_online_planner_3rd_envelope(online_planner_3rd_case<Scalar> const& c, online_planner_3rd_state<Scalar> const& s,
                                       Scalar time_resolution_floor, bool after_retarget)
{
    constexpr Scalar eps = std::numeric_limits<Scalar>::epsilon();
    constexpr auto margin = static_cast<Scalar>(online_planner_3rd_rounding_op_margin);
    constexpr auto samples = static_cast<Scalar>(online_planner_3rd_samples);
    const Scalar v_env_tol = c.a_max * c.dt + eps * margin * c.v_max * samples + time_resolution_floor;
    const Scalar a_env_tol = c.j_max * c.dt + eps * margin * c.a_max * samples + time_resolution_floor;
    if (std::abs(s.v) > c.v_max + v_env_tol)
        std::abort();
    if (std::abs(s.a) > c.a_max + a_env_tol)
        std::abort();
    if (after_retarget && (std::abs(s.v) > c.v_max || std::abs(s.a) > c.a_max))
        std::abort();
}

template <typename Scalar>
void check_online_planner_3rd_differences(online_planner_3rd_case<Scalar> const& c, online_planner_3rd_state<Scalar> const& prev,
                                          online_planner_3rd_state<Scalar> const& s, Scalar time_resolution_floor)
{
    constexpr Scalar eps = std::numeric_limits<Scalar>::epsilon();
    constexpr auto margin = static_cast<Scalar>(online_planner_3rd_rounding_op_margin);
    constexpr auto samples = static_cast<Scalar>(online_planner_3rd_samples);
    // Continuity: position cannot move faster than the kinematic bound
    // set by bounded velocity and acceleration over one sample step
    // (this bound only uses the v/a limits, so it holds regardless of
    // how the acceleration itself varies within the step).
    const Scalar continuity_bound = c.v_max * c.dt + Scalar{0.5} * c.a_max * c.dt * c.dt;
    const Scalar continuity_tol = eps * margin * (std::abs(s.q) + std::abs(prev.q) + continuity_bound);
    if (std::abs(s.q - prev.q) > continuity_bound + continuity_tol)
        std::abort();

    // Finite-difference velocity from the position trace: the average
    // velocity over one step differs from an instantaneous sample by
    // at most a_max*dt, since velocity's own rate of change (the
    // acceleration) never exceeds a_max. A sample pair straddling a
    // phase boundary is evaluated by two distinct closed-form
    // branches that meet exactly only in real arithmetic, so an
    // extra term scaled by the position magnitude in play (the same
    // rounding-op margin, divided back out by dt to land in velocity
    // units) absorbs that seam's rounding.
    const Scalar v_fd = (s.q - prev.q) / c.dt;
    const Scalar seam_v_tol = eps * margin * (std::abs(s.q) + std::abs(prev.q)) / c.dt;
    const Scalar v_fd_tol = c.a_max * c.dt + eps * margin * c.v_max * samples + time_resolution_floor + seam_v_tol;
    if (std::abs(v_fd) > c.v_max + v_fd_tol)
        std::abort();

    // Finite-difference acceleration from the velocity trace: the
    // average acceleration over one step differs from an
    // instantaneous sample by at most j_max*dt, since acceleration's
    // own rate of change (the jerk) never exceeds j_max, plus the
    // same phase-boundary seam term one derivative order up (scaled
    // by the velocity magnitude in play).
    const Scalar a_fd = (s.v - prev.v) / c.dt;
    const Scalar seam_a_tol = eps * margin * (std::abs(s.v) + std::abs(prev.v)) / c.dt;
    const Scalar a_fd_tol = c.j_max * c.dt + eps * margin * c.a_max * samples + time_resolution_floor + seam_a_tol;
    if (std::abs(a_fd) > c.a_max + a_fd_tol)
        std::abort();

    // Finite-difference jerk from the acceleration trace: the
    // profile's jerk is piecewise constant, so the average jerk over
    // one sample step cannot exceed j_max beyond rounding, plus the
    // same phase-boundary seam term one more derivative order up
    // (scaled by the acceleration magnitude in play).
    const Scalar j_fd = (s.a - prev.a) / c.dt;
    const Scalar seam_j_tol = eps * margin * (std::abs(s.a) + std::abs(prev.a)) / c.dt;
    const Scalar j_fd_tol = eps * margin * c.j_max * samples + time_resolution_floor + seam_j_tol;
    if (std::abs(j_fd) > c.j_max + j_fd_tol)
        std::abort();
}

// Reject configurations whose position scale swamps the per-step motion:
// if a handful of ULPs of pos/target's own magnitude (the same
// rounding-op margin used below) already exceeds the kinematic bound on
// how far the position can move in a single sample step, adding that
// step's motion to a position of this magnitude rounds away completely,
// and no continuity signal survives to check -- this is a floating-point
// representation limit, not a planner defect.
template <typename Scalar>
auto online_planner_3rd_position_resolvable(online_planner_3rd_case<Scalar> const& c) -> bool
{
    constexpr Scalar eps = std::numeric_limits<Scalar>::epsilon();
    const Scalar position_scale = std::max({std::abs(c.pos), std::abs(c.target), std::abs(c.retarget_target), Scalar{1}});
    const Scalar continuity_bound_reject = c.v_max * c.dt + Scalar{0.5} * c.a_max * c.dt * c.dt;
    return eps * static_cast<Scalar>(online_planner_3rd_rounding_op_margin) * position_scale < continuity_bound_reject;
}

template <typename Scalar>
auto run_online_planner_3rd_case(const std::uint8_t* data, std::size_t size) -> int
{
    online_planner_3rd_case<Scalar> c{};
    if (!decode_online_planner_3rd_case(data, size, c) || !online_planner_3rd_position_resolvable(c))
        return 0;

    // Sample-time floor: advancing through num_samples discrete sample times
    // carries its own fixed rounding floor of about one part in num_samples of
    // machine epsilon once divided back out by dt in a finite difference,
    // independent of which physical limit (v_max/a_max/j_max) is being
    // checked; it only matters once that limit is itself very small.
    const Scalar time_resolution_floor = std::numeric_limits<Scalar>::epsilon()
                                         * static_cast<Scalar>(online_planner_3rd_rounding_op_margin)
                                         * static_cast<Scalar>(online_planner_3rd_samples);

    // Construction goes through create: a limit outside the finite,
    // strictly positive domain is rejected by design (each limit divides in
    // the planner math), and that rejection is the correct behavior for such
    // an input, not a crash. The clamps above keep the fuzzed limits inside
    // the domain, so the scan below runs for every input that reaches here.
    auto planner_result = ctrlpp::online_planner_3rd<Scalar>::create({.v_max = c.v_max, .a_max = c.a_max, .j_max = c.j_max});
    if (!planner_result.has_value())
        return 0;
    auto& planner = *planner_result;

    // The retarget replans mid-motion from the state the planner last reported, as a servo loop does.
    planner.reset(c.pos);
    if (!planner.update(c.target).has_value())
        return 0;

    online_planner_3rd_state<Scalar> prev{};
    Scalar t = c.t_start;
    for (std::int32_t i = 0; i < online_planner_3rd_samples; ++i) {
        t += c.dt;
        auto const pt = planner.sample(t);
        online_planner_3rd_state<Scalar> const s{pt.position(0), pt.velocity(0), pt.acceleration(0)};
        if (!std::isfinite(s.q) || !std::isfinite(s.v) || !std::isfinite(s.a))
            return 0;
        check_online_planner_3rd_envelope(c, s, time_resolution_floor, i > c.retarget_index);
        if (i > 0)
            check_online_planner_3rd_differences(c, prev, s, time_resolution_floor);
        prev = s;
        if (i == c.retarget_index && !planner.update(c.retarget_target).has_value())
            return 0;
    }
    return 0;
}

}

#endif
