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

// The planner stores at most eleven phases, so a sample's phase start is a sum of at most ten durations.
inline constexpr std::int32_t online_planner_3rd_phase_starts = 10;

template <typename Scalar>
inline constexpr Scalar online_planner_3rd_unit_roundoff = std::numeric_limits<Scalar>::epsilon() / Scalar{2};

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
    Scalar t;
    Scalar q;
    Scalar v;
    Scalar a;
};

template <typename Scalar>
struct online_planner_3rd_plan
{
    Scalar t_ref;
    Scalar position_scale;
};

template <typename Scalar>
auto online_planner_3rd_within_limits(online_planner_3rd_case<Scalar> const& c, online_planner_3rd_state<Scalar> const& s) -> bool
{
    return std::abs(s.v) <= c.v_max && std::abs(s.a) <= c.a_max;
}

// Two samples in one phase share its rounded start time, so a pair adds one rounding of the elapsed time per phase
// start between them, at most ten, and three for the comparison itself.
template <typename Scalar>
auto online_planner_3rd_difference_bound(Scalar rate, Scalar dt, Scalar elapsed, Scalar previous_rounding,
                                         Scalar current_rounding) -> Scalar
{
    constexpr auto phase_starts = static_cast<Scalar>(online_planner_3rd_phase_starts);
    return rate * dt + previous_rounding + current_rounding
           + online_planner_3rd_unit_roundoff<Scalar> * rate * (phase_starts * elapsed + Scalar{3} * dt);
}

// An acceleration carries one rounding of the elapsed time and at most 9 of a_max, since a jerk phase lasts at most
// a_max / j_max: the in-phase and time-to-end subtractions, its product and sum, a carried phase's product and sum, and
// a stored duration's quotient and its one-ulp trim, worth two (Higham 2002, Lemma 3.1).
template <typename Scalar>
auto online_planner_3rd_acceleration_rounding(Scalar a_max, Scalar j_max, Scalar elapsed) -> Scalar
{
    return online_planner_3rd_unit_roundoff<Scalar> * (j_max * elapsed + Scalar{9} * a_max);
}

// A velocity carries one rounding of the elapsed time and at most 42 of v_max, since a phase changes the velocity by at
// most 2 v_max and so lasts at most 4 v_max / |a|: 5 forming it, 7 for a carried phase, 8 for the in-phase and
// time-to-end subtractions, and 22 for a ramp's stored durations reproducing its velocity change (Higham 2002, Lemma 3.1).
template <typename Scalar>
auto online_planner_3rd_velocity_rounding(Scalar v_max, Scalar a_max, Scalar elapsed) -> Scalar
{
    return online_planner_3rd_unit_roundoff<Scalar> * (a_max * elapsed + Scalar{42} * v_max);
}

// A position carries one rounding of the elapsed time d, 11 sums of positions within P0 + v_max tau, 51 roundings of
// v_max tau for travel, carried velocity and acceleration, and subtractions, and 20 of |h| + v_max tau for the stored
// durations sweeping |h| <= 2 P0 (Higham 2002, Lemma 3.1). tau = 2 d + v_max / a_max + 2 a_max / j_max bounds the time
// traversed, since a target-anchored position lies past the cruise midpoint and a deceleration lasts at most
// v_max / a_max + a_max / j_max, and the duration of any brake and nulling overshoot.
template <typename Scalar>
auto online_planner_3rd_position_rounding(Scalar v_max, Scalar a_max, Scalar j_max, Scalar elapsed, Scalar position_scale)
    -> Scalar
{
    auto const traversed = Scalar{2} * elapsed + v_max / a_max + Scalar{2} * a_max / j_max;
    return online_planner_3rd_unit_roundoff<Scalar>
           * (v_max * elapsed + Scalar{51} * position_scale + Scalar{82} * v_max * traversed);
}

// The mean of a derivative over a step stays within the derivative's limit, so each difference is held to the limit
// times the step and the counted rounding alone.
template <typename Scalar>
auto online_planner_3rd_differences_hold(online_planner_3rd_case<Scalar> const& c, online_planner_3rd_plan<Scalar> const& plan,
                                         online_planner_3rd_state<Scalar> const& prev, online_planner_3rd_state<Scalar> const& s) -> bool
{
    auto const dt = s.t - prev.t;
    auto const elapsed = std::max(s.t - plan.t_ref, Scalar{0});
    auto const elapsed_prev = std::max(prev.t - plan.t_ref, Scalar{0});
    auto const position = online_planner_3rd_difference_bound(
        c.v_max, dt, elapsed, online_planner_3rd_position_rounding(c.v_max, c.a_max, c.j_max, elapsed_prev, plan.position_scale),
        online_planner_3rd_position_rounding(c.v_max, c.a_max, c.j_max, elapsed, plan.position_scale));
    auto const velocity = online_planner_3rd_difference_bound(c.a_max, dt, elapsed,
                                                              online_planner_3rd_velocity_rounding(c.v_max, c.a_max, elapsed_prev),
                                                              online_planner_3rd_velocity_rounding(c.v_max, c.a_max, elapsed));
    auto const acceleration = online_planner_3rd_difference_bound(
        c.j_max, dt, elapsed, online_planner_3rd_acceleration_rounding(c.a_max, c.j_max, elapsed_prev),
        online_planner_3rd_acceleration_rounding(c.a_max, c.j_max, elapsed));
    return std::abs(s.q - prev.q) <= position && std::abs(s.v - prev.v) <= velocity && std::abs(s.a - prev.a) <= acceleration;
}

// Reject configurations whose position scale swamps the per-step motion:
// if a handful of ULPs of pos/target's own magnitude (the
// rounding-op margin) already exceeds the kinematic bound on
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

    online_planner_3rd_plan<Scalar> plan{Scalar{0}, std::max(std::abs(c.pos), std::abs(c.target))};
    online_planner_3rd_state<Scalar> prev{};
    Scalar t = c.t_start;
    for (std::int32_t i = 0; i < online_planner_3rd_samples; ++i) {
        t += c.dt;
        auto const pt = planner.sample(t);
        online_planner_3rd_state<Scalar> const s{t, pt.position(0), pt.velocity(0), pt.acceleration(0)};
        if (!std::isfinite(s.q) || !std::isfinite(s.v) || !std::isfinite(s.a))
            return 0;
        if (!online_planner_3rd_within_limits(c, s) || (i > 0 && !online_planner_3rd_differences_hold(c, plan, prev, s)))
            std::abort();
        prev = s;
        if (i == c.retarget_index) {
            if (!planner.update(c.retarget_target).has_value())
                return 0;
            plan = {t, std::max(std::abs(s.q), std::abs(c.retarget_target))};
        }
    }
    return 0;
}

}

#endif
