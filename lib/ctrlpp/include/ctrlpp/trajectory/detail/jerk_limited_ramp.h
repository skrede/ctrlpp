#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_JERK_LIMITED_RAMP_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_JERK_LIMITED_RAMP_H

#include <cmath>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

/// @brief Shape of a jerk-limited ramp between two velocities whose endpoint
/// accelerations are both zero.
template <typename Scalar>
struct velocity_ramp
{
    Scalar T;     ///< total ramp duration
    Scalar T_j;   ///< duration of each constant-jerk sub-segment
    Scalar a_lim; ///< acceleration magnitude actually reached
};

// The longest jerk time whose product with j_max does not round above a_max. One step suffices for a finite
// quotient, because round-to-nearest leaves the quotient's predecessor below a_max / j_max.
template <typename Scalar>
auto jerk_time_within(Scalar a_max, Scalar j_max) -> Scalar
{
    auto const T_j = a_max / j_max;
    if (j_max * T_j > a_max) {
        return std::nextafter(T_j, Scalar{0});
    }
    return T_j;
}

// The longest jerk time from sqrt(v_max) / sqrt(j_max) down whose peak j_max T_j^2 does not round above v_max. With normal operands the
// rounded peak is at most v_max (1+u)^6 / (1-u)^2 at the start and each step scales it by at most (1-u)^2, so four steps suffice.
template <typename Scalar>
auto jerk_time_within_velocity(Scalar v_max, Scalar j_max) -> Scalar
{
    auto T_j = std::sqrt(v_max) / std::sqrt(j_max);
    for (std::int32_t step = 0; step < 4 && j_max * T_j * T_j > v_max; ++step) {
        T_j = std::nextafter(T_j, Scalar{0});
    }
    return T_j;
}

/// @brief Jerk time of the rest-to-rest profile that reaches neither limit, h = 2 j_max T_j^3.
/// @cite biagiotti2009 -- Sec. 3.4.3, p.90-91
template <typename Scalar>
auto doubly_degenerate_jerk_time(Scalar h, Scalar v_max, Scalar a_max, Scalar j_max) -> Scalar
{
    auto const T_j = std::min(std::cbrt(h / (Scalar{2} * j_max)), jerk_time_within(a_max, j_max));
    return std::min(T_j, jerk_time_within_velocity(v_max, j_max));
}

// Rest-to-rest with a_max reached and no cruise: c = v_max T_a - h >= 0, and the peak's deficit d below v_max is the small root of
// d^2 / a_max - d (2 v_max / a_max + a_max / j_max) + c = 0. The peak is formed from the nearer of v_max and zero.
template <typename Scalar>
auto no_cruise_peak(Scalar c, Scalar h, Scalar v_max, Scalar a_max, Scalar j_max) -> Scalar
{
    auto const a = Scalar{1} / a_max;
    auto const b = a_max / j_max;
    auto const root = std::sqrt(b * b + Scalar{4} * a * h);
    auto const deficit = Scalar{2} * c / (Scalar{2} * v_max * a + b + root);
    if (deficit + deficit <= v_max) {
        return v_max - deficit;
    }
    return (-b + root) / (Scalar{2} * a);
}

/// @brief Fastest admissible ramp spanning a velocity change.
///
/// The ramp is jerk, constant acceleration, jerk. When the velocity change is
/// too small to build the acceleration up to a_max the constant-acceleration
/// sub-segment vanishes, the ramp is triangular in acceleration, and it peaks
/// at j_max * T_j instead. Both branches are exact closed forms; there is no
/// sub-case left over for which the acceleration bound has to be lowered and
/// the solve retried.
///
/// @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.20), p.80
template <typename Scalar>
auto make_velocity_ramp(Scalar dv, Scalar a_max, Scalar j_max) -> velocity_ramp<Scalar>
{
    auto const d = std::abs(dv);
    if (d * j_max < a_max * a_max) {
        // Analytically the triangular jerk time is below a_max / j_max, so the bound trims only the rounding that
        // carried it across the branch test.
        auto const T_j = std::min(std::sqrt(d / j_max), jerk_time_within(a_max, j_max));
        return {Scalar{2} * T_j, T_j, j_max * T_j};
    }
    auto const T_j = a_max / j_max;
    return {d / a_max + T_j, T_j, a_max};
}

/// @brief Distance swept by the two ramps of a cruise-free profile whose peak
/// velocity rises by `rise` above the larger boundary velocity.
///
/// Each ramp's acceleration is point-symmetric about the ramp's own midpoint,
/// so the velocity it sweeps has mean exactly the average of its two
/// endpoints and the swept distance is that mean times the ramp duration.
/// Derived in place from that symmetry, which is the same argument the double-S
/// profile's compute_q_at_Ta() rests on. Both mean velocities are formed from
/// the rise rather than from the peak so that a rise far below the boundary
/// velocities keeps every digit it has.
template <typename Scalar>
auto no_cruise_displacement(Scalar rise, Scalar pv_hi, Scalar pv_lo, Scalar a_max, Scalar j_max) -> Scalar
{
    auto const far_rise = rise + (pv_hi - pv_lo);
    auto const near_ramp = make_velocity_ramp(rise, a_max, j_max);
    auto const far_ramp = make_velocity_ramp(far_rise, a_max, j_max);
    return (Scalar{2} * pv_hi + rise) * near_ramp.T / Scalar{2}
         + (Scalar{2} * pv_lo + far_rise) * far_ramp.T / Scalar{2};
}

}

#endif
