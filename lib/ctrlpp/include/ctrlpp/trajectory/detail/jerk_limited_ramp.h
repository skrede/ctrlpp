#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_JERK_LIMITED_RAMP_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_JERK_LIMITED_RAMP_H

#include <cmath>

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
        auto const T_j = std::sqrt(d / j_max);
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
