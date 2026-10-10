#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_RAMP_ANCHOR_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_RAMP_ANCHOR_H

namespace ctrlpp::detail
{

// From the nearer end the sample steps back toward the other end, and monotone rounding cannot carry it past the
// end it steps from; the far end is out of reach because the exact value lies short of the ramp's midpoint.
template <typename Scalar>
auto nearer_anchor_velocity(Scalar v_start, Scalar v_end, Scalar accel, Scalar t_from_start, Scalar t_to_end) -> Scalar
{
    return (t_from_start < t_to_end) ? v_start + accel * t_from_start : v_end - accel * t_to_end;
}

template <typename Scalar>
auto nearer_anchor_position(Scalar q_start, Scalar v_start, Scalar q_end, Scalar v_end, Scalar accel,
                            Scalar t_from_start, Scalar t_to_end) -> Scalar
{
    if (t_from_start < t_to_end) {
        return q_start + v_start * t_from_start + Scalar{0.5} * accel * t_from_start * t_from_start;
    }
    return q_end - v_end * t_to_end + Scalar{0.5} * accel * t_to_end * t_to_end;
}

}

#endif
