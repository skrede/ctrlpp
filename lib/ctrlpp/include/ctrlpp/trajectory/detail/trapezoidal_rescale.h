#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRAPEZOIDAL_RESCALE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRAPEZOIDAL_RESCALE_H

#include "ctrlpp/trajectory/detail/trapezoidal_state.h"
#include "ctrlpp/trajectory/detail/trapezoidal_rescale_shapes.h"

#include "ctrlpp/trajectory/trajectory_error.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
class trapezoidal_rescale : protected trapezoidal_rescale_shapes<Scalar>
{
protected:
    using trapezoidal_rescale_shapes<Scalar>::trapezoidal_rescale_shapes;

    // Works in the positive-displacement frame, with a ramp's swept distance taken as the exact trapezoid area,
    // since |v^2 - v0^2| / (2a) equals it only when the two velocities sum to a positive value. The total
    // duration is strictly decreasing in the cruise velocity, so the shape is identified by evaluating it at
    // the two shape boundaries rather than by trying each shape in turn.
    // @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.13a)-(3.13c), p.71 -- the
    //   three-phase parametrization these expressions solve for the cruise velocity
    ctrlpp::expected<trapezoidal_rescaled_state<Scalar>, trajectory_error> solve_rescale(Scalar T_new) const
    {
        // Exact equality first: a synchronized set passes the slowest axis a bit-exact copy of its own
        // duration, so that axis lands here rather than in the shortening rejection below.
        if(T_new == this->m_duration)
        {
            return this->current_state();
        }
        if(!std::isfinite(T_new) || T_new <= Scalar{0})
        {
            return ctrlpp::unexpected(trajectory_error::non_positive_duration);
        }
        if(T_new < this->m_duration)
        {
            return ctrlpp::unexpected(trajectory_error::duration_shorter_than_current);
        }

        auto const h = std::abs(this->m_q1 - this->m_q0);
        auto const a = this->m_a;
        auto const v_lo = std::min(this->m_v0, this->m_v1);
        auto const v_hi = std::max(this->m_v0, this->m_v1);
        auto const v_sum_sq = this->m_v0 * this->m_v0 + this->m_v1 * this->m_v1;

        if(this->peak_below_boundary(h, a, v_sum_sq, v_hi))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        auto const residual_scale = std::max(h, v_sum_sq / (Scalar{2} * a));
        auto const ramp_residual = h - this->ramp_through_distance();
        auto const ramp_resolved = this->residual_resolved(ramp_residual, residual_scale);

        if(this->beyond_supremum(h, a, v_lo, v_hi, ramp_resolved, T_new))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        auto const at_hi = (v_hi > Scalar{0}) ? this->duration_and_scale_at(v_hi) : duration_with_scale<Scalar>{};
        auto const at_lo = (v_lo > Scalar{0}) ? this->duration_and_scale_at(v_lo) : duration_with_scale<Scalar>{};

        auto const selected = (v_hi <= Scalar{0} || T_new <= at_hi.T) ? rescale_shape::plateau
                : (v_lo <= Scalar{0} || T_new <= at_lo.T)             ? rescale_shape::ramp_through
                                                                      : rescale_shape::valley;

        // The two cruise-velocity forms stay in this function: a * h is shared with the peak test, and a compiler
        // that contracts across statements fuses it differently once it is used alone.
        Scalar v{};
        if(selected == rescale_shape::plateau && v_hi > Scalar{0})
        {
            auto const solved = this->plateau_rise_velocity(a, v_hi, residual_scale, ramp_residual, ramp_resolved, at_hi, T_new);
            if(!solved.has_value())
            {
                return ctrlpp::unexpected(solved.error());
            }
            v = *solved;
        }
        else if(selected == rescale_shape::plateau)
        {
            // A non-positive larger boundary velocity is no admissible cruise velocity, so the rise is never small
            // and the cruise velocity is solved directly, as the smaller root in the form that avoids a cancellation.
            auto const b = (this->m_v0 + this->m_v1) + a * T_new;
            auto const c = a * h + v_sum_sq / Scalar{2};
            auto const disc = b * b - Scalar{4} * c;
            if(!std::isfinite(disc))
            {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            if(!(disc >= Scalar{0}))
            {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            auto const root = std::sqrt(disc);
            v = (b > Scalar{0}) ? (Scalar{2} * c / (b + root)) : ((b - root) / Scalar{2});
        }
        else if(selected == rescale_shape::ramp_through)
        {
            auto const solved = this->ramp_through_velocity(a, ramp_residual, ramp_resolved, T_new);
            if(!solved.has_value())
            {
                return ctrlpp::unexpected(solved.error());
            }
            v = *solved;
        }
        else
        {
            // v = v_lo - d gives d^2 - B d + C = 0 with B = a (dT + r / v_lo) and C = a v_lo dT; the smaller
            // positive root is the decrement nearest the boundary.
            auto const dT = T_new - at_lo.T;
            if(!this->boundary_resolved(ramp_resolved, dT, at_lo.scale))
            {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            auto const B = a * (dT + ramp_residual / v_lo);
            auto const C = a * v_lo * dT;
            if(!(B > Scalar{0}))
            {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            auto const disc = B * B - Scalar{4} * C;
            // A non-finite discriminant drives the root onto the shape boundary, a success that realizes the
            // boundary's duration, so every solve reports it as unrepresentable rather than as reachability.
            if(!std::isfinite(disc))
            {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            if(!(disc >= Scalar{0}))
            {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            auto const d = Scalar{2} * C / (B + std::sqrt(disc));

            // Carry whichever of the decrement and the cruise velocity is the smaller: forming the smaller from
            // the boundary velocity and the larger is the cancellation this rearrangement exists to avoid. Deep in
            // the valley the cruise-velocity form's discriminant is a sum and carries no cancellation.
            if(d + d <= v_lo)
            {
                v = v_lo - d;
            }
            else
            {
                auto const b = a * T_new - (this->m_v0 + this->m_v1);
                auto const c = v_sum_sq / Scalar{2} - a * h;
                auto const disc_v = b * b - Scalar{4} * c;
                if(!std::isfinite(disc_v))
                {
                    return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
                }
                if(!(disc_v >= Scalar{0}))
                {
                    return ctrlpp::unexpected(trajectory_error::unreachable_duration);
                }
                auto const root = std::sqrt(disc_v);
                v = (b > Scalar{0}) ? (Scalar{-2} * c / (b + root)) : ((root - b) / Scalar{2});
            }
        }

        if(!velocity_in_shape(selected, v, v_lo, v_hi))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        return this->rescaled_state_at(h, a, v);
    }

    // The valley interval is open at its upper end: a root exactly on the boundary is what the closed form
    // returns when its discriminant cancels to zero, not a solution of the request.
    bool velocity_in_shape(rescale_shape selected, Scalar v, Scalar v_lo, Scalar v_hi) const
    {
        auto const in_shape = (selected == rescale_shape::plateau) ? (v >= v_hi) : (selected == rescale_shape::ramp_through) ? (v >= v_lo && v <= v_hi) : (v < v_lo);
        return std::isfinite(v) && v > Scalar{0} && !(v > this->m_v_max) && in_shape;
    }
};

}

#endif
