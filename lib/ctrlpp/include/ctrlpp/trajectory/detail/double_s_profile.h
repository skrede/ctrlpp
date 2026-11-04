#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_PROFILE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_PROFILE_H

#include "ctrlpp/trajectory/detail/ramp_anchor.h"
#include "ctrlpp/trajectory/detail/double_s_solve.h"

#include "ctrlpp/trajectory/trajectory_types.h"

#include <cmath>

namespace ctrlpp::detail
{

template <typename Scalar>
class double_s_profile : protected double_s_solve<Scalar>
{
protected:
    using double_s_solve<Scalar>::double_s_solve;

    trajectory_point<Scalar, 1> make_point(Scalar q, Scalar dq, Scalar ddq) const
    {
        return {.position = Vector<Scalar, 1>{q}, .velocity = Vector<Scalar, 1>{dq}, .acceleration = Vector<Scalar, 1>{ddq}};
    }

    static void store_point(Scalar& q, Scalar& dq, Scalar& ddq, Scalar q_value, Scalar dq_value, Scalar ddq_value)
    {
        q = q_value;
        dq = dq_value;
        ddq = ddq_value;
    }

    // @cite biagiotti2009 -- Sec. 3.4, eq. (3.30a)-(3.30g), p.85-86
    void evaluate_positive_frame(Scalar tc, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        auto const t3 = this->m_accel_duration;
        auto const t4 = this->m_accel_duration + this->m_cruise_duration;

        if(tc < t3)
        {
            eval_accel_phase(tc, q, dq, ddq);
        }
        else if(tc < t4)
        {
            auto const q_at_Ta = compute_q_at_Ta();
            auto const t = tc - this->m_accel_duration;
            q = q_at_Ta + this->m_v_lim * t;
            dq = this->m_v_lim;
            ddq = Scalar{0};
        }
        else
        {
            eval_decel_phase(tc, q, dq, ddq);
        }
    }

    // @cite biagiotti2009 -- Sec. 3.4, eq. (3.30a)-(3.30c), p.85
    void eval_accel_phase(Scalar tc, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        auto const j = this->m_j_max;
        auto const dt = this->m_accel_duration - tc;

        // The stored switch times are rounded sums, so each jerk segment's formula is used only strictly
        // within that segment, timed from the end it is anchored at; between them the stored peak is held.
        if(tc < this->m_accel_jerk_duration)
        {
            q = this->m_pv0 * tc + j * tc * tc * tc / Scalar{6};
            dq = this->m_pv0 + j * tc * tc / Scalar{2};
            ddq = j * tc;
        }
        else if(dt < this->m_accel_jerk_duration)
        {
            q = compute_q_at_Ta() - this->m_v_lim * dt + j * dt * dt * dt / Scalar{6};
            dq = this->m_v_lim - j * dt * dt / Scalar{2};
            ddq = j * dt;
        }
        else
        {
            eval_constant_accel(tc - this->m_accel_jerk_duration, dt - this->m_accel_jerk_duration, q, dq);
            ddq = this->m_a_lim_a;
        }
    }

    void eval_constant_accel(Scalar t_from_start, Scalar t_to_end, Scalar& q, Scalar& dq) const
    {
        auto const j = this->m_j_max;
        auto const q_start = this->m_pv0 * this->m_accel_jerk_duration + j * this->m_accel_jerk_duration * this->m_accel_jerk_duration * this->m_accel_jerk_duration / Scalar{6};
        auto const v_start = this->m_pv0 + j * this->m_accel_jerk_duration * this->m_accel_jerk_duration / Scalar{2};
        auto const q_end = compute_q_at_Ta() - this->m_v_lim * this->m_accel_jerk_duration +
                j * this->m_accel_jerk_duration * this->m_accel_jerk_duration * this->m_accel_jerk_duration / Scalar{6};
        auto const v_end = this->m_v_lim - j * this->m_accel_jerk_duration * this->m_accel_jerk_duration / Scalar{2};
        q = nearer_anchor_position(q_start, v_start, q_end, v_end, this->m_a_lim_a, t_from_start, t_to_end);
        dq = nearer_anchor_velocity(v_start, v_end, this->m_a_lim_a, t_from_start, t_to_end);
    }

    // @cite biagiotti2009 -- Sec. 3.4, eq. (3.30e)-(3.30g), p.86
    void eval_decel_phase(Scalar tc, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        auto const j = this->m_j_max;
        auto const dt_end = this->m_duration - tc;
        auto const t = tc - (this->m_accel_duration + this->m_cruise_duration);

        if(t < this->m_decel_jerk_duration)
        {
            q = compute_q_at_Ta() + this->m_v_lim * this->m_cruise_duration + this->m_v_lim * t - j * t * t * t / Scalar{6};
            dq = this->m_v_lim - j * t * t / Scalar{2};
            ddq = -j * t;
        }
        else if(dt_end < this->m_decel_jerk_duration)
        {
            q = std::abs(this->m_q1 - this->m_q0) - this->m_pv1 * dt_end - j * dt_end * dt_end * dt_end / Scalar{6};
            dq = this->m_pv1 + j * dt_end * dt_end / Scalar{2};
            ddq = -j * dt_end;
        }
        else
        {
            // Timed from the selection tests' own differences, so neither is negative; tc - t5 could be, since t5
            // is a separately rounded sum.
            eval_constant_decel(t - this->m_decel_jerk_duration, dt_end - this->m_decel_jerk_duration, q, dq);
            ddq = -this->m_a_lim_d;
        }
    }

    void eval_constant_decel(Scalar t_from_start, Scalar t_to_end, Scalar& q, Scalar& dq) const
    {
        auto const j = this->m_j_max;
        auto const q_start = compute_q_at_Ta() + this->m_v_lim * this->m_cruise_duration + this->m_v_lim * this->m_decel_jerk_duration -
                j * this->m_decel_jerk_duration * this->m_decel_jerk_duration * this->m_decel_jerk_duration / Scalar{6};
        auto const v_start = this->m_v_lim - j * this->m_decel_jerk_duration * this->m_decel_jerk_duration / Scalar{2};
        auto const q_end = std::abs(this->m_q1 - this->m_q0) - this->m_pv1 * this->m_decel_jerk_duration -
                j * this->m_decel_jerk_duration * this->m_decel_jerk_duration * this->m_decel_jerk_duration / Scalar{6};
        auto const v_end = this->m_pv1 + j * this->m_decel_jerk_duration * this->m_decel_jerk_duration / Scalar{2};
        q = nearer_anchor_position(q_start, v_start, q_end, v_end, -this->m_a_lim_d, t_from_start, t_to_end);
        dq = nearer_anchor_velocity(v_start, v_end, -this->m_a_lim_d, t_from_start, t_to_end);
    }

    // The acceleration profile is symmetric about the ramp's midpoint, so the velocity ramp is point-symmetric
    // and the trapezoidal area is exact for any initial velocity.
    Scalar compute_q_at_Ta() const { return (this->m_pv0 + this->m_v_lim) * this->m_accel_duration / Scalar{2}; }
};

}

#endif
