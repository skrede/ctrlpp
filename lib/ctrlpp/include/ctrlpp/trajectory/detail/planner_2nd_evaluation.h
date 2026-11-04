#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_2ND_EVALUATION_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_2ND_EVALUATION_H

#include "ctrlpp/trajectory/detail/ramp_anchor.h"
#include "ctrlpp/trajectory/detail/planner_2nd_profile.h"

namespace ctrlpp::detail
{

template <typename Scalar>
class planner_2nd_evaluation : protected planner_2nd_profile<Scalar>
{
protected:
    using planner_2nd_profile<Scalar>::planner_2nd_profile;

    // @cite biagiotti2009 -- Sec. 4.6.2
    void evaluate_profile(Scalar dt, Scalar& q, Scalar& v, Scalar& a) const
    {
        if(this->m_duration <= Scalar{0} || dt >= this->m_duration)
        {
            this->hold_target(q, v, a);
            return;
        }

        if(this->m_needs_brake)
        {
            if(dt < this->m_brake_duration)
            {
                auto const brake_sign = (this->m_v_ref >= Scalar{0}) ? Scalar{-1} : Scalar{1};
                a = brake_sign * this->m_a_max;
                v = this->m_v_ref + a * dt;
                q = this->m_q_ref + this->m_v_ref * dt + Scalar{0.5} * a * dt * dt;
                return;
            }
            dt -= this->m_brake_duration;
            evaluate_rest_to_rest(dt, this->m_q_after_brake, q, v, a);
            return;
        }

        evaluate_with_velocity(dt, q, v, a);
    }

    void evaluate_rest_to_rest(Scalar dt, Scalar q0, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const total_duration = this->m_accel_duration + this->m_cruise_duration + this->m_decel_duration;
        if(dt >= total_duration)
        {
            this->hold_target(q, v, a);
            return;
        }

        auto const dt_end = total_duration - dt;
        if(dt < this->m_accel_duration)
        {
            accelerate_from_rest(dt, q0, q, v, a);
        }
        else if(dt_end >= this->m_decel_duration)
        {
            cruise_from_rest(dt, q0, q, v, a);
        }
        else
        {
            decelerate(dt_end, q, v, a);
        }
    }

    void accelerate_from_rest(Scalar dt, Scalar q0, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const s = this->m_sigma;
        auto const aa = this->m_a_a;
        a = s * aa;
        v = s * aa * dt;
        q = q0 + s * Scalar{0.5} * aa * dt * dt;
    }

    void cruise_from_rest(Scalar dt, Scalar q0, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const s = this->m_sigma;
        auto const aa = this->m_a_a;
        auto const dt_c = dt - this->m_accel_duration;
        auto const q_end_accel = q0 + s * Scalar{0.5} * aa * this->m_accel_duration * this->m_accel_duration;
        a = Scalar{0};
        v = s * this->m_v_v;
        q = q_end_accel + s * this->m_v_v * dt_c;
    }

    void evaluate_with_velocity(Scalar dt, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const s = this->m_sigma;
        auto const aa = this->m_a_a;
        auto const total_duration = this->m_accel_duration + this->m_cruise_duration + this->m_decel_duration;

        if(dt >= total_duration)
        {
            this->hold_target(q, v, a);
            return;
        }
        auto const dt_end = total_duration - dt;
        auto const q_end_accel = this->m_q_ref + this->m_v_ref * this->m_accel_duration + s * Scalar{0.5} * aa * this->m_accel_duration * this->m_accel_duration;

        if(dt < this->m_accel_duration)
        {
            accelerate_carried(dt, q_end_accel, q, v, a);
        }
        else if(dt_end >= this->m_decel_duration)
        {
            cruise_carried(dt, q_end_accel, q, v, a);
        }
        else
        {
            decelerate(dt_end, q, v, a);
        }
    }

    void accelerate_carried(Scalar dt, Scalar q_end_accel, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const s = this->m_sigma;
        auto const aa = this->m_a_a;
        a = s * aa;
        v = detail::nearer_anchor_velocity(this->m_v_ref, s * this->m_v_v, s * aa, dt, this->m_accel_duration - dt);
        q = detail::nearer_anchor_position(this->m_q_ref, this->m_v_ref, q_end_accel, s * this->m_v_v, s * aa, dt, this->m_accel_duration - dt);
    }

    void cruise_carried(Scalar dt, Scalar q_end_accel, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const s = this->m_sigma;
        a = Scalar{0};
        v = s * this->m_v_v;
        q = q_end_accel + s * this->m_v_v * (dt - this->m_accel_duration);
    }

    // The total is a rounded sum, so deceleration is used only strictly inside its phase, timed back from the end.
    void decelerate(Scalar dt_end, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const s = this->m_sigma;
        a = -s * this->m_a_d;
        v = s * this->m_a_d * dt_end;
        q = this->m_target - s * Scalar{0.5} * this->m_a_d * dt_end * dt_end;
    }
};

}

#endif
