#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_2ND_PROFILE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_2ND_PROFILE_H

#include "ctrlpp/trajectory/detail/planner_2nd_state.h"

#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <cmath>
#include <limits>
#include <cstdint>

namespace ctrlpp::detail
{

template <typename Scalar>
class planner_2nd_profile : protected planner_2nd_state<Scalar>
{
protected:
    using planner_2nd_state<Scalar>::planner_2nd_state;

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2009, Sec. 3.2.7, p.73-74
    void compute_profile()
    {
        auto const h_signed = this->m_target - this->m_q_ref;
        // A command at the target is a no-op only from rest: removing any nonzero speed takes |v| / a_max, and the
        // planner is never told the sample period, so no magnitude exists below which dropping it respects a_max.
        if(h_signed == Scalar{0} && this->m_v_ref == Scalar{0})
        {
            this->settle_at_target();
            return;
        }

        this->m_settled = false;
        auto const v0 = this->m_v_ref;

        // One braking duration serves the overshoot decision and the stop placement, so the two cannot disagree.
        auto const brake_duration = std::abs(v0) / this->m_a_max;
        if(must_brake(h_signed, v0, brake_duration))
        {
            brake_and_replan(v0, brake_duration);
            return;
        }

        this->m_needs_brake = false;
        this->m_brake_duration = Scalar{0};
        compute_with_initial_velocity(this->m_q_ref, v0, this->m_target);
        this->record(online_planner_disposition::commanded_profile, online_planner_substitution_reason::none, this->m_duration, Scalar{0}, this->m_q_ref);
    }

    bool must_brake(Scalar h_signed, Scalar v0, Scalar brake_duration) const
    {
        // The mean-speed form, because it forms no square of the speed.
        auto const brake_dist = std::abs(v0) * brake_duration / Scalar{2};

        // One rounding at the scale of the velocity limit: no jerk is bounded, so the speed is tested as sampled, with
        // no acceleration-nulling arithmetic before it (the jerk-limited planner counts seven).
        constexpr std::int32_t speed_rounding_ops = 1;
        auto const speed_floor = Scalar{speed_rounding_ops} * std::numeric_limits<Scalar>::epsilon() * this->m_v_max;

        // At the scale of the stopping distance from full speed, v_max^2 / (2 a_max), which needs no sample period:
        // one rounding forms the displacement, and three the scale (the square, the doubling, the division).
        constexpr std::int32_t displacement_rounding_ops = 1;
        constexpr std::int32_t length_scale_rounding_ops = 3;
        constexpr std::int32_t length_rounding_ops = displacement_rounding_ops + length_scale_rounding_ops;
        auto const length_scale = this->m_v_max * this->m_v_max / (Scalar{2} * this->m_a_max);
        auto const length_floor = Scalar{length_rounding_ops} * std::numeric_limits<Scalar>::epsilon() * length_scale;

        bool const wrong_direction =
                (h_signed > Scalar{0} && v0 < Scalar{0}) || (h_signed < Scalar{0} && v0 > Scalar{0}) || (std::abs(h_signed) < length_floor && std::abs(v0) > speed_floor);
        bool const overshoot = !wrong_direction && brake_dist > std::abs(h_signed) && std::abs(v0) > speed_floor;
        return wrong_direction || overshoot;
    }

    void brake_and_replan(Scalar v0, Scalar brake_duration)
    {
        this->m_needs_brake = true;
        this->m_brake_duration = brake_duration;
        this->m_q_after_brake = this->m_q_ref + v0 * this->m_brake_duration / Scalar{2};
        compute_rest_to_rest(this->m_q_after_brake, this->m_target);
        this->record(online_planner_disposition::braked_and_replanned, online_planner_substitution_reason::reversal_or_overshoot, this->m_duration, this->m_brake_duration,
                this->m_q_after_brake);
    }

    void compute_rest_to_rest(Scalar q0, Scalar target)
    {
        auto const h = target - q0;
        this->m_sigma = (h >= Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const abs_h = std::abs(h);

        if(abs_h == Scalar{0})
        {
            this->set_standstill_after_brake();
            return;
        }
        plan_rest_to_rest(abs_h);
    }

    void plan_rest_to_rest(Scalar abs_h)
    {
        auto const v = this->m_v_max;
        auto const a = this->m_a_max;
        auto const v_tri_sq = a * abs_h;
        auto const v_tri = std::sqrt(v_tri_sq);
        select_peak(v_tri, v);

        this->m_a_a = a;
        this->m_a_d = a;
        this->m_accel_duration = this->m_v_v / this->m_a_a;
        this->m_decel_duration = this->m_v_v / this->m_a_d;
        if(v_tri >= v)
        {
            this->m_cruise_duration = abs_h / this->m_v_v - this->m_accel_duration / Scalar{2} - this->m_decel_duration / Scalar{2};
            if(this->m_cruise_duration < Scalar{0})
            {
                this->m_cruise_duration = Scalar{0};
            }
        }
        else
        {
            this->m_cruise_duration = Scalar{0};
        }
        this->m_duration = this->m_brake_duration + this->m_accel_duration + this->m_cruise_duration + this->m_decel_duration;
    }

    void compute_with_initial_velocity(Scalar q0, Scalar v0, Scalar target)
    {
        auto const h = target - q0;
        this->m_sigma = (h >= Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const abs_h = std::abs(h);
        // v0 in the positive frame; compute_profile brakes every speed pointing away from the target.
        auto const sv0 = this->m_sigma * v0;

        auto const a = this->m_a_max;
        auto const v = this->m_v_max;
        // The triangular peak, from h = (v_v^2 - sv0^2) / (2 a) + v_v^2 / (2 a).
        auto const v_tri_sq = (Scalar{2} * a * abs_h + sv0 * sv0) / Scalar{2};
        auto const v_tri = std::sqrt(v_tri_sq);
        select_peak(v_tri, v);
        // A peak below sv0 would need a deceleration first, which compute_profile handles by braking.
        if(this->m_v_v < sv0)
        {
            this->m_v_v = sv0;
        }

        this->m_a_a = a;
        this->m_a_d = a;
        this->m_accel_duration = (this->m_v_v - sv0) / this->m_a_a;
        this->m_decel_duration = this->m_v_v / this->m_a_d;
        set_carried_cruise(abs_h, sv0, a, v_tri, v);
        this->m_duration = this->m_accel_duration + this->m_cruise_duration + this->m_decel_duration;
    }

    void select_peak(Scalar v_tri, Scalar v)
    {
        if(v_tri < v)
        {
            this->m_v_v = v_tri;
        }
        else
        {
            this->m_v_v = v;
        }
    }

    void set_carried_cruise(Scalar abs_h, Scalar sv0, Scalar a, Scalar v_tri, Scalar v)
    {
        if(v_tri >= v)
        {
            auto const area_accel = sv0 * this->m_accel_duration + Scalar{0.5} * a * this->m_accel_duration * this->m_accel_duration;
            auto const area_decel = this->m_v_v * this->m_decel_duration - Scalar{0.5} * a * this->m_decel_duration * this->m_decel_duration;
            this->m_cruise_duration = (abs_h - area_accel - area_decel) / this->m_v_v;
            if(this->m_cruise_duration < Scalar{0})
            {
                this->m_cruise_duration = Scalar{0};
            }
        }
        else
        {
            this->m_cruise_duration = Scalar{0};
        }
    }
};

}

#endif
