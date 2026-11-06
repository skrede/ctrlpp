#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_PROFILE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_PROFILE_H

#include "ctrlpp/trajectory/detail/jerk_limited_ramp.h"
#include "ctrlpp/trajectory/detail/planner_3rd_double_s.h"

#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
class planner_3rd_profile : protected planner_3rd_double_s<Scalar>
{
protected:
    struct substitution_outcome
    {
        substitution_outcome()
            : reason{online_planner_substitution_reason::none}
            , brake_duration{}
            , replan_start_position{}
        {
        }

        online_planner_substitution_reason reason;
        Scalar brake_duration;
        Scalar replan_start_position;
    };

    // Upper bound on the roundings in the nulling-end speed: the duration's division, the product a T, the three products of
    // the jerk term (its exact halving counted anyway) and the two sums. A plan straight from the snapshot performs none.
    static constexpr std::int32_t speed_rounding_ops = 7;

    using planner_3rd_double_s<Scalar>::planner_3rd_double_s;

    // Variable-structure approach: a nonzero acceleration is first nulled by a jerk phase, and a double-S is then
    // planned from the resulting state at zero acceleration.
    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2009, Sec. 3.4, p.79-93
    void compute_profile()
    {
        this->m_phase_count = 0;
        this->m_target_anchor_first = this->max_phases;
        this->m_phase_end_pinned.fill(false);

        // A command at the target is a no-op only from rest: removing any nonzero velocity or acceleration takes time
        // at the limits, and no magnitude exists below which dropping it would respect them.
        if(this->m_target == this->m_q_ref && this->m_v_ref == Scalar{0} && this->m_a_ref == Scalar{0})
        {
            this->settle_at_target();
            return;
        }

        this->m_settled = false;
        plan_from_snapshot();
    }

    void plan_from_snapshot()
    {
        Scalar q_start = this->m_q_ref;
        Scalar v_start = this->m_v_ref;
        if(this->m_a_ref != Scalar{0})
        {
            append_nulling_phase(q_start, v_start);
        }

        auto const outcome = plan_from_zero_accel(q_start, v_start);
        this->anchor_target_positions();

        this->m_duration = Scalar{0};
        for(std::int32_t i = 0; i < this->m_phase_count; ++i)
        {
            this->m_duration += this->m_phase_durations[static_cast<std::size_t>(i)];
        }
        if(outcome.reason == online_planner_substitution_reason::none)
        {
            this->record(online_planner_disposition::commanded_profile, online_planner_substitution_reason::none, this->m_duration, Scalar{0}, this->m_q_ref);
        }
        else
        {
            this->record(online_planner_disposition::braked_and_replanned, outcome.reason, this->m_duration, outcome.brake_duration, outcome.replan_start_position);
        }
    }

    // The evaluation starts from the snapshot's acceleration, so every nonzero one is nulled by its own phase,
    // however short. A phase too short to represent still pins its end state.
    void append_nulling_phase(Scalar& q_start, Scalar& v_start)
    {
        auto const a_start = this->m_a_ref;
        auto const T_az = std::min(std::abs(a_start) / this->m_j_max, jerk_time_within(this->m_a_max, this->m_j_max));
        auto const j_az = (a_start > Scalar{0}) ? -this->m_j_max : this->m_j_max;
        auto const v_after = snap_nulling_end_speed(v_start + a_start * T_az + Scalar{0.5} * j_az * T_az * T_az);
        q_start = q_start + v_start * T_az + Scalar{0.5} * a_start * T_az * T_az + j_az * T_az * T_az * T_az / Scalar{6};
        auto const first = this->m_phase_count;
        this->append_phase(T_az, j_az);
        this->pin_group_end(first, v_after, Scalar{0});
        v_start = v_after;
    }

    // Analytically the nulling-end speed of a state this planner produced is at most the ramp's peak, so an
    // excess within the counted roundings is that rounding and is placed on the limit; a larger excess is
    // left alone for the brake path.
    Scalar snap_nulling_end_speed(Scalar v_after) const
    {
        auto const excess = std::abs(v_after) - this->m_v_max;
        auto const rounding = static_cast<Scalar>(speed_rounding_ops) * std::numeric_limits<Scalar>::epsilon() * this->m_v_max;
        if(excess > Scalar{0} && excess <= rounding)
        {
            return std::copysign(this->m_v_max, v_after);
        }
        return v_after;
    }

    substitution_outcome plan_from_zero_accel(Scalar q0, Scalar v0)
    {
        auto const h_signed = this->m_target - q0;

        // The evaluation integrates from any nonzero speed, so only an exactly zero one is planned as rest.
        if(v0 == Scalar{0})
        {
            this->plan_rest_to_rest(q0);
            return {};
        }

        auto reason = online_planner_substitution_reason::reversal_or_overshoot;
        if(!must_brake(h_signed, v0))
        {
            if(this->append_incorporate_velocity(q0, v0))
            {
                return {};
            }
            // With the shape's own displacement test passed and the carried speed exact, the shape refuses only through
            // an arithmetic rejection; braking to rest and replanning is then the always-admissible fallback.
            reason = online_planner_substitution_reason::carry_velocity_shape_unavailable;
        }
        return brake_and_replan(q0, v0, reason);
    }

    bool must_brake(Scalar h_signed, Scalar v0) const
    {
        // Each floor has the scale of the quantity it tests and a count of the roundings that formed it, so a length
        // and a speed never meet inside one expression. Below the speed floor a speed's sign says nothing about direction.
        auto const speed_floor = static_cast<Scalar>(speed_rounding_ops) * std::numeric_limits<Scalar>::epsilon() * this->m_v_max;

        // The length scale is the stopping distance from full speed, intrinsic to the limits and independent of the
        // sample period, which the planner is never told. The nulling phase's position takes twelve roundings (the
        // duration's division, two in the speed term, four in the squared term, five in the cubed term), the
        // displacement one and the scale three; a plan that nulled no acceleration performs fewer.
        constexpr std::int32_t nulling_position_rounding_ops = 12;
        constexpr std::int32_t displacement_rounding_ops = 1;
        constexpr std::int32_t length_scale_rounding_ops = 3;
        constexpr std::int32_t length_rounding_ops = nulling_position_rounding_ops + displacement_rounding_ops + length_scale_rounding_ops;
        auto const length_scale = this->m_v_max * this->m_v_max / (Scalar{2} * this->m_a_max);
        auto const length_floor = Scalar{length_rounding_ops} * std::numeric_limits<Scalar>::epsilon() * length_scale;

        bool const wrong_way = (h_signed > length_floor && v0 < -speed_floor) || (h_signed < -length_floor && v0 > speed_floor) || (std::abs(h_signed) < length_floor);

        // The carry shape's own displacement test on the operands create() sees, so the verdict and the
        // shape's refusal agree bit for bit.
        bool const overshoot = !wrong_way && std::abs(h_signed) < no_cruise_displacement(Scalar{0}, std::abs(v0), Scalar{0}, this->m_a_max, this->m_j_max);
        return wrong_way || overshoot;
    }

    substitution_outcome brake_and_replan(Scalar q0, Scalar v0, online_planner_substitution_reason reason)
    {
        auto const stop_info = this->compute_stop(v0);
        this->append_brake_phases(v0, stop_info);
        auto const q_after = q0 + stop_info.dist;
        this->plan_rest_to_rest(q_after);
        substitution_outcome outcome;
        outcome.reason = reason;
        outcome.brake_duration = stop_info.brake_duration;
        outcome.replan_start_position = q_after;
        return outcome;
    }
};

}

#endif
