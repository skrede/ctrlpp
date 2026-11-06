#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_EVALUATION_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_EVALUATION_H

#include "ctrlpp/trajectory/detail/planner_3rd_profile.h"

#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

template <typename Scalar>
class planner_3rd_evaluation : protected planner_3rd_profile<Scalar>
{
protected:
    using planner_3rd_profile<Scalar>::planner_3rd_profile;

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.4.1, eq. (3.30a)-(3.30g), p.88-89
    void evaluate_profile(Scalar dt, Scalar& q, Scalar& v, Scalar& a) const
    {
        if(this->m_phase_count == 0 || dt >= this->m_duration)
        {
            this->hold_target(q, v, a);
            return;
        }

        q = this->m_q_ref;
        v = this->m_v_ref;
        a = this->m_a_ref;
        Scalar t_elapsed = Scalar{0};
        for(std::int32_t i = 0; i < this->m_phase_count; ++i)
        {
            auto const T_i = this->m_phase_durations[static_cast<std::size_t>(i)];
            auto const t_remaining = dt - t_elapsed;
            if(t_remaining <= Scalar{0})
                break;
            if(t_remaining < T_i)
            {
                advance_phase(i, t_remaining, q, v, a);
                break;
            }
            advance_phase(i, T_i, q, v, a);
            t_elapsed += T_i;
        }
    }

    // Inside a pinned phase the sample is formed backward from the pinned end state, because a forward sum can round past the end value just before the boundary.
    // A pinned zero-jerk phase is sampled from its nearer end, since near its start the pin's magnitude would quantize a much smaller velocity.
    // Positions take the backward form only where they are anchored to the target; before that they are summed forward.
    void advance_phase(std::int32_t i, Scalar tau, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const index = static_cast<std::size_t>(i);
        auto const j_i = this->m_phase_jerks[index];
        auto const s = this->m_phase_durations[index] - tau;
        bool const forward = !this->m_phase_end_pinned[index] || (j_i == Scalar{0} && tau < s);
        if(forward || i < this->m_target_anchor_first)
        {
            q += this->travel(v, a, j_i, tau);
        }
        else
        {
            q = this->m_phase_end_position[index] - this->travel(this->m_phase_end_velocity[index], -this->m_phase_end_acceleration[index], j_i, s);
        }
        if(forward)
        {
            v += a * tau + Scalar{0.5} * j_i * tau * tau;
            a += j_i * tau;
            return;
        }
        v = this->m_phase_end_velocity[index] - this->m_phase_end_acceleration[index] * s + Scalar{0.5} * j_i * s * s;
        a = this->m_phase_end_acceleration[index] - j_i * s;
    }
};

}

#endif
