#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_STATE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_3RD_STATE_H

#include "ctrlpp/trajectory/detail/jerk_limited_ramp.h"

#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
class planner_3rd_state
{
protected:
    Scalar m_v_max;
    Scalar m_a_max;
    Scalar m_j_max;
    // Read only by sample(), never a divisor and never a loop bound, so no value can make the planner run unboundedly.
    Scalar m_position_settle_tol;
    Scalar m_velocity_settle_tol;
    Scalar m_acceleration_settle_tol;
    mutable Scalar m_q;
    mutable Scalar m_v;
    mutable Scalar m_a;
    mutable Scalar m_t_last;
    mutable bool m_settled;
    Scalar m_q_ref;
    Scalar m_v_ref;
    Scalar m_a_ref;
    Scalar m_target;
    Scalar m_t_ref;
    // Every branch that completes a plan assigns every field through record(), so none survives from another branch.
    online_planner_diagnostics<Scalar> m_diagnostics;
    // At most nulling 1 + brake 3 + ramp 3 + cruise 1 + deceleration 3; an emptied group stores one pinned phase.
    // A phase whose end state is known analytically records it, and evaluation takes that state in place of the integrated one.
    // From the cruise end on, or over the last phase when there is no double-S, end positions are formed backward from the target.
    static constexpr std::int32_t max_phases = 11;
    std::int32_t m_phase_count;
    std::int32_t m_target_anchor_first;
    std::array<Scalar, max_phases> m_phase_durations;
    std::array<Scalar, max_phases> m_phase_jerks;
    std::array<Scalar, max_phases> m_phase_end_position;
    std::array<Scalar, max_phases> m_phase_end_velocity;
    std::array<Scalar, max_phases> m_phase_end_acceleration;
    std::array<bool, max_phases> m_phase_end_pinned;
    Scalar m_duration;

    template <typename Config>
    explicit planner_3rd_state(Config const& cfg)
        : m_v_max{cfg.v_max}
        , m_a_max{cfg.a_max}
        , m_j_max{cfg.j_max}
        , m_position_settle_tol{cfg.position_settle_tol}
        , m_velocity_settle_tol{cfg.velocity_settle_tol}
        , m_acceleration_settle_tol{cfg.acceleration_settle_tol}
        , m_q{}
        , m_v{}
        , m_a{}
        , m_t_last{}
        , m_settled{true}
        , m_q_ref{}
        , m_v_ref{}
        , m_a_ref{}
        , m_target{}
        , m_t_ref{}
        , m_diagnostics{}
        , m_phase_count{0}
        , m_target_anchor_first{0}
        , m_phase_durations{}
        , m_phase_jerks{}
        , m_phase_end_position{}
        , m_phase_end_velocity{}
        , m_phase_end_acceleration{}
        , m_phase_end_pinned{}
        , m_duration{}
    {
    }

    void record(online_planner_disposition disposition, online_planner_substitution_reason reason, Scalar planned_duration, Scalar brake_duration, Scalar replan_start_position)
    {
        m_diagnostics.disposition = disposition;
        m_diagnostics.substitution_reason = reason;
        m_diagnostics.commanded_target = m_target;
        m_diagnostics.initial_velocity = m_v_ref;
        m_diagnostics.planned_duration = planned_duration;
        m_diagnostics.brake_duration = brake_duration;
        m_diagnostics.replan_start_position = replan_start_position;
    }

    void hold_target(Scalar& q, Scalar& v, Scalar& a) const
    {
        q = m_target;
        v = Scalar{0};
        a = Scalar{0};
    }

    void settle_at_target()
    {
        m_duration = Scalar{0};
        m_settled = true;
        record(online_planner_disposition::settled, online_planner_substitution_reason::none, Scalar{0}, Scalar{0}, m_q_ref);
    }

    static Scalar travel(Scalar v, Scalar a, Scalar j, Scalar t) { return v * t + Scalar{0.5} * a * t * t + j * t * t * t / Scalar{6}; }

    // A jerk phase starts or pinned ends at zero acceleration, so bounding its duration, which analytically never
    // exceeds a_max / j_max, bounds every acceleration.
    void append_phase(Scalar duration, Scalar jerk)
    {
        if(duration > Scalar{0})
        {
            auto const index = static_cast<std::size_t>(m_phase_count);
            m_phase_durations[index] = (jerk == Scalar{0}) ? duration : std::min(duration, jerk_time_within(m_a_max, m_j_max));
            m_phase_jerks[index] = jerk;
            ++m_phase_count;
        }
    }

    // Stores a zero-length phase to carry the end state when the group since `first` stored none.
    void pin_group_end(std::int32_t first, Scalar v, Scalar a)
    {
        if(m_phase_count == first)
        {
            m_phase_durations[static_cast<std::size_t>(m_phase_count)] = Scalar{0};
            m_phase_jerks[static_cast<std::size_t>(m_phase_count)] = Scalar{0};
            ++m_phase_count;
        }
        auto const last = static_cast<std::size_t>(m_phase_count - 1);
        m_phase_end_velocity[last] = v;
        m_phase_end_acceleration[last] = a;
        m_phase_end_pinned[last] = true;
    }

    // Walks back from (target, 0, 0) through the pinned or back-integrated end states, so on the final approach no
    // forward sum of earlier roundings can carry a position past the target.
    void anchor_target_positions()
    {
        m_target_anchor_first = std::max(std::min(m_target_anchor_first, m_phase_count - 1), std::int32_t{0});
        auto q = m_target;
        Scalar v{};
        Scalar a{};
        for(auto i = m_phase_count - 1; i >= m_target_anchor_first; --i)
        {
            auto const index = static_cast<std::size_t>(i);
            if(m_phase_end_pinned[index])
            {
                v = m_phase_end_velocity[index];
                a = m_phase_end_acceleration[index];
            }
            m_phase_end_position[index] = q;
            auto const T = m_phase_durations[index];
            auto const j = m_phase_jerks[index];
            q -= travel(v, -a, j, T);
            v = v - a * T + Scalar{0.5} * j * T * T;
            a = a - j * T;
        }
    }
};

}

#endif
