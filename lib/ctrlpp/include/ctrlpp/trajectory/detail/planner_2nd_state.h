#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_2ND_STATE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PLANNER_2ND_STATE_H

#include "ctrlpp/trajectory/online_planner_diagnostics.h"

namespace ctrlpp::detail
{

template <typename Scalar>
class planner_2nd_state
{
protected:
    Scalar m_v_max;
    Scalar m_a_max;
    // Read only by sample(), never a divisor and never a loop bound, so no value can make the planner run unboundedly.
    Scalar m_position_settle_tol;
    Scalar m_velocity_settle_tol;
    mutable Scalar m_q;
    mutable Scalar m_v;
    mutable Scalar m_t_last;
    mutable bool m_settled;
    Scalar m_q_ref;
    Scalar m_v_ref;
    Scalar m_target;
    Scalar m_t_ref;
    // Every branch that completes a plan assigns every field through record(), so none survives from another branch.
    online_planner_diagnostics<Scalar> m_diagnostics;
    Scalar m_sigma;
    Scalar m_v_v; // cruise velocity, positive frame
    Scalar m_a_a; // acceleration magnitude
    Scalar m_a_d; // deceleration magnitude
    Scalar m_accel_duration;
    Scalar m_cruise_duration;
    Scalar m_decel_duration;
    Scalar m_duration;
    Scalar m_q_stop;
    bool m_needs_brake;
    Scalar m_brake_duration;
    Scalar m_q_after_brake;
    Scalar m_v_after_brake;

    template <typename Config>
    explicit planner_2nd_state(Config const& cfg)
        : m_v_max{cfg.v_max}
        , m_a_max{cfg.a_max}
        , m_position_settle_tol{cfg.position_settle_tol}
        , m_velocity_settle_tol{cfg.velocity_settle_tol}
        , m_q{}
        , m_v{}
        , m_t_last{}
        , m_settled{true}
        , m_q_ref{}
        , m_v_ref{}
        , m_target{}
        , m_t_ref{}
        , m_diagnostics{}
        , m_sigma{1}
        , m_v_v{}
        , m_a_a{}
        , m_a_d{}
        , m_accel_duration{}
        , m_cruise_duration{}
        , m_decel_duration{}
        , m_duration{}
        , m_q_stop{}
        , m_needs_brake{false}
        , m_brake_duration{}
        , m_q_after_brake{}
        , m_v_after_brake{}
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

    void settle_at_target()
    {
        m_sigma = Scalar{1};
        m_v_v = Scalar{0};
        m_a_a = Scalar{0};
        m_a_d = Scalar{0};
        m_accel_duration = Scalar{0};
        m_cruise_duration = Scalar{0};
        m_decel_duration = Scalar{0};
        m_duration = Scalar{0};
        m_needs_brake = false;
        m_settled = true;
        record(online_planner_disposition::settled, online_planner_substitution_reason::none, Scalar{0}, Scalar{0}, m_q_ref);
    }

    void set_standstill_after_brake()
    {
        m_v_v = Scalar{0};
        m_a_a = Scalar{0};
        m_a_d = Scalar{0};
        m_accel_duration = Scalar{0};
        m_cruise_duration = Scalar{0};
        m_decel_duration = Scalar{0};
        m_duration = m_brake_duration;
    }

    void hold_target(Scalar& q, Scalar& v, Scalar& a) const
    {
        q = m_target;
        v = Scalar{0};
        a = Scalar{0};
    }
};

}

#endif
