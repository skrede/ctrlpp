#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_DIAGNOSTICS_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_DIAGNOSTICS_H

namespace ctrlpp
{

enum class online_planner_disposition
{
    commanded_profile,    // the commanded shape was planned as asked
    braked_and_replanned, // braked to rest, then replanned from the stopping point
    settled,              // already exactly at the target and at rest; the move became a zero-duration profile
};

// `online_planner_2nd` bounds no jerk, so it has no carry-velocity shape and never reports
// `carry_velocity_shape_unavailable`.
enum class online_planner_substitution_reason
{
    none,                             // nothing was substituted
    reversal_or_overshoot,            // the commanded motion reverses direction, or would overshoot the target
    carry_velocity_shape_unavailable, // the carry-velocity shape does not exist for the commanded state
};

template <typename Scalar>
struct online_planner_diagnostics
{
    online_planner_diagnostics()
        : disposition{online_planner_disposition::settled}
        , substitution_reason{online_planner_substitution_reason::none}
        , commanded_target{}
        , initial_velocity{}
        , planned_duration{}
        , brake_duration{}
        , replan_start_position{}
    {
    }

    online_planner_disposition disposition;
    online_planner_substitution_reason substitution_reason;
    Scalar commanded_target;
    Scalar initial_velocity;
    Scalar planned_duration;
    Scalar brake_duration;
    Scalar replan_start_position;
};

}

#endif
