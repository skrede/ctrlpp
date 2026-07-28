#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_DIAGNOSTICS_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_DIAGNOSTICS_H

namespace ctrlpp
{

/// Which profile the last update built, against the one it was commanded.
enum class online_planner_disposition
{
    commanded_profile,     ///< the commanded shape was planned as asked
    braked_and_replanned,  ///< braked to rest, then replanned from the stopping point
    settled,               ///< already within the settle tolerance; the move became a zero-duration profile
};

/// Why the commanded profile shape was replaced by a brake-then-replan.
///
/// `carry_velocity_shape_unavailable` is unreachable from `online_planner_2nd`:
/// that planner bounds no jerk, so it has no carry-velocity shape whose domain
/// the commanded state can fall outside of. It reports `none` or
/// `reversal_or_overshoot` and nothing else.
enum class online_planner_substitution_reason
{
    none,                              ///< nothing was substituted
    reversal_or_overshoot,             ///< the commanded motion reverses direction, or would overshoot the target
    carry_velocity_shape_unavailable,  ///< the carry-velocity shape does not exist for the commanded state
};

/// Disposition of the last update: what was commanded, against what was planned.
template <typename Scalar>
struct [[nodiscard]] online_planner_diagnostics
{
    online_planner_disposition          disposition{online_planner_disposition::settled};
    online_planner_substitution_reason  substitution_reason{online_planner_substitution_reason::none};
    Scalar commanded_target{};       ///< target position the update was given
    Scalar initial_velocity{};       ///< velocity the planner was carrying when it planned
    Scalar planned_duration{};       ///< total duration the plan realizes
    Scalar brake_duration{};         ///< duration spent braking before the replan; zero when nothing was substituted
    Scalar replan_start_position{};  ///< position the replan starts from; the commanded start position when nothing was substituted
};

}

#endif
