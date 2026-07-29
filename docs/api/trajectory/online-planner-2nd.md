# online_planner_2nd

2nd-order online trajectory planner that generates trapezoidal-like velocity profiles in real time. On each `update(target)`, the planner computes a time-optimal profile from the current state to the target position, respecting `v_max` and `a_max` constraints. Suitable for real-time control loops where targets change dynamically.

Unlike pre-computed trajectory segments, online planners are stateful filters with no fixed duration and do not satisfy `trajectory_segment`.

## Header

```cpp
#include "ctrlpp/trajectory/online_planner_2nd.h"
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (`float`, `double`, `long double`) |

## Config

```cpp
struct config {
    Scalar v_max;   // Maximum velocity magnitude
    Scalar a_max;   // Maximum acceleration magnitude
};
```

Both limits divide in the planner math (stopping distance `v^2 / (2 * a_max)`, phase durations `v_v / a_max` and `h / v_v`), so the domain of each is finite and strictly positive.

## Construction

```cpp
static auto create(config const& cfg)
    -> ctrlpp::expected<online_planner_2nd, trajectory_error>;
```

`create` is the only construction path. There is no non-fallible constructor: a rejected configuration is a value the caller has to inspect, never an object that quietly stands in for one. It validates the kinematic limits and constructs a planner with initial state at rest at q = 0. Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| NaN/Inf or non-positive `v_max` | `trajectory_error::non_positive_velocity_limit` |
| NaN/Inf or non-positive `a_max` | `trajectory_error::non_positive_acceleration_limit` |

## Methods

### update

```cpp
auto update(Scalar target) -> ctrlpp::expected<void, trajectory_error>;
```

Set a new target position and replan from the current state. Computes a time-optimal trapezoidal profile from (q, v) to (target, 0) respecting `v_max` and `a_max`. A velocity pointing away from the target, or too large to stop in the available distance, is braked to rest first and the move is replanned from the stopping point.

The command is rejected with `trajectory_error::non_finite_input` before any
planner state changes when the target is NaN or infinite. For an accepted
target, which profile was built is read back from
[`diagnostics()`](#diagnostics). The motion respects every limit either way;
what changes is the time it takes.

### diagnostics

```cpp
auto diagnostics() const -> online_planner_diagnostics<Scalar> const&;
```

Report what the last `update` (or `reset`) planned, against what it was commanded. This describes the plan, not the state reached since; `is_settled()` answers that.

### sample

```cpp
auto sample(Scalar t) const -> trajectory_point<Scalar, 1>;
```

Evaluate the planned trajectory at absolute time `t`. Returns position, velocity, and acceleration. Updates internal state for future `update()` calls. Call this once per control loop timestep.

### is_settled

```cpp
auto is_settled() const -> bool;
```

Returns true when the planner has reached the target with zero velocity.

### reset

```cpp
void reset(Scalar q0);
```

Reset state to position `q0` with zero velocity.

## Profile Behavior

The planner produces trapezoidal velocity profiles with three phases:

1. **Acceleration**<br/>ramp velocity toward `v_max` at rate `a_max`
2. **Cruise**<br/>hold at `v_max` (may be zero duration for short moves)
3. **Deceleration**<br/>ramp velocity to zero at rate `a_max`

For short displacements where `v_max` cannot be reached, the profile degenerates to a triangular velocity shape. Mid-motion target changes trigger replanning from the current state, with automatic brake-and-replan for reversal and overshoot scenarios.

## Substitution Reporting

The brake-and-replan is a substitution: the caller commanded a move from the current velocity, and the planner realized a different, longer profile that respects the same limits and reaches the same target. The motion alone does not distinguish the two. The disposition does.

```cpp
enum class online_planner_disposition
{
    commanded_profile,     // the commanded shape was planned as asked
    braked_and_replanned,  // braked to rest, then replanned from the stopping point
    settled,               // already within the settle tolerance; a zero-duration profile
};

enum class online_planner_substitution_reason
{
    none,
    reversal_or_overshoot,             // the commanded motion reverses, or would overshoot
    carry_velocity_shape_unavailable,  // reported only by the 3rd-order planner
};

template <typename Scalar>
struct [[nodiscard]] online_planner_diagnostics
{
    online_planner_disposition          disposition;
    online_planner_substitution_reason  substitution_reason;
    Scalar commanded_target;
    Scalar initial_velocity;
    Scalar planned_duration;
    Scalar brake_duration;
    Scalar replan_start_position;
};
```

| Field | Meaning |
|-------|---------|
| `disposition` | Which profile was built |
| `substitution_reason` | Which condition selected the substitution; `none` when nothing was substituted |
| `commanded_target` | Target position the update was given |
| `initial_velocity` | Velocity the planner was carrying when it planned |
| `planned_duration` | Total duration the plan realizes |
| `brake_duration` | Duration spent braking before the replan; zero when nothing was substituted |
| `replan_start_position` | Position the replan starts from; the commanded start position when nothing was substituted |

The type is shared with [online-planner-3rd](online-planner-3rd.md). This planner bounds no jerk, so it has no carry-velocity shape whose domain a commanded state can fall outside of: it reports `none` or `reversal_or_overshoot` and nothing else.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include <ctrlpp/trajectory/online_planner_2nd.h>

#include <iostream>

int main()
{
    auto result = ctrlpp::online_planner_2nd<double>::create({.v_max = 1.0, .a_max = 5.0});
    if (!result.has_value())
        return 1;
    auto& planner = *result;
    if (!planner.update(10.0).has_value())
        return 1;

    constexpr double dt = 0.001;
    for (double t = 0.0; !planner.is_settled(); t += dt) {
        auto pt = planner.sample(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [online-planner-3rd](online-planner-3rd.md)<br/> 3rd-order variant with jerk limiting
- [trapezoidal-trajectory](trapezoidal-trajectory.md)<br/> Pre-computed trapezoidal profile
- [Trajectory Generation Theory](../../background/trajectory-generation.md)<br/> Online trajectory generation algorithms
