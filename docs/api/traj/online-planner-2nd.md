# online_planner_2nd

2nd-order online trajectory planner that generates trapezoidal-like velocity profiles in real time. On each `update(target)`, the planner computes a time-optimal profile from the current state to the target position, respecting `v_max` and `a_max` constraints. Suitable for real-time control loops where targets change dynamically.

Unlike pre-computed trajectory segments, online planners are stateful filters with no fixed duration and do not satisfy `trajectory_segment`.

## Header

```cpp
#include "ctrlpp/traj/online_planner_2nd.h"
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

## Constructor

```cpp
explicit online_planner_2nd(config const& cfg);
```

Constructs a planner with kinematic limits. Initial state is at rest at q = 0.

## Methods

### update

```cpp
void update(Scalar target);
```

Set a new target position and replan from the current state. Computes a time-optimal trapezoidal profile from (q, v) to (target, 0) respecting `v_max` and `a_max`. Handles overshoot recovery when current velocity requires braking before replanning.

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

## Profile Behaviour

The planner produces trapezoidal velocity profiles with three phases:

1. **Acceleration** -- ramp velocity toward `v_max` at rate `a_max`
2. **Cruise** -- hold at `v_max` (may be zero duration for short moves)
3. **Deceleration** -- ramp velocity to zero at rate `a_max`

For short displacements where `v_max` cannot be reached, the profile degenerates to a triangular velocity shape. Mid-motion target changes trigger replanning from the current state, with automatic brake-and-replan for overshoot scenarios.

## Usage Example

```cpp
#include "ctrlpp/traj/online_planner_2nd.h"

ctrlpp::online_planner_2nd<double> planner({.v_max = 1.0, .a_max = 5.0});

planner.update(10.0);  // move to position 10

double dt = 0.001;
for (double t = 0.0; !planner.is_settled(); t += dt) {
    auto pt = planner.sample(t);
    // Send pt.position to servo loop
}
```

## See Also

- [online-planner-3rd](online-planner-3rd.md) -- 3rd-order variant with jerk limiting
- [trapezoidal-trajectory](trapezoidal-trajectory.md) -- Pre-computed trapezoidal profile
- [Trajectory Generation Theory](../../background/trajectory-generation.md) -- Online trajectory generation algorithms
