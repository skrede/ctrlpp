# Real-Time Replanning

Generating trajectories on-the-fly with online planners.

## When to Use What

| Planner | Profile Shape | Constraints | Best For |
|---------|--------------|-------------|----------|
| [2nd-order](../../api/trajectory/online-planner-2nd.md) | Trapezoidal-like | v_max + a_max | Simpler real-time motion with acceleration limits |
| [3rd-order](../../api/trajectory/online-planner-3rd.md) | Double-S-like | v_max + a_max + j_max | Smoother real-time motion with jerk limits |

Both planners are designed for control-loop integration: call `update()` each
cycle and get the current reference position, velocity, and acceleration.

## Quick Start: 2nd-Order Planner

The 2nd-order online planner generates trapezoidal-like profiles that can be
updated with new targets at any time:

```cpp
#include <ctrlpp/trajectory/online_planner_2nd.h>

#include <iostream>

// create validates the kinematic limits (finite and strictly positive,
// because they divide in the planner math) and reports rejections through
// ctrlpp::expected; unwrap after checking.
auto planner_result =
    ctrlpp::online_planner_2nd<double>::create({.v_max = 5.0, .a_max = 2.0});
if (!planner_result.has_value()) {
    std::cerr << "invalid planner limits\n";
    return 1;
}
auto& planner = *planner_result;

planner.update(10.0);  // set a new target at any time; recomputes the profile

double dt = 0.01;  // 100 Hz control loop
for (int i = 0; i < 500; ++i) {
    double t = i * dt;
    auto point = planner.sample(t);  // evaluate the profile at time t
    std::cout << point.position[0] << "," << point.velocity[0] << "\n";
}
```

See [example 08](../../../examples/trajectory/ctrlpp_trajectory_08_online_planner.cpp)
for a runnable version.

## Target Changes Mid-Motion

A new target set while the planner is still moving does not always produce the
profile it was asked for, and `update()` returns nothing to say so. There are
two outcomes:

1. **The commanded profile.** The current velocity is carried straight through
   into the new move, with no full-stop dip. This is what a same-direction
   retarget with room to stop gets.
2. **A brake-then-replan.** The planner decelerates to zero (respecting the
   kinematic limits), then plans a new profile from the stopping point. This is
   what a reversal, an overshoot, or (on the 3rd-order planner) a state outside
   the carry-velocity shape's domain gets.

Both respect every limit and both reach the target. What differs is the time
taken, and nothing about the sampled motion says which one you got. Read the
disposition:

```cpp
// Change target while moving; the planner handles it safely
planner.update(10.0);
// ... some time later, before reaching 10.0 ...
planner.update(-5.0);  // a reversal at speed: brakes to zero, then heads to -5.0

auto const& diagnostics = planner.diagnostics();
if (diagnostics.disposition
    == ctrlpp::online_planner_disposition::braked_and_replanned)
{
    // Not the commanded profile. diagnostics.substitution_reason names the
    // condition, diagnostics.brake_duration is the time the braking costs, and
    // diagnostics.replan_start_position is where the new move begins.
}
```

A supervisory layer synchronizing several axes needs this: an axis on a
brake-then-replan is on a longer profile than the one commanded, and a
supervisor that assumes otherwise desynchronizes with nothing anywhere
reporting an error.

## Integration with Controllers

Online planner output feeds directly into PID or MPC controllers as the
reference signal. Each cycle, sample the planner at the loop's current time and
use the result as the setpoint:

```cpp
auto const point = planner.sample(t);

// Use point.position[0] as the PID setpoint, or position, velocity and
// acceleration together as feedforward terms in an MPC cost function.
```

See [example 10](../../../examples/trajectory/ctrlpp_trajectory_10_mpc_tracking.cpp)
for a complete trajectory-tracking MPC example.

## Links

**API reference:**
[2nd-order planner](../../api/trajectory/online-planner-2nd.md) |
[3rd-order planner](../../api/trajectory/online-planner-3rd.md)

**Examples:**
[08 online planner](../../../examples/trajectory/ctrlpp_trajectory_08_online_planner.cpp) |
[10 MPC tracking](../../../examples/trajectory/ctrlpp_trajectory_10_mpc_tracking.cpp)
