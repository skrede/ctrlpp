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

ctrlpp::online_planner_2nd_config<double> cfg{};
cfg.v_max = 5.0;
cfg.a_max = 2.0;

ctrlpp::online_planner_2nd<double> planner(cfg);
planner.set_target(10.0);

double dt = 0.01;  // 100 Hz control loop
for (int i = 0; i < 500; ++i) {
    planner.update(dt);
    auto [pos, vel, acc] = planner.sample();
    std::cout << pos << "," << vel << "\n";
}
```

See [example 08](../../../examples/trajectory/ctrlpp_trajectory_08_online_planner.cpp)
for a runnable version.

## Target Changes Mid-Motion

Online planners handle target changes safely using a brake-to-zero-then-replan
strategy. When a new target is set while the planner is still moving:

1. The planner decelerates to zero velocity (respecting kinematic limits).
2. Once stopped, it plans a new profile toward the updated target.

This approach is robust and predictable -- the system never attempts to
reverse direction at speed. For applications that need smoother transitions,
the 3rd-order planner additionally limits jerk during the deceleration phase.

```cpp
// Change target while moving -- planner handles it safely
planner.set_target(10.0);
// ... some time later, before reaching 10.0 ...
planner.set_target(-5.0);  // brakes to zero, then heads to -5.0
```

## Integration with Controllers

Online planner output feeds directly into PID or MPC controllers as the
reference signal. Each cycle, sample the planner state and use it as the
setpoint:

```cpp
planner.update(dt);
auto [ref_pos, ref_vel, ref_acc] = planner.sample();

// Use ref_pos as PID setpoint, or ref_pos + ref_vel + ref_acc as
// feedforward terms in an MPC cost function.
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
