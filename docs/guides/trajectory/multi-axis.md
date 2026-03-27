# Multi-Axis Coordination

Synchronizing multiple trajectory axes to move together.

## The Problem

When multiple axes need to reach their targets simultaneously (e.g., a 3-axis
CNC machine cutting a straight diagonal line), each axis typically has a
different distance to travel. If each axis plans its own trajectory
independently, they finish at different times -- causing curved tool paths
instead of straight lines in Cartesian space.

Synchronization solves this by stretching the faster axes so that all axes
finish at the same time, while still respecting each axis's kinematic limits.

## Quick Start: synchronize()

Create individual per-axis trajectories, then synchronize them so they all
share the same duration:

```cpp
#include <ctrlpp/trajectory/trapezoidal_trajectory.h>
#include <ctrlpp/trajectory/synchronize.h>

#include <iostream>

// Three axes with different distances
ctrlpp::trapezoidal_config<double> cfg_x{.q0=0, .q1=10, .v_max=5, .a_max=2};
ctrlpp::trapezoidal_config<double> cfg_y{.q0=0, .q1=3,  .v_max=5, .a_max=2};
ctrlpp::trapezoidal_config<double> cfg_z{.q0=0, .q1=1,  .v_max=5, .a_max=2};

auto ax_x = ctrlpp::trapezoidal_trajectory(cfg_x);
auto ax_y = ctrlpp::trapezoidal_trajectory(cfg_y);
auto ax_z = ctrlpp::trapezoidal_trajectory(cfg_z);

// Synchronize: all axes now finish at the same time
ctrlpp::synchronize(ax_x, ax_y, ax_z);

// Evaluate at any time -- all axes are coordinated
double t = 2.0;
auto [px, vx, ax] = ax_x.evaluate(t);
auto [py, vy, ay] = ax_y.evaluate(t);
auto [pz, vz, az] = ax_z.evaluate(t);
std::cout << "pos: " << px << ", " << py << ", " << pz << "\n";
```

See [example 09](../../../examples/trajectory/ctrlpp_trajectory_09_multi_axis_sync.cpp)
for a runnable version.

## How It Works

The `synchronize()` function:

1. Finds the maximum duration among all provided trajectories.
2. Rescales each shorter trajectory to match that maximum duration via
   `rescale_to()`.
3. Each axis's kinematic limits are still respected -- the rescaled profile
   simply uses lower velocities and accelerations for shorter moves.

The result is that all axes start and finish together, with coordinated
velocity profiles throughout the motion.

## Supported Profile Types

Any trajectory type satisfying the `syncable_profile` concept can be
synchronized. A type is syncable if it provides:

- `duration()` -- returns the trajectory duration
- `rescale_to(new_duration)` -- stretches the profile to a new duration
- `scalar_type` -- the scalar type alias

Currently, [trapezoidal](../../api/trajectory/trapezoidal-trajectory.md) and
[double-S](../../api/trajectory/double-s-trajectory.md) profiles support
synchronization out of the box.

## Links

**API reference:**
[synchronize](../../api/trajectory/synchronize.md) |
[trapezoidal](../../api/trajectory/trapezoidal-trajectory.md) |
[double-S](../../api/trajectory/double-s-trajectory.md)

**Examples:**
[09 multi-axis sync](../../../examples/trajectory/ctrlpp_trajectory_09_multi_axis_sync.cpp)
