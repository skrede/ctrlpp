# Multi-Axis Coordination

Synchronizing multiple trajectory axes to move together.

## The Problem

When multiple axes need to reach their targets simultaneously (e.g., a 3-axis
CNC machine cutting a straight diagonal line), each axis typically has a
different distance to travel. If each axis plans its own trajectory
independently, they finish at different times &mdash; causing curved tool paths
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

// Synchronize: all axes now finish at the same time. The call is fallible --
// an axis that cannot reach the slowest duration is reported, and nothing is
// retimed unless every axis can be.
auto const synced = ctrlpp::synchronize(ax_x, ax_y, ax_z);
if (!synced) {
    std::cerr << "Synchronization rejected\n";
    return 1;
}

// Evaluate at any time &mdash; all axes are coordinated
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
2. Asks every axis, through `can_rescale_to()`, whether it can reach that
   duration.
3. Only once every axis has answered yes does it retime any of them, through
   `rescale_to()`. A rejection returns that axis's `trajectory_error` with
   nothing mutated, so a partially synchronized set is never handed back.

Each axis is rebuilt at the longer duration rather than patched, so its
commanded displacement and its commanded boundary velocities still hold
afterwards and its kinematic limits are still respected: the rebuilt profile
simply runs at lower velocities and accelerations.

Retiming is only ever a slowdown. An axis whose displacement, limits, and
boundary velocities cannot reach the target duration is reported rather than
retimed to something else -- with nonzero boundary velocities the reachable
durations are bounded, because the velocity limit cannot fall below the speed
the axis was commanded to enter or leave with.

The result is that all axes start and finish together, with coordinated
velocity profiles throughout the motion.

## Runtime-Sized Sets

For a runtime-sized set of identical profiles, pass a contiguous view. The view
carries no ownership, so the same call serves a vector, a plain array, or a
statically allocated buffer, and the call allocates nothing:

```cpp
#include <span>
#include <vector>

std::vector<ctrlpp::trapezoidal_trajectory<double>> axes = /* ... */;
auto const synced = ctrlpp::synchronize(std::span{axes});
```

## Supported Profile Types

Any trajectory type satisfying the `syncable_profile` concept can be
synchronized. A type is syncable if it provides:

- `duration()`, which returns the trajectory duration
- `rescale_to(new_duration)`, which retimes the profile and returns
  `expected<void, trajectory_error>`
- `can_rescale_to(new_duration) const`, which answers the same question without
  mutating and has the same return
- `scalar_type`, which is the scalar type alias

Currently, [trapezoidal](../../api/trajectory/trapezoidal-trajectory.md) and [double-S](../../api/trajectory/double-s-trajectory.md) profiles support synchronization out of the box.

## Links

**API reference:**
[synchronize](../../api/trajectory/synchronize.md) |
[trapezoidal](../../api/trajectory/trapezoidal-trajectory.md) |
[double-S](../../api/trajectory/double-s-trajectory.md)

**Examples:**
[09 multi-axis sync](../../../examples/trajectory/ctrlpp_trajectory_09_multi_axis_sync.cpp)
