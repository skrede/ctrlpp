# Point-to-Point Motion

Choosing the right trajectory profile for single-segment moves.

## When to Use What

| Profile | Continuity | Constraints | Best For |
|---------|-----------|-------------|----------|
| [Cubic](../../api/trajectory/cubic-trajectory.md) | C1 | pos + vel BCs | Simple moves without acceleration control |
| [Quintic](../../api/trajectory/quintic-trajectory.md) | C2 | pos + vel + acc BCs | Smooth motion with zero-acceleration endpoints |
| [Septic](../../api/trajectory/septic-trajectory.md) | C3 | pos + vel + acc + jerk BCs | Ultra-smooth when jerk matters |
| [Trapezoidal](../../api/trajectory/trapezoidal-trajectory.md) | C1 | v_max + a_max | Industrial moves with clear speed limits |
| [Double-S](../../api/trajectory/double-s-trajectory.md) | C2 | v_max + a_max + j_max | Vibration-sensitive applications |
| [Modified trapezoidal](../../api/trajectory/modified-trap-trajectory.md) | C1 | duration-based | Smoother than standard trapezoidal |
| [Modified sinusoidal](../../api/trajectory/modified-sin-trajectory.md) | C1 | duration-based | Minimum residual vibration |

## Quick Start: Polynomial

A quintic trajectory guarantees zero velocity and acceleration at both
endpoints. Specify boundary conditions and evaluate at any time:

```cpp
#include <ctrlpp/trajectory/quintic_trajectory.h>

#include <Eigen/Core>
#include <iostream>

Eigen::Vector3d p0{0.0, 0.0, 0.0};
Eigen::Vector3d p1{1.0, 2.0, 0.5};

auto traj = ctrlpp::make_quintic_trajectory(
    0.0, p0, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(),
    2.0, p1, Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero());

auto [pos, vel, acc] = traj.evaluate(1.0);
std::cout << "position at t=1: " << pos.transpose() << "\n";
```

See [example 02](../../../examples/trajectory/ctrlpp_trajectory_02_quintic.cpp)
for a runnable version.

## Quick Start: Velocity Profile

Trapezoidal profiles work with kinematic limits (maximum velocity and
acceleration) rather than polynomial boundary conditions:

```cpp
#include <ctrlpp/trajectory/trapezoidal_trajectory.h>

#include <iostream>

ctrlpp::trapezoidal_config<double> cfg{};
cfg.q0    = 0.0;
cfg.q1    = 10.0;
cfg.v_max = 5.0;
cfg.a_max = 2.0;

auto traj = ctrlpp::trapezoidal_trajectory(cfg);
auto [pos, vel, acc] = traj.evaluate(1.5);
std::cout << "position at t=1.5: " << pos << "\n";
```

See [example 03](../../../examples/trajectory/ctrlpp_trajectory_03_trapezoidal.cpp)
for a runnable version.

## Choosing by Smoothness

- **Position control only** -- cubic is sufficient (C1 continuity).
- **Smooth velocity needed** -- quintic guarantees continuous acceleration (C2).
- **Acceleration matters** (e.g., force control) -- double-S or septic provide
  jerk-limited or C3 profiles.
- **Duration-based profiling** -- modified trapezoidal or modified sinusoidal
  when you know the total move time and want reduced residual vibration.

## Kinematic Time Scaling

If you need to compute the minimum move duration that respects velocity,
acceleration, and jerk limits before creating a trajectory, use the
[time scaling](../../api/trajectory/time-scaling.md) utilities. These compute
feasible durations from kinematic constraints, which you can then pass to
any profile constructor.

## Links

**API reference:**
[cubic](../../api/trajectory/cubic-trajectory.md) |
[quintic](../../api/trajectory/quintic-trajectory.md) |
[septic](../../api/trajectory/septic-trajectory.md) |
[trapezoidal](../../api/trajectory/trapezoidal-trajectory.md) |
[double-S](../../api/trajectory/double-s-trajectory.md) |
[modified trapezoidal](../../api/trajectory/modified-trap-trajectory.md) |
[modified sinusoidal](../../api/trajectory/modified-sin-trajectory.md) |
[time scaling](../../api/trajectory/time-scaling.md)

**Examples:**
[01 cubic](../../../examples/trajectory/ctrlpp_trajectory_01_cubic.cpp) |
[02 quintic](../../../examples/trajectory/ctrlpp_trajectory_02_quintic.cpp) |
[03 trapezoidal](../../../examples/trajectory/ctrlpp_trajectory_03_trapezoidal.cpp) |
[04 double-S](../../../examples/trajectory/ctrlpp_trajectory_04_double_s.cpp) |
[05 modified profiles](../../../examples/trajectory/ctrlpp_trajectory_05_modified_profiles.cpp)
