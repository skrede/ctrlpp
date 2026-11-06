# Point-to-Point Motion

Choosing the right trajectory profile for single-segment moves.

## When to Use What

| Profile | Continuity | Constraints | Best For |
|---------|-----------|-------------|----------|
| [Cubic](../../api/trajectory/cubic-trajectory.md) | C1 | pos + vel BCs | Simple moves without acceleration control |
| [Quintic](../../api/trajectory/quintic-trajectory.md) | C2 | pos + vel + acc BCs | Smooth motion with set boundary accelerations |
| [Septic](../../api/trajectory/septic-trajectory.md) | C3 | pos + vel + acc + jerk BCs | Ultra-smooth when jerk matters |
| [Trapezoidal](../../api/trajectory/trapezoidal-trajectory.md) | C1 | v_max + a_max | Industrial moves with clear speed limits |
| [Double-S](../../api/trajectory/double-s-trajectory.md) | C2 | v_max + a_max + j_max | Vibration-sensitive applications |
| [Modified trapezoidal](../../api/trajectory/modified-trap-trajectory.md) | C2 | q0, q1, duration | Smoother than standard trapezoidal |
| [Modified sinusoidal](../../api/trajectory/modified-sin-trajectory.md) | C2 | q0, q1, duration | Rest-to-rest moves over a set duration |

## Quick Start: Polynomial

A quintic trajectory takes position, velocity and acceleration at both
endpoints and a duration. With zero boundary velocities and accelerations it
starts and ends at rest. Build it and evaluate at any time:

```cpp
#include <ctrlpp/trajectory/quintic_trajectory.h>

#include <Eigen/Core>

#include <iostream>

int main()
{
    Eigen::Vector3d const p0{0.0, 0.0, 0.0};
    Eigen::Vector3d const p1{1.0, 2.0, 0.5};
    Eigen::Vector3d const zero = Eigen::Vector3d::Zero();

    auto const created = ctrlpp::make_quintic_trajectory(p0, p1, zero, zero, zero, zero, 2.0);
    if (!created) {
        std::cerr << "invalid trajectory configuration\n";
        return 1;
    }
    auto [pos, vel, acc] = created->evaluate(1.0);
    std::cout << "position at t=1: " << pos.transpose() << "\n";
}
```

See [example 02](../../../examples/trajectory/ctrlpp_trajectory_02_quintic.cpp)
for a runnable version.

## Quick Start: Velocity Profile

Trapezoidal profiles work with kinematic limits (maximum velocity and
acceleration) rather than polynomial boundary conditions:

```cpp
#include <ctrlpp/trajectory/trapezoidal_trajectory.h>

#include <iostream>

int main()
{
    // create() is the only construction path. It returns the profile or the
    // reason the command has none, so an unrealizable command can never be
    // mistaken for a profile that quietly does nothing.
    auto const created = ctrlpp::trapezoidal_trajectory<double>::create({
        .q0 = 0.0, .q1 = 10.0, .v_max = 5.0, .a_max = 2.0});
    if (!created) {
        std::cerr << "The commanded move has no trapezoidal profile\n";
        return 1;
    }
    auto [pos, vel, acc] = created->evaluate(1.5);
    std::cout << "position at t=1.5: " << pos << "\n";
}
```

See [example 03](../../../examples/trajectory/ctrlpp_trajectory_03_trapezoidal.cpp)
for a runnable version.

## Choosing by Smoothness

- **Position control only**<br/>cubic is sufficient (C1 continuity).
- **Smooth velocity needed**<br/>quintic guarantees continuous acceleration (C2).
- **Acceleration matters** (e.g., force control)<br/> double-S or septic provide
  jerk-limited or C3 profiles.
- **Duration-based profiling**<br/>modified trapezoidal or modified sinusoidal
  when you know the total move time and want reduced residual vibration.

## Kinematic Time Scaling

If you need the minimum duration at which a normalized path respects
velocity, acceleration, and jerk limits, use the
[time scaling](../../api/trajectory/time-scaling.md) utilities. They take the
path's peak normalized derivatives and the limits, and return a duration to
pass to [`make_trajectory`](../../api/trajectory/trajectory.md) with the same
path. The trapezoidal and double-S profiles take no duration: they solve their
own timing from the limits.

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
