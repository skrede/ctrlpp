# double_s_trajectory

Double-S (7-segment) velocity profile with jerk-limited motion. Computes a time-optimal S-curve that respects velocity, acceleration, and jerk constraints simultaneously. The 7 segments are: jerk(+), const-accel, jerk(-), cruise, jerk(-), const-decel, jerk(+).

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/double_s_trajectory.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

Output dimension is fixed at ND=1 (scalar trajectory).

## Config

```cpp
struct config {
    Scalar q0, q1;              // start/end positions
    Scalar v_max, a_max, j_max; // kinematic limits
    Scalar v0{}, v1{};          // initial/final velocities (default 0)
};
```

## Constructor

```cpp
explicit double_s_trajectory(config const& cfg);
```

Construction follows the B&M flowchart (Fig 3.18) to solve all phase durations. Handles negative displacement via sigma transformation. Zero displacement produces a stationary profile.

## 7-Segment Structure

```
velocity
  ^
  |     /------\
  |    / |    | \
  |   /  |    |  \
  |  /   |    |   \
  | /    |    |    \
  +--+---+----+---+---> time
   1  2   3  4  5  6  7

1: jerk(+)  2: const-accel  3: jerk(-)
4: cruise
5: jerk(-)  6: const-decel  7: jerk(+)
```

## Degenerate Cases

The profile degenerates when kinematic limits cannot all be reached:

- **v_max not reached:** Cruise phase (segment 4) collapses to zero
- **a_max not reached:** Constant-accel/decel phases (segments 2,6) collapse to zero
- **Doubly degenerate:** Both v_max and a_max unreachable; purely jerk-limited profile with `T_j = cbrt(h / (2*j_max))`

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, 1> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Total duration |
| `is_degenerate` | `bool is_degenerate() const` | True if v_max or a_max not reached |
| `peak_velocity` | `Scalar peak_velocity() const` | Actual peak velocity achieved |
| `phase_durations` | `std::array<Scalar, 7> phase_durations() const` | Per-segment durations |
| `rescale_to` | `void rescale_to(Scalar T_new)` | Extend cruise phase for multi-axis sync |

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <iostream>

int main()
{
    ctrlpp::double_s_trajectory<double> traj({
        .q0 = 0.0, .q1 = 10.0,
        .v_max = 5.0, .a_max = 10.0, .j_max = 50.0
    });
    for (double t = 0; t <= traj.duration(); t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [trapezoidal-trajectory](trapezoidal-trajectory.md) -- simpler 3-segment alternative
- [modified-trap-trajectory](modified-trap-trajectory.md) -- smooth acceleration variant
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
