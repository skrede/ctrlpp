# modified_trap_trajectory

Modified trapezoidal motion law: the acceleration of a trapezoidal velocity profile, with its steps replaced by sinusoidal (cycloidal) ramps, so the acceleration is continuous. The motion is rest to rest, and the duration is user-specified (not solved from kinematic limits).

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/modified_trap_trajectory.h` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type |

Output dimension is fixed at ND=1 (scalar trajectory).

## Config

```cpp
struct config {
    Scalar q0;  // start position
    Scalar q1;  // end position
    Scalar T;   // total duration (user-specified)
};
```

`config` has no defaults; every field is required.

## Constructor

```cpp
explicit modified_trap_trajectory(config const& cfg);
```

The constructor validates nothing and there is no factory. `T` must be finite and strictly positive: `evaluate` divides by it, and `T = 0` evaluates to NaN. The direction of travel is the sign of `q1 - q0`; the profile is built for `h = abs(q1 - q0)` and mirrored for a negative displacement.

## Profile Structure

The profile has 6 phases with fixed T/8 boundaries, symmetric about T/2:

| Phase | Interval | Description |
|-------|----------|-------------|
| A-B | [0, T/8) | Cycloidal acceleration ramp-up |
| B-C | [T/8, 3T/8) | Constant acceleration |
| C-D | [3T/8, T/2] | Cycloidal acceleration ramp-down |
| D-E | [T/2, 5T/8) | Symmetric to C-D |
| E-F | [5T/8, 7T/8) | Symmetric to B-C |
| F-G | [7T/8, T] | Symmetric to A-B |

`t = T/2` belongs to the first half. The second half is computed via symmetry, in the positive-displacement frame: `q(t) = h - q(T-t)`, with the acceleration negated. `evaluate` clamps `t` to `[0, T]`.

## Key Constants

- **Peak velocity:** `2*h/T`, reached at `t = T/2`
- **Peak acceleration:** `8*pi/(2+pi) * h/T^2`, about `4.888*h/T^2`, held over the constant-acceleration phases

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, 1> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Total duration `T` |
| `peak_velocity` | `Scalar peak_velocity() const` | `2*h/T`, a nonnegative speed for either direction of travel |

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include "ctrlpp/trajectory/modified_trap_trajectory.h"

#include <iostream>

int main()
{
    ctrlpp::modified_trap_trajectory<double> traj({
        .q0 = 0.0, .q1 = 5.0, .T = 2.0
    });
    for (double t = 0; t <= traj.duration(); t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [trapezoidal-trajectory](trapezoidal-trajectory.md)<br/> standard trapezoidal (discontinuous acceleration)
- [modified-sin-trajectory](modified-sin-trajectory.md)<br/> alternative smooth profile
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
