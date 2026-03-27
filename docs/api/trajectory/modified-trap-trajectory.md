# modified_trap_trajectory

Modified trapezoidal velocity profile with cycloidal acceleration phases. Replaces the constant-acceleration ramps of a standard trapezoidal profile with sinusoidal transitions, producing continuous acceleration and reduced vibration excitation. Duration is user-specified (not solved from kinematic limits).

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/modified_trap_trajectory.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

Output dimension is fixed at ND=1 (scalar trajectory).

## Config

```cpp
struct config {
    Scalar q0;  // start position
    Scalar q1;  // end position
    Scalar T;   // total duration (user-specified)
};
```

## Constructor

```cpp
explicit modified_trap_trajectory(config const& cfg);
```

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

The second half is computed via symmetry: `q(t) = h - q(T-t)`.

## Key Constants

- **Peak velocity:** `2*h/T`
- **Peak acceleration:** `2*h*(2+pi) / (pi*T^2)`

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, 1> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Total duration `T` |
| `peak_velocity` | `Scalar peak_velocity() const` | `2*h/T` |

## Usage Example

```cpp
#include "ctrlpp/traj/modified_trap_trajectory.h"

ctrlpp::modified_trap_trajectory<double> traj({
    .q0 = 0.0, .q1 = 5.0, .T = 2.0
});

auto pt = traj.evaluate(1.0);  // midpoint
// Continuous acceleration profile -- no discontinuities
```

## See Also

- [trapezoidal-trajectory](trapezoidal-trajectory.md) -- standard trapezoidal (discontinuous acceleration)
- [modified-sin-trajectory](modified-sin-trajectory.md) -- alternative smooth profile
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
