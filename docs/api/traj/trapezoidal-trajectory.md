# trapezoidal_trajectory

Trapezoidal velocity profile (LSPB) with three phases: acceleration, cruise, deceleration. Solves phase durations from kinematic limits `v_max` and `a_max`. Handles the triangular degenerate case when displacement is too short for full cruise phase, and non-null initial/final velocities with feasibility adjustment.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/trapezoidal_trajectory.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

Output dimension is fixed at ND=1 (scalar trajectory).

## Config

```cpp
struct config {
    Scalar q0, q1;         // start/end positions
    Scalar v_max, a_max;   // kinematic limits
    Scalar v0{}, v1{};     // initial/final velocities (default 0)
};
```

## Constructor

```cpp
explicit trapezoidal_trajectory(config const& cfg);
```

Construction solves phase durations from the kinematic constraints. Negative displacement is handled via sigma transformation. When `v_max` cannot be reached, the profile degenerates to a triangular shape with `v_peak = sqrt((2*a*h + v0^2 + v1^2) / 2)`.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, 1> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Total duration `T = T_a + T_v + T_d` |
| `is_triangular` | `bool is_triangular() const` | True if cruise phase duration is zero |
| `peak_velocity` | `Scalar peak_velocity() const` | Signed peak velocity in original frame |
| `phase_durations` | `std::array<Scalar, 3> phase_durations() const` | `{T_accel, T_cruise, T_decel}` |
| `rescale_to` | `void rescale_to(Scalar T_new)` | Extend to longer duration for multi-axis sync |

## Triangular Degenerate Case

When the displacement is too small for the profile to reach `v_max`, the cruise phase collapses to zero duration. The profile becomes triangular: accelerate to `v_peak < v_max`, then immediately decelerate. This is detected and handled automatically during construction.

## Usage Example

```cpp
#include "ctrlpp/traj/trapezoidal_trajectory.h"

ctrlpp::trapezoidal_trajectory<double> traj({
    .q0 = 0.0, .q1 = 10.0,
    .v_max = 5.0, .a_max = 2.0
});

auto pt = traj.evaluate(1.0);
// pt.position, pt.velocity, pt.acceleration

if (traj.is_triangular()) {
    // cruise phase was too short -- triangular profile used
}

auto [T_a, T_v, T_d] = traj.phase_durations();
```

## See Also

- [double-s-trajectory](double-s-trajectory.md) -- jerk-limited alternative (7-segment)
- [modified-trap-trajectory](modified-trap-trajectory.md) -- smooth acceleration variant
- [time-scaling](time-scaling.md) -- duration computation for elementary paths
- [Trajectory Generation Theory](../../reference/trajectory-generation.md)
