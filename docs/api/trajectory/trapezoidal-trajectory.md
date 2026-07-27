# trapezoidal_trajectory

Trapezoidal velocity profile (LSPB) with three phases: acceleration, cruise, deceleration. Solves phase durations from kinematic limits `v_max` and `a_max`. Handles the triangular degenerate case when displacement is too short for full cruise phase, and non-null initial/final velocities with feasibility adjustment.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/trapezoidal_trajectory.h` |

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

When the two boundary velocities are not feasible over the commanded displacement at the commanded acceleration, the acceleration is raised to the smallest value that makes them feasible together (B&M eq. (3.15)). That raise is a division by the commanded displacement, so it has no representable answer once the displacement is small enough, and none at all when it is zero.

## Realizability

Both ramps of a three-phase profile run toward one cruise velocity lying at or above each boundary velocity, so the profile sweeps at least the ground the transition between those two velocities already sweeps. A command below that is not realizable within this shape at any acceleration the scalar type can hold; the clearest instance is a zero commanded displacement with two boundary speeds that differ, which asks the axis to change speed while covering no ground. `try_create` reports exactly that case:

```cpp
auto profile = ctrlpp::trapezoidal_trajectory<double>::try_create(cfg);
if (!profile) {
    // profile.error() == ctrlpp::trajectory_error::unreachable_boundary_velocity
}
```

The non-fallible constructor stays available for callers that have already established their command is realizable. Given one that is not, it yields a stationary zero-duration profile rather than one whose duration, phase durations, and evaluation are all NaN. Retiming that stand-in is `trajectory_error::unreachable_duration`, since it has no traversal to slow down.

Two zero-displacement commands are realizable and are not rejected: equal boundary velocities, where there is no speed change to cover, and opposed boundary velocities of equal magnitude, where the ramp between them sweeps exactly zero ground.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `try_create` | `static expected<trapezoidal_trajectory, trajectory_error> try_create(config const&)` | Construct, reporting an unrealizable command |
| `evaluate` | `trajectory_point<Scalar, 1> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Total duration `T = T_a + T_v + T_d` |
| `is_triangular` | `bool is_triangular() const` | True if cruise phase duration is zero |
| `peak_velocity` | `Scalar peak_velocity() const` | Signed peak velocity in original frame |
| `phase_durations` | `std::array<Scalar, 3> phase_durations() const` | `{T_accel, T_cruise, T_decel}` |
| `rescale_to` | `expected<void, trajectory_error> rescale_to(Scalar T_new)` | Rebuild at a longer duration for multi-axis sync |
| `can_rescale_to` | `expected<void, trajectory_error> can_rescale_to(Scalar T_new) const` | Whether `rescale_to(T_new)` would succeed, without mutating |

## Time Rescaling

`rescale_to()` rebuilds the profile at a lower cruise velocity rather than patching the one it has. The commanded displacement, both boundary velocities, and the acceleration magnitude are held fixed and the cruise velocity that realizes the requested duration is solved in closed form, so the traversed displacement and the terminal velocity hold by construction. The stored duration stays the sum of the three realized phase durations and is never assigned the requested value; it lands within a few units in the last place of it.

The solve covers all three shapes the three-phase parametrization admits, and picks between them by monotonicity: the total duration falls as the cruise velocity grows, so exactly one shape can contain the root.

| Shape | Validity | Solve |
|-------|----------|-------|
| plateau | cruise velocity at or above both boundary velocities | quadratic, smaller root |
| ramp-through | cruise velocity strictly between the two boundary velocities | linear in the reciprocal of the cruise velocity, single root |
| valley | cruise velocity at or below both boundary velocities | quadratic, larger root |

The valley shape is emitted, not rejected: with both boundary velocities above the cruise velocity a long duration needs, the profile decelerates away from the initial velocity, holds a low cruise velocity, and accelerates back up to the final one. A root is accepted only inside its own shape's validity interval, with all three phase durations nonnegative and the cruise velocity within the velocity limit; a root failing any of those is a rejection rather than a clamped value.

Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| a duration below the current one | `trajectory_error::duration_shorter_than_current` |
| NaN, infinite, or non-positive `T_new` | `trajectory_error::non_positive_duration` |
| a duration the displacement, limits, and boundary velocities cannot realize together | `trajectory_error::unreachable_duration` |

A request equal to the current duration succeeds and changes nothing, which is the path the slowest axis of a synchronized set always takes.

Reachability is derived rather than assumed. The cruise duration is what runs out: the shape whose validity interval reaches down toward a vanishing cruise velocity fixes the supremum of the reachable durations, and whether that supremum is finite is the sign of that shape's own residual displacement term. With both boundary velocities positive, the durations grow without bound exactly when the displacement exceeds `(v0^2 + v1^2) / (2 a)`; below that the cruise duration reaches zero at a strictly positive cruise velocity and the reachable durations stop at `(v0 + v1 - 2 sqrt((v0^2 + v1^2) / 2 - a h)) / a`. A stationary profile therefore reaches its own duration and nothing longer, with no epsilon taking part in the decision.

`can_rescale_to()` runs the identical solve and discards the result, so the two cannot disagree. That is what lets [`synchronize()`](synchronize.md) check every axis before it commits any of them.

## Triangular Degenerate Case

When the displacement is too small for the profile to reach `v_max`, the cruise phase collapses to zero duration. The profile becomes triangular: accelerate to `v_peak < v_max`, then immediately decelerate. This is detected and handled automatically during construction.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include "ctrlpp/trajectory/trapezoidal_trajectory.h"

#include <iostream>

int main()
{
    ctrlpp::trapezoidal_trajectory<double> traj({
        .q0 = 0.0, .q1 = 10.0,
        .v_max = 5.0, .a_max = 2.0
    });
    for (double t = 0; t <= traj.duration(); t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [double-s-trajectory](double-s-trajectory.md)<br/> jerk-limited alternative (7-segment)
- [modified-trap-trajectory](modified-trap-trajectory.md)<br/> smooth acceleration variant
- [time-scaling](time-scaling.md)<br/> duration computation for elementary paths
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
