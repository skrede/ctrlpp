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

## Construction

```cpp
static expected<double_s_trajectory, trajectory_error> create(config const& cfg);
```

`create` is the **only** construction path. There is no public non-fallible constructor: a profile object either satisfies its contract or it was never built, so an unrealizable command is a value the caller has to inspect rather than a profile that quietly stands in for one.

```cpp
auto profile = ctrlpp::double_s_trajectory<double>::create(cfg);
if (!profile) {
    // profile.error() names which part of the contract the command failed
}
```

Construction follows the B&M flowchart (Fig 3.18) to solve all phase durations. Handles negative displacement via sigma transformation. A zero displacement produces the standstill, and only when both boundary velocities are zero as well; the zero-displacement branch is keyed on exact equality, so a displacement small enough to underflow the rest of the algebra is a rejection rather than a command rounded down to a resting axis.

Nonzero initial and final velocities `v0`, `v1` are supported through the general B&M Sec 3.4.1 formulation: the acceleration phase ramps from `v0` and the deceleration phase ramps to `v1`, so `evaluate(0)` reports `v0` and `evaluate(duration())` reports `v1`, both with zero acceleration. When `v0 = v1 = 0` the solver reduces exactly to the symmetric Sec 3.4.3 special case.

In the no-cruise sub-case the peak velocity is solved directly against the commanded displacement, and each ramp takes whichever of its two closed forms applies -- constant-acceleration segment present, or triangular in acceleration when the velocity change is too small to build up to `a_max`. All three limits stay respected in both shapes. Where the two ramps carry different square roots of the peak the solve falls back to a bracketed search that terminates by exhaustion of the floating-point bracket, not on an iteration count.

## Rejections

Checked in order:

| Condition | Error |
|-----------|-------|
| NaN or infinite `q0`, `q1`, `v0`, or `v1` | `trajectory_error::non_finite_input` |
| NaN, infinite, or non-positive `v_max` | `trajectory_error::non_positive_velocity_limit` |
| NaN, infinite, or non-positive `a_max` | `trajectory_error::non_positive_acceleration_limit` |
| NaN, infinite, or non-positive `j_max` | `trajectory_error::non_positive_jerk_limit` |
| `abs(v0) > v_max` or `abs(v1) > v_max` | `trajectory_error::boundary_velocity_exceeds_limit` |
| a displacement below the distance the transition between the boundary velocities already sweeps | `trajectory_error::unreachable_boundary_velocity` |
| a negative segment duration or total duration | `trajectory_error::unreachable_boundary_velocity` |
| a duration outside the representable range, or a total that underflowed to zero on a nonzero displacement | `trajectory_error::unrepresentable_duration` |

### The velocity limit is a precondition

`abs(v0) <= v_max` and `abs(v1) <= v_max` are preconditions of the type, not invariants construction restores. Raising the limit to fit a boundary velocity would return a profile that violates a bound the caller stated, which is the same silent-success defect a rejection exists to prevent. Reaching the limit exactly is inside the domain: the ramp attached to that boundary velocity simply vanishes.

### Realizability

A seven-segment profile cannot sweep less ground than the fastest admissible transition from the larger of the two boundary velocities to the smaller one, so a command shorter than that distance is not realizable within this shape, because reaching it would require overshooting the target and returning. That minimum is the swept distance of the cruise-free profile whose peak rises by nothing above the larger boundary velocity, computed from the same ramp expressions the profile is built out of rather than from a separate approximation of them.

A zero commanded displacement is realizable by exactly one member of the family, the standstill, and only when both boundary velocities are zero. Anything else asks the axis to change speed while covering no ground.

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
| `create` | `static expected<double_s_trajectory, trajectory_error> create(config const&)` | The only construction path; reports an unrealizable command |
| `evaluate` | `trajectory_point<Scalar, 1> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Total duration |
| `is_degenerate` | `bool is_degenerate() const` | True if v_max or a_max not reached |
| `peak_velocity` | `Scalar peak_velocity() const` | Actual peak velocity achieved |
| `phase_durations` | `std::array<Scalar, 7> phase_durations() const` | Per-segment durations |
| `rescale_to` | `expected<void, trajectory_error> rescale_to(Scalar T_new)` | Rebuild under scaled limits for multi-axis sync |
| `can_rescale_to` | `expected<void, trajectory_error> can_rescale_to(Scalar T_new) const` | Whether `rescale_to(T_new)` would succeed, without mutating |

## Time Rescaling

`rescale_to()` rebuilds the profile under scaled kinematic limits rather than patching the one it has. A time scaling that slows a profile by a factor divides its velocity by that factor, its acceleration by its square, and its jerk by its cube, so the three limits are multiplied by the first, second, and third power of one scale factor and the profile is constructed again from the same command.

The rebuild goes back through `create`, so the scaled command is held to the same contract the original was and a scale that produces no profile is reported rather than applied.

The two boundary velocities are left **unscaled**. They are what the caller commanded the axis to enter and leave with, not limits, and scaling them would land a synchronized axis at the wrong terminal velocity. Displacement, terminal velocity, and continuity then hold by construction. The stored duration stays whatever the rebuilt profile realizes and is never assigned the requested value; it lands within a few units in the last place of it.

With both boundary velocities at rest the duration is exactly proportional to the reciprocal of the scale, so the scale is the ratio of the current duration to the requested one and a single rebuild settles it. With a nonzero boundary velocity that proportionality fails, and the scale is found by halving a bracket that runs from the scale at which the commanded boundary velocities themselves reach the scaled velocity limit up to the profile's own scale. Termination is bracket exhaustion, when the midpoint falls on an endpoint: no tolerance and no iteration cap. The worst case is derived rather than chosen. Halving an interval whose endpoints share a binary exponent reaches the spacing of the representable values after one step more than the significand width, which is 25 evaluations at single precision and 54 at double precision; a bracket spanning several exponents costs one further step per exponent spanned. A derivative step is not used: the duration carries real kinks where the segment shape flips, so a Newton step can leave the bracket, and bounding a safeguarded variant would require an iteration cap.

Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| a duration below the current one | `trajectory_error::duration_shorter_than_current` |
| NaN, infinite, or non-positive `T_new` | `trajectory_error::non_positive_duration` |
| a duration no admissible scale realizes | `trajectory_error::unreachable_duration` |

A request equal to the current duration succeeds and changes nothing, which is the path the slowest axis of a synchronized set always takes. Because the boundary velocities stay fixed while the velocity limit falls, the reachable durations are bounded: the limit cannot drop below the speed the caller commanded the axis to enter or leave with.

`can_rescale_to()` runs the identical solve and discards the result. No closed-form reachability predicate exists for this family, so replaying the same deterministic computation is what makes the check and the commit agree, and it is what lets [`synchronize()`](synchronize.md) check every axis before it commits any of them.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include "ctrlpp/trajectory/double_s_trajectory.h"

#include <iostream>

int main()
{
    auto const created = ctrlpp::double_s_trajectory<double>::create({
        .q0 = 0.0, .q1 = 10.0,
        .v_max = 5.0, .a_max = 10.0, .j_max = 50.0
    });
    if (!created) {
        std::cerr << "The commanded move has no double-S profile\n";
        return 1;
    }
    auto const& traj = created.value();
    for (double t = 0; t <= traj.duration(); t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [trapezoidal-trajectory](trapezoidal-trajectory.md)<br/> simpler 3-segment alternative
- [modified-trap-trajectory](modified-trap-trajectory.md)<br/> smooth acceleration variant
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
