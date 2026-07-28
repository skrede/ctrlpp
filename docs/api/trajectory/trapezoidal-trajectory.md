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

## Construction

```cpp
static expected<trapezoidal_trajectory, trajectory_error> create(config const& cfg);
```

`create` is the **only** construction path. There is no public non-fallible constructor: a profile object either satisfies its contract or it was never built, so an unrealizable command is a value the caller has to inspect rather than a profile that quietly stands in for one.

```cpp
auto profile = ctrlpp::trapezoidal_trajectory<double>::create(cfg);
if (!profile) {
    // profile.error() names which part of the contract the command failed
}
```

Construction solves phase durations from the kinematic constraints. Negative displacement is handled via sigma transformation. When `v_max` cannot be reached, the profile degenerates to a triangular shape with `v_peak = sqrt((2*a*h + v0^2 + v1^2) / 2)`; that is one of the shapes this family covers and is reported by `is_triangular()`, not a rejection.

When the two boundary velocities are not feasible over the commanded displacement at the commanded acceleration, the acceleration is raised to the smallest value that makes them feasible together (B&M eq. (3.15)). That raise is a division by the commanded displacement, so it has no representable answer once the displacement is small enough, and none at all when it is zero.

## The acceleration raise is a disposition, not a rejection

A raised acceleration is a **success whose realized limit differs from the commanded one**, and it is reported through `disposition()` rather than through the failure channel.

```cpp
template <typename Scalar>
struct trapezoidal_disposition
{
    Scalar commanded_acceleration{};  // the a_max passed to create
    Scalar realized_acceleration{};   // the magnitude every ramp actually uses
};
```

```cpp
auto const& d = traj.disposition();
if(d.realized_acceleration > d.commanded_acceleration)
{
    // the profile is valid and respects d.realized_acceleration in every
    // phase; it does NOT respect the limit that was commanded
}
```

The profile that comes back is correct and limit-respecting **under the realized limit**, and it is the profile the command asks for at the only acceleration that can deliver it. So it is not a failure, and it does not become one: forcing a caller to handle a non-failure through the failure path teaches the caller that the failure path is usually noise.

It is not a flag either. A supervisory layer that learns only that something was substituted cannot decide anything with it; the commanded and realized values are what let a caller whose acceleration limit is physical rather than advisory compare the two and act. `realized_acceleration == commanded_acceleration` exactly when nothing was raised.

The disposition is fixed when the profile is built and is never recomputed: the raise is decided once, and `rescale_to` holds the acceleration magnitude fixed by construction. No enumerator was added to `trajectory_error` for it.

## Rejections

Checked in order:

| Condition | Error |
|-----------|-------|
| NaN or infinite `q0`, `q1`, `v0`, or `v1` | `trajectory_error::non_finite_input` |
| NaN, infinite, or non-positive `v_max` | `trajectory_error::non_positive_velocity_limit` |
| NaN, infinite, or non-positive `a_max` | `trajectory_error::non_positive_acceleration_limit` |
| `abs(v0) > v_max` or `abs(v1) > v_max` | `trajectory_error::boundary_velocity_exceeds_limit` |
| a displacement the two boundary velocities cannot be reconciled with at a representable acceleration | `trajectory_error::unreachable_boundary_velocity` |
| a negative phase duration or total duration | `trajectory_error::unreachable_boundary_velocity` |
| a duration outside the representable range, or a total that underflowed to zero on a nonzero displacement | `trajectory_error::unrepresentable_duration` |

### The velocity limit is a precondition

`abs(v0) <= v_max` and `abs(v1) <= v_max` are preconditions of the type, not invariants construction restores. Raising the limit to fit a boundary velocity would return a profile that violates a bound the caller stated, which is the same silent-success defect a rejection exists to prevent; honoring the limit while accepting the command would need a ramp that runs backwards in time, since the acceleration phase spans `(v_v - v0) / a` and a cruise velocity held under the limit with `v0` above it makes that span negative. Reaching the limit exactly is inside the domain: the ramp attached to that boundary velocity simply vanishes.

### Realizability

Both ramps of a three-phase profile run toward one cruise velocity lying at or above each boundary velocity, so the profile sweeps at least the ground the transition between those two velocities already sweeps. A command below that is not realizable within this shape at any acceleration the scalar type can hold; the clearest instance is a zero commanded displacement with two boundary speeds that differ, which asks the axis to change speed while covering no ground.

Two zero-displacement commands are realizable and are not rejected: equal boundary velocities, where there is no speed change to cover, and opposed boundary velocities of equal magnitude, where the ramp between them sweeps exactly zero ground.

The realized phase durations are checked directly rather than inferred from the feasibility test that is supposed to guarantee them. Below the square root of the smallest normal value a boundary velocity squares into the subnormal range, and below the square root of the smallest subnormal it squares to exactly zero; both the feasibility test and the triangular peak are built from those squares, so the test can read as satisfied on a command that does not satisfy it and the peak can land below the boundary velocity it is analytically bounded by. The ramp duration formed from that difference is negative, and `evaluate` clamps into `[0, T]`, where a lower bound above the upper one is undefined behavior rather than an odd clamp. Checking the durations themselves closes that route without introducing a constant of its own.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `create` | `static expected<trapezoidal_trajectory, trajectory_error> create(config const&)` | The only construction path; reports an unrealizable command |
| `evaluate` | `trajectory_point<Scalar, 1> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Total duration `T = T_a + T_v + T_d` |
| `is_triangular` | `bool is_triangular() const` | True if cruise phase duration is zero |
| `peak_velocity` | `Scalar peak_velocity() const` | Signed peak velocity in original frame |
| `phase_durations` | `std::array<Scalar, 3> phase_durations() const` | `{T_accel, T_cruise, T_decel}` |
| `disposition` | `trapezoidal_disposition<Scalar> const& disposition() const` | The commanded acceleration limit against the one this profile realizes |
| `rescale_to` | `expected<void, trajectory_error> rescale_to(Scalar T_new)` | Rebuild at a longer duration for multi-axis sync |
| `can_rescale_to` | `expected<void, trajectory_error> can_rescale_to(Scalar T_new) const` | Whether `rescale_to(T_new)` would succeed, without mutating |

## Time Rescaling

`rescale_to()` rebuilds the profile at a lower cruise velocity rather than patching the one it has. The commanded displacement, both boundary velocities, and the acceleration magnitude are held fixed and the cruise velocity that realizes the requested duration is solved in closed form, so the traversed displacement and the terminal velocity hold by construction.

**The duration the profile reports is the duration it realizes.** The stored duration stays the sum of the three realized phase durations and is never assigned the requested value after the solve. When the request cannot be met the call reports a typed failure; it does not accept the request and quietly realize something else.

The solve covers all three shapes the three-phase parametrization admits, and picks between them by monotonicity: the total duration falls as the cruise velocity grows, so exactly one shape can contain the root.

| Shape | Validity | Solve |
|-------|----------|-------|
| plateau | cruise velocity at or above both boundary velocities | quadratic in the cruise velocity's **rise above** the larger boundary velocity |
| ramp-through | cruise velocity strictly between the two boundary velocities | linear in the reciprocal of the cruise velocity, single root |
| valley | cruise velocity at or below both boundary velocities | quadratic in the cruise velocity's **decrement below** the smaller boundary velocity |

The valley shape is emitted, not rejected: with both boundary velocities above the cruise velocity a long duration needs, the profile decelerates away from the initial velocity, holds a low cruise velocity, and accelerates back up to the final one. A root is accepted only inside its own shape's validity interval, with all three phase durations nonnegative and the cruise velocity within the velocity limit; a root failing any of those is a rejection rather than a clamped value.

### How accurate the realized duration is

The guarantee is two-tier.

**Exact** on the branches whose closed form inverts cleanly. A request equal to the current duration succeeds and changes nothing, which is the path the slowest axis of a synchronized set always takes, and the ramp-through shape is linear in the reciprocal of the cruise velocity with a single root and no selection.

**A derived bound, in units in the last place at the input scale**, on the two solved branches. Both are solved for the distance from their own shape boundary rather than for the cruise velocity itself, because that distance is routinely orders of magnitude smaller than the velocity it sits on. Carrying the velocity instead would discard most of the distance's significant digits before the phase durations ever saw them, and it would make the solve's conditioning proportional to the square of the boundary velocity rather than to the answer.

Measured through this API over 400,000 randomly drawn retimings per branch, spanning four decades of boundary velocity and eight of acceleration:

| Branch | Requests near the shape boundary | Requests across the shape's whole range |
|--------|----------------------------------|-----------------------------------------|
| plateau | 4 units in the last place, worst case | 27 units in the last place, worst case |
| valley | 4 units in the last place, worst case | see the caveat below |

One caveat is stated rather than hidden. As the cruise velocity approaches zero, which happens when the displacement exceeds `(v0^2 + v1^2) / (2 a)` and arbitrarily long durations become reachable, the cruise phase's duration is a residual displacement divided by that vanishing velocity. That division amplifies the residual's own rounding without bound, and it does so for any parametrization of this profile family rather than for this solve in particular. In the sweep above, six of 390,595 accepted valley retimings had a cruise velocity 25 or more orders of magnitude below the boundary velocities, and those reached 5.4e5 units in the last place. Requests that do not drive the cruise velocity to nothing are not affected.

The same caveat applies, for the same reason, to the ramp-through shape when the two boundary velocities straddle zero: its cruise velocity is the residual displacement divided by the requested duration less the fixed duration of its two ramps, so a request far longer than those ramps drives that velocity to nothing. A sweep of 1.5 million accepted retimings across twelve decades of each limit found the realized duration within a part in 1e9 of the request everywhere except there, where it reached a fifth of the request. A caller that needs a duration guarantee on this shape should keep the request within a few multiples of the ramp duration, which `phase_durations()` reports.

### Rejections

Checked in order:

| Condition | Error |
|-----------|-------|
| a duration below the current one | `trajectory_error::duration_shorter_than_current` |
| NaN, infinite, or non-positive `T_new` | `trajectory_error::non_positive_duration` |
| a request below the arithmetic resolution of the expression that would answer it | `trajectory_error::unrepresentable_duration` |
| a duration the displacement, limits, and boundary velocities cannot realize together | `trajectory_error::unreachable_duration` |

The last two are kept apart on purpose, because they tell the caller to change different things. `unreachable_duration` means no admissible shape exists for the request: change the limits or the command. `unrepresentable_duration` means a shape may well exist and the arithmetic cannot locate it: change the request.

Four forward conditions on the inputs decide the second, and all are stated against the inputs rather than by comparing a realized duration back against the request. No tolerance anywhere in the library compares the two.

- **The ramp residual must be resolved.** The residual is the commanded displacement less the distance the two ramps sweep between the boundary velocities. It collapses whenever the constructor had to raise the acceleration to make the boundary velocities feasible, because that raised value is defined by making those two equal. Both solved branches are governed by it: the plateau's whole admissible velocity range and the valley's whole admissible duration window are functions of it alone. Its floor is six chained rounding operations, at the scale of the largest operand that entered.
- **The increment or decrement against the boundary duration must be resolved.** That boundary duration is itself accurate only to some units in the last place, and the difference measured against it inherits that error. Its floor is fifteen chained rounding operations, at the scale of the largest operand that entered the boundary duration rather than at the scale of the duration alone.
- **The plateau's linear coefficient must be resolved.** That coefficient is the residual over the larger boundary velocity less the duration decrement, a difference where the valley's mirror is a sum. Both of its terms grow without bound as that boundary velocity vanishes against the commanded displacement, and deep in the plateau they agree to the full width of the significand. Its floor is twenty-three chained rounding operations, at the larger of the two terms' own scales. Only a coefficient that cannot be told apart from zero is reported here; one that is decisively negative means no such shape exists and is reported as unreachable.
- **The discriminant must be representable.** A squared linear coefficient leaves the finite range on a long move at a small cruise velocity, well before anything about the request is unreasonable. An infinite discriminant drives the root selection's denominator to infinity and its root to zero, which yields a cruise velocity sitting exactly on the shape boundary -- inside its own validity interval, with every phase duration nonnegative, and realizing the boundary's duration for a request that asked for something else.

None of the coefficients is fitted. Each is a count of the arithmetic operations chained to produce the quantity, multiplied by the scalar type's epsilon.

Reachability is derived rather than assumed, and no epsilon takes part in the decision. The cruise duration is what runs out: the shape whose validity interval reaches down toward a vanishing cruise velocity fixes the supremum of the reachable durations, and whether that supremum is finite is the sign of that shape's own residual displacement term. With both boundary velocities positive, the durations grow without bound exactly when the displacement exceeds `(v0^2 + v1^2) / (2 a)`. Below that the supremum is finite, and rather than precompute it from a square root of a difference of two nearly equal quantities, the solve decides it where it is exact: a non-positive linear coefficient means the shape's quadratic has no positive root, so no cruise velocity on the far side of the boundary answers the request, and a negative discriminant means the request lies past the shape's reachable extreme. A stationary profile therefore reaches its own duration and nothing longer.

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
    auto const created = ctrlpp::trapezoidal_trajectory<double>::create({
        .q0 = 0.0, .q1 = 10.0,
        .v_max = 5.0, .a_max = 2.0
    });
    if (!created) {
        std::cerr << "The commanded move has no trapezoidal profile\n";
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

- [double-s-trajectory](double-s-trajectory.md)<br/> jerk-limited alternative (7-segment)
- [modified-trap-trajectory](modified-trap-trajectory.md)<br/> smooth acceleration variant
- [time-scaling](time-scaling.md)<br/> duration computation for elementary paths
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
