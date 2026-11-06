# Multi-Waypoint Paths

Interpolating sequences of waypoints with splines and B-splines.

## When to Use What

| Spline Type | Endpoints | Best For |
|-------------|----------|----------|
| [Natural cubic spline](../../api/trajectory/cubic-spline.md) | Zero acceleration at both ends | General-purpose waypoint interpolation |
| [Clamped cubic spline](../../api/trajectory/cubic-spline.md) | Specified endpoint velocities | Known start/end velocity requirements |
| [Periodic cubic spline](../../api/trajectory/cubic-spline.md) | Cyclic (wraps around) | Repeating motions and closed paths |
| [Smoothing spline](../../api/trajectory/smoothing-spline.md) | Zero acceleration at both ends | Noisy waypoints needing approximation |
| [B-spline](../../api/trajectory/bspline-trajectory.md) | Configurable degree | Local control through control points |

Every spline type is scalar: it interpolates one coordinate and evaluates to a
one-dimensional `trajectory_point`. A multi-axis path is one spline per axis
over the same knot times.

## Quick Start: Cubic Spline

A natural cubic spline interpolates through all waypoints with C2
continuity. Provide time stamps and positions:

```cpp
#include <ctrlpp/trajectory/cubic_spline.h>

#include <iostream>

int main()
{
    auto const spline = ctrlpp::cubic_spline<double>::create({
        .times = {0.0, 1.0, 2.5, 4.0},
        .positions = {0.0, 1.0, 2.0, 3.0},
        .bc = ctrlpp::boundary_condition::natural,
    });
    if (!spline.has_value()) {
        return 1;
    }

    auto const pt = spline->evaluate(1.5);
    std::cout << "position at t=1.5: " << pt.position[0] << "\n";
}
```

See [example 06](../../../examples/trajectory/ctrlpp_trajectory_06_cubic_spline.cpp)
for a runnable version.

## Quick Start: B-spline

B-splines offer local control &mdash; moving one control point only affects
nearby segments. The degree is a template parameter. With an empty knot
vector `create` generates a uniform clamped one over `[0, 1]`, so the
B-spline's parameter runs over `[0, 1]` and its `duration()` is 1:

```cpp
#include <ctrlpp/trajectory/bspline_trajectory.h>

#include <iostream>

int main()
{
    auto const bspline = ctrlpp::bspline_trajectory<double, 3>::create({
        .control_points = {0.0, 1.0, 3.0, 4.0},
    });
    if (!bspline.has_value()) {
        return 1;
    }

    auto const pt = bspline->evaluate(0.5);
    std::cout << "position at t=0.5: " << pt.position[0] << "\n";
}
```

The control points shape the curve but are not passed through, apart from the
two ends of a clamped knot vector. To pass through waypoints at given times,
use `make_bspline_interpolation`. See
[example 07](../../../examples/trajectory/ctrlpp_trajectory_07_bspline.cpp)
for a runnable version.

## Smoothing Noisy Data

The [smoothing spline](../../api/trajectory/smoothing-spline.md) balances
interpolation fidelity against curve smoothness via a single parameter `mu`
in `(0, 1]`:

- `mu = 1` gives exact interpolation (passes through every waypoint)
- smaller values smooth more, and as `mu` approaches 0 the curve approaches
  the least-squares straight line through the waypoints
- `mu <= 0`, `mu > 1` and NaN are refused with `spline_error::mu_out_of_range`

The value that suits a data set depends on the noise level relative to the
signal; the library has no rule for choosing it.

## Piecewise Composition

When different stretches of a motion need different trajectory types (e.g. a
quintic, then a spline, then a cycloidal-law move), use
[piecewise_trajectory](../../api/trajectory/piecewise-trajectory.md) to place
them one after another in a single evaluable trajectory. It offsets each
segment's time by the durations of the segments before it and evaluates each
segment from local time zero. It does not join positions or velocities: each
segment's end state must match the next segment's start state, and a spline
composed this way must have its first knot at time zero.

## Limitations

### Chaining point-to-point profiles through via points is not supported

The supported multi-waypoint routes are the ones this guide shows: the cubic
and smoothing splines and the B-spline, which compute the velocity at every
waypoint themselves, and `piecewise_trajectory` over segments whose boundary
states the caller already knows. Chaining
[trapezoidal](../../api/trajectory/trapezoidal-trajectory.md) or
[double-S](../../api/trajectory/double-s-trajectory.md) profiles through via
points is not a supported workflow. The library computes no velocity at a via
point, so the join velocity is the caller's to compute and pass as `v1` of one
profile and `v0` of the next. Both profiles refuse a boundary velocity above
`v_max` with `boundary_velocity_exceeds_limit`, and the tolerance for rounding
above it, `boundary_velocity_rounding_ops`, defaults to zero: a join velocity
meant to land on `v_max` that rounds one unit in the last place above it is
refused at every displacement.

## Links

**API reference:**
[cubic spline](../../api/trajectory/cubic-spline.md) |
[smoothing spline](../../api/trajectory/smoothing-spline.md) |
[B-spline](../../api/trajectory/bspline-trajectory.md) |
[piecewise trajectory](../../api/trajectory/piecewise-trajectory.md)

**Examples:**
[06 cubic spline](../../../examples/trajectory/ctrlpp_trajectory_06_cubic_spline.cpp) |
[07 B-spline](../../../examples/trajectory/ctrlpp_trajectory_07_bspline.cpp)
