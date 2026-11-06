# piecewise_trajectory

Variadic composition of heterogeneous trajectory segments into a single trajectory. Segments are evaluated in sequence with automatic time offsets at the breakpoints. The composed trajectory itself satisfies `trajectory_segment`, enabling recursive composition.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/piecewise_trajectory.h` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | none (`typename`) | The time and value type |
| `ND` | `std::size_t` | Number of dimensions |
| `Segments...` | `trajectory_segment<Segment, Scalar, ND>` for each | The segment types, in evaluation order |

## Construction

```cpp
explicit piecewise_trajectory(Segments... segs);
```

Segments are taken by value and moved into internal storage. Breakpoints are the cumulative sums of the segments' `duration()`, computed once here. Nothing is validated: each duration must be finite and nonnegative.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, ND> evaluate(Scalar t) const` | Evaluate at time `t`, dispatched to the owning segment |
| `duration` | `Scalar duration() const` | Sum of all segment durations |

## Static Members

| Member | Type | Description |
|--------|------|-------------|
| `segment_count` | `std::size_t` | Number of segments (compile-time) |

## Dispatch Logic

`t` is clamped to `[0, duration()]`. Segment `i` owns `[b_(i-1), b_i)`, where `b_i` is the sum of the first `i + 1` durations and `b_(-1) = 0`, and is evaluated at the local time `t - b_(i-1)`. A time exactly on a breakpoint belongs to the later segment, at its local time zero; the last segment also owns the end, `t = duration()`. A zero-duration segment other than the last is never evaluated.

Every segment is evaluated from local time zero. A [`cubic_spline`](cubic-spline.md) or [`smoothing_spline`](smoothing-spline.md) evaluates over its own knot interval `[t_0, t_n]` and clamps to it, so as a segment it must have `t_0 = 0`; with `t_0 > 0` its start is held for the first `t_0` of its slot and its last `t_0` is never reached. The same holds for a [`bspline_trajectory`](bspline-trajectory.md) whose knot vector does not start at zero.

Only time is offset. The values are each segment's own: `piecewise_trajectory` does not shift positions or match velocities at a breakpoint, so a continuous composition is one whose segments end where the next begins, as in the example below.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include "ctrlpp/trajectory/cubic_trajectory.h"
#include "ctrlpp/trajectory/piecewise_trajectory.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    auto const seg1 = ctrlpp::make_cubic_trajectory(Vec1{0.0}, Vec1{1.0}, Vec1{0.0}, Vec1{0.5}, 1.0);
    auto const seg2 = ctrlpp::make_cubic_trajectory(Vec1{1.0}, Vec1{3.0}, Vec1{0.5}, Vec1{0.0}, 1.5);
    if (!seg1.has_value() || !seg2.has_value())
        return 1;

    using segment = ctrlpp::cubic_trajectory<double, 1>;
    ctrlpp::piecewise_trajectory<double, 1, segment, segment> pw{*seg1, *seg2};
    for (double t = 0; t <= pw.duration(); t += 0.01) {
        auto pt = pw.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [piecewise-path](piecewise-path.md)<br/> analogous composition for normalized paths
- [trajectory-segment](trajectory-segment.md)<br/> concept each segment must satisfy
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
