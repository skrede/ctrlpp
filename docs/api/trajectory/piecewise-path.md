# piecewise_path

Variadic composition of heterogeneous path segments into a single path. Segments are evaluated in sequence with automatic time offset at breakpoints. The composed path itself satisfies `path_segment`, enabling recursive composition.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/piecewise_path.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |
| `Segments...` | Pack of types, each satisfying `path_segment<Scalar>` |

## Construction

```cpp
piecewise_path(Segments... segs);
```

Segments are moved into internal storage. Breakpoints are computed from cumulative `duration()` calls.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `path_point<Scalar> evaluate(Scalar t) const` | Evaluate at time `t`, dispatched to correct segment |
| `duration` | `Scalar duration() const` | Sum of all segment durations |

## Static Members

| Member | Type | Description |
|--------|------|-------------|
| `segment_count` | `std::size_t` | Number of segments (compile-time) |

## Dispatch Logic

Time `t` is clamped to `[0, total_duration]`. The correct segment is found by comparing `t` against cumulative breakpoints. Local time is computed as `t - breakpoint[i-1]`.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include "ctrlpp/trajectory/piecewise_trajectory.h"
#include "ctrlpp/trajectory/cubic_trajectory.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    // piecewise_path requires types satisfying path_segment concept.
    // Elementary paths are free functions, so demonstrate piecewise_trajectory instead.
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    auto seg1 = ctrlpp::make_cubic_trajectory(Vec1{0.0}, Vec1{0.5}, Vec1{0.0}, Vec1{0.3}, 1.0);
    auto seg2 = ctrlpp::make_cubic_trajectory(Vec1{0.5}, Vec1{1.0}, Vec1{0.3}, Vec1{0.0}, 1.0);
    ctrlpp::piecewise_trajectory<double, 1, decltype(seg1), decltype(seg2)> pw(seg1, seg2);
    for (double t = 0; t <= pw.duration(); t += 0.005) {
        auto pt = pw.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [piecewise-trajectory](piecewise-trajectory.md) -- analogous composition for trajectories
- [path-segment](path-segment.md) -- concept each segment must satisfy
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
