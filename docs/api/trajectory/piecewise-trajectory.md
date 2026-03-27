# piecewise_trajectory

Variadic composition of heterogeneous trajectory segments into a single trajectory. Segments are evaluated in sequence with automatic time offset at breakpoints. The composed trajectory itself satisfies `trajectory_segment`, enabling recursive composition.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/piecewise_trajectory.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |
| `ND` | Number of dimensions (`std::size_t`) |
| `Segments...` | Pack of types, each satisfying `trajectory_segment<Scalar, ND>` |

## Construction

```cpp
piecewise_trajectory(Segments... segs);
```

Segments are moved into internal storage. Breakpoints are computed from cumulative `duration()` calls.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, ND> evaluate(Scalar t) const` | Evaluate at time `t`, dispatched to correct segment |
| `duration` | `Scalar duration() const` | Sum of all segment durations |

## Static Members

| Member | Type | Description |
|--------|------|-------------|
| `segment_count` | `std::size_t` | Number of segments (compile-time) |

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include "ctrlpp/traj/piecewise_trajectory.h"
#include "ctrlpp/traj/cubic_trajectory.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    auto seg1 = ctrlpp::make_cubic_trajectory(Vec1{0.0}, Vec1{1.0}, Vec1{0.0}, Vec1{0.5}, 1.0);
    auto seg2 = ctrlpp::make_cubic_trajectory(Vec1{1.0}, Vec1{3.0}, Vec1{0.5}, Vec1{0.0}, 1.5);
    ctrlpp::piecewise_trajectory<double, 1, decltype(seg1), decltype(seg2)> pw(seg1, seg2);
    for (double t = 0; t <= pw.duration(); t += 0.01) {
        auto pt = pw.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [piecewise-path](piecewise-path.md) -- analogous composition for normalized paths
- [trajectory-segment](trajectory-segment.md) -- concept each segment must satisfy
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
