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
#include "ctrlpp/traj/piecewise_trajectory.h"
#include "ctrlpp/traj/cubic_trajectory.h"

using Vec1 = Eigen::Matrix<double, 1, 1>;
Vec1 q0{0.0}, q1{1.0}, q2{3.0}, v0{0.0}, v1{0.5}, v2{0.0};

auto seg1 = ctrlpp::make_cubic_trajectory(q0, q1, v0, v1, 1.0);
auto seg2 = ctrlpp::make_cubic_trajectory(q1, q2, v1, v2, 1.5);

ctrlpp::piecewise_trajectory<double, 1, decltype(seg1), decltype(seg2)> pw(seg1, seg2);
auto pt = pw.evaluate(1.5);  // in second segment
// pw.duration() == 2.5
```

## See Also

- [piecewise-path](piecewise-path.md) -- analogous composition for normalized paths
- [trajectory-segment](trajectory-segment.md) -- concept each segment must satisfy
- [Trajectory Generation Theory](../../reference/trajectory-generation.md)
