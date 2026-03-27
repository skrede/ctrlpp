# piecewise_path

Variadic composition of heterogeneous path segments into a single path. Segments are evaluated in sequence with automatic time offset at breakpoints. The composed path itself satisfies `path_segment`, enabling recursive composition.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/piecewise_path.h` |

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
#include "ctrlpp/traj/piecewise_path.h"
#include "ctrlpp/traj/cubic_path.h"

// Compose two path segments (requires types satisfying path_segment concept)
// For free-function paths, use the trajectory adapter instead.
```

## See Also

- [piecewise-trajectory](piecewise-trajectory.md) -- analogous composition for trajectories
- [path-segment](path-segment.md) -- concept each segment must satisfy
- [Trajectory Generation Theory](../../reference/trajectory-generation.md)
