# trajectory_segment

Concept constraining types that represent a trajectory segment with timed evaluation. A trajectory segment provides `evaluate(t)` returning a `trajectory_point` and `duration()` returning the total time span.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/trajectory_segment.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `S` | The type being constrained |
| `Scalar` | Floating-point type |
| `ND` | Number of dimensions (`std::size_t`) |

## Concept Definition

```cpp
template <typename S, typename Scalar, std::size_t ND>
concept trajectory_segment = requires(const S& seg, Scalar t) {
    { seg.evaluate(t) } -> std::convertible_to<trajectory_point<Scalar, ND>>;
    { seg.duration() } -> std::convertible_to<Scalar>;
};
```

## Required Expressions

| Expression | Return type | Description |
|------------|-------------|-------------|
| `seg.evaluate(t)` | convertible to `trajectory_point<Scalar, ND>` | Position, velocity, acceleration at time `t` |
| `seg.duration()` | convertible to `Scalar` | Total segment duration in seconds |

## Satisfying Types

All trajectory types in `ctrlpp/traj/` satisfy this concept:

- `trajectory<Law, Scalar, ND>`
- `cubic_trajectory<Scalar, ND>`
- `quintic_trajectory<Scalar, ND>`
- `septic_trajectory<Scalar, ND>`
- `trapezoidal_trajectory<Scalar>` (ND=1)
- `double_s_trajectory<Scalar>` (ND=1)
- `modified_trap_trajectory<Scalar>` (ND=1)
- `modified_sin_trajectory<Scalar>` (ND=1)
- `piecewise_trajectory<Scalar, ND, Segments...>`

## Usage Example

```cpp
#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/cubic_trajectory.h"

template <ctrlpp::trajectory_segment<double, 1> Traj>
void sample(Traj const& traj, double dt) {
    for (double t = 0; t <= traj.duration(); t += dt) {
        auto pt = traj.evaluate(t);
    }
}
```

## See Also

- [trajectory-types](trajectory-types.md) -- `trajectory_point` struct
- [path-segment](path-segment.md) -- analogous concept for normalized paths
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
