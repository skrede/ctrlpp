# trajectory_segment

Concept constraining types that represent a trajectory segment with timed evaluation. A trajectory segment provides `evaluate(t)` returning a `trajectory_point` and `duration()` returning the total time span.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/trajectory_segment.h` |

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

All trajectory types in `ctrlpp/trajectory/` satisfy this concept:

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
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include "ctrlpp/trajectory/cubic_trajectory.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    auto seg = ctrlpp::make_cubic_trajectory(Vec1{0.0}, Vec1{1.0}, Vec1{0.0}, Vec1{0.0}, 2.0);
    for (double t = 0; t <= seg.duration(); t += 0.01) {
        auto pt = seg.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [trajectory-types](trajectory-types.md) -- `trajectory_point` struct
- [path-segment](path-segment.md) -- analogous concept for normalized paths
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
