# smoothing_spline

Smoothing spline approximation with configurable mu tradeoff parameter. Constructs a C2-continuous spline that balances data fidelity against smoothness. The mu parameter controls the tradeoff: mu=1 yields exact interpolation (passes through all waypoints), while mu near 0 maximizes smoothness at the cost of data fit.

## Header

```cpp
#include "ctrlpp/trajectory/smoothing_spline.h"
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (`float`, `double`, `long double`) |

## Config

```cpp
struct config {
    std::vector<Scalar> times;       // Knot times t_0 ... t_n (n+1 entries)
    std::vector<Scalar> positions;   // Waypoint positions q_0 ... q_n (n+1 entries)
    Scalar mu{0.5};                  // Tradeoff in (0, 1]: 1 = interpolation, near 0 = max smoothness
};
```

## The mu Parameter

The `mu` parameter maps to an internal regularization weight lambda = 2(1-mu) / (3*mu):

| mu | Behavior |
|----|-----------|
| 1.0 | Exact interpolation (lambda = 0, passes through all waypoints) |
| 0.5 | Balanced smoothness and data fidelity (lambda = 2/3) |
| near 0 | Maximum smoothness (nearly straight line, ignores data) |

The domain of `mu` is the half-open interval (0, 1]: mu = 1 is the exact interpolation limit and lambda diverges as mu approaches 0, so mu <= 0 has no defined weight. Values outside the domain (including NaN) are rejected by `create` with `spline_error::mu_out_of_range`.

## Factory

```cpp
[[nodiscard]] static auto create(config const& cfg)
    -> ctrlpp::expected<smoothing_spline, spline_error>;
```

`create` is the only construction path. There is no non-fallible constructor: a rejected configuration is a value the caller has to inspect, never an object that quietly stands in for one. It validates the configuration and constructs a smoothing spline. Solves the regularized system (R + lambda * Q^T * Q) * d = Q^T * q for interior second derivatives using dense QR factorization. Requires at least 2 waypoints. For 2 waypoints, degenerates to a linear segment.

Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| Fewer than 2 waypoints | `spline_error::too_few_points` |
| `times` and `positions` differ in length | `spline_error::size_mismatch` |
| Knot times not strictly increasing | `spline_error::non_increasing_times` |
| `mu` outside (0, 1] or NaN | `spline_error::mu_out_of_range` |

## Methods

### evaluate

```cpp
auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>;
```

Evaluate position, velocity, and acceleration at time `t`. Time is clamped to [t_0, t_n].

### duration

```cpp
auto duration() const -> Scalar;
```

Returns total spline duration: t_n - t_0.

## Concept Satisfaction

`smoothing_spline<Scalar>` satisfies `trajectory_segment<smoothing_spline<Scalar>, Scalar, 1>`.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include <ctrlpp/trajectory/smoothing_spline.h>

#include <iostream>

int main()
{
    // Noisy waypoints, smoothing removes noise while preserving shape
    auto spline = ctrlpp::smoothing_spline<double>::create({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.1, 0.4, 1.6, 2.0},  // noisy measurements
        .mu = 0.7,  // moderate smoothing
    });
    if (!spline.has_value()) {
        return 1;
    }

    double T = spline->duration();
    constexpr double dt = 0.01;
    for (double t = 0.0; t <= T; t += dt) {
        auto pt = spline->evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [cubic-spline](cubic-spline.md)<br/> Exact interpolation with natural, clamped, or periodic BCs
- [trajectory-types](trajectory-types.md)<br/> `spline_error` enumerators returned by `create`
- [Trajectory Generation Theory](../../background/trajectory-generation.md)<br/> Smoothing spline regularization formulation
