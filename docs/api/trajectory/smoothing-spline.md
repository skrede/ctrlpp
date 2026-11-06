# smoothing_spline

Smoothing spline approximation of scalar waypoints with configurable mu tradeoff parameter. Constructs a C2-continuous spline, with zero acceleration at both ends, that balances data fidelity against smoothness. The mu parameter controls the tradeoff: mu=1 yields exact interpolation (passes through all waypoints), while mu near 0 maximizes smoothness at the cost of data fit.

## Header

```cpp
#include "ctrlpp/trajectory/smoothing_spline.h"
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type (`float`, `double`, `long double`) |

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
| near 0 | Maximum smoothness (approaches the least-squares straight line through the data) |

The domain of `mu` is the half-open interval (0, 1]: mu = 1 is the exact interpolation limit and lambda diverges as mu approaches 0, so mu <= 0 has no defined weight. Values outside the domain (including NaN) are rejected by `create` with `spline_error::mu_out_of_range`.

## Factory

```cpp
static ctrlpp::expected<smoothing_spline, spline_error> create(config const& cfg);
```

`create` is the only construction path. There is no non-fallible constructor: a rejected configuration is a value the caller has to inspect, never an object that quietly stands in for one. It validates the configuration and constructs a smoothing spline. Solves the regularized system (R + lambda * Q^T * Q) * d = Q^T * q for interior second derivatives using dense QR factorization. Requires at least 2 waypoints. With 2 waypoints there is no interior knot, and the spline is the line through them whatever `mu` is.

Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| Fewer than 2 waypoints | `spline_error::too_few_points` |
| `times` and `positions` differ in length | `spline_error::size_mismatch` |
| A non-finite knot time or waypoint | `spline_error::non_finite_input` |
| Knot times not strictly increasing | `spline_error::non_increasing_times` |
| A knot span `t_(i+1) - t_i` that overflows `Scalar` (checked span by span with the row above) | `spline_error::unrepresentable_spline` |
| `mu` outside (0, 1] or NaN | `spline_error::mu_out_of_range` |
| Regularized system, or resulting coefficients, outside the representable range of `Scalar` | `spline_error::unrepresentable_spline` |

A `mu` inside its domain can still ask for a weight the arithmetic cannot carry, and it fails two different ways. The system matrix `R + lambda * Q^T * Q` is solved by a decomposition that forms sums of squares of its entries, so an entry whose square overflows turns the pivots into infinities and every coefficient into a NaN. On the smallest system -- three waypoints, a single interior knot -- the same weight instead overflows that one entry, the solve divides by the infinity, the interior second derivatives come back as exact zeros, and the smoothed positions collapse onto the raw waypoints: a finite, plausible straight-line-through-the-data answer that is not the least-squares limit the weight asked for, and that nothing downstream can distinguish from a correct one. The entries are therefore bounded before the solve, by the square root of the largest representable value, which is the condition that the squaring survives. Both terms of the sum are counted, because each reaches that bound on its own and from opposite directions: the regularization term through a small `mu` or a small knot spacing, and the smoothness term `R` through a large knot spacing. At a spacing of 1e200 `R`'s diagonal is 6.7e199, whose square is not representable, and the solve returns NaN for every coefficient on data as ordinary as a straight line. `Q^T * Q` is a Gram matrix and `R` is diagonally dominant, so neither exceeds its own largest diagonal entry, and both diagonals follow from the knot spacing and `mu` alone.

The band of knot spacings that bound admits is therefore bounded at both ends, and for `mu = 0.5` and equally spaced knots it runs from about 2e-77 to about 2e154. A **constant-velocity segment is accepted across the whole of it**, and down to about 4e-125 for the interpolating spline (`mu = 1`). Curved data can still be refused inside the band, by the coefficient check: the waypoints `0, 1, 0.5, 2` at `mu = 0.5` are refused above a spacing of about 3e102. The constant-velocity segment survives there because the cubic coefficient's numerator is a difference of second derivatives, so on collinear waypoints it is mathematically zero and arrives as the rounding residual of that cancellation rather than as an exact zero. That residual reaches a small fraction of the resolution the evaluated position is known at, so the term cannot move the answer and the coefficient's precision does not matter. What decides representability is how far a term reaches across its own span against that resolution, not the magnitude of any intermediate on its own.

## Methods

### evaluate

```cpp
trajectory_point<Scalar, 1> evaluate(Scalar t) const;
```

Evaluate position, velocity, and acceleration at time `t`. Time is clamped to [t_0, t_n]. A time exactly on an interior knot `t_i` is evaluated by the span that starts there; `t_n` is evaluated by the last span.

### duration

```cpp
Scalar duration() const;
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
