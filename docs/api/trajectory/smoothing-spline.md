# smoothing_spline

Smoothing spline approximation with configurable mu tradeoff parameter. Constructs a C2-continuous spline that balances data fidelity against smoothness. The mu parameter controls the tradeoff: mu=1 yields exact interpolation (passes through all waypoints), while mu near 0 maximises smoothness at the cost of data fit.

## Header

```cpp
#include "ctrlpp/traj/smoothing_spline.h"
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
    Scalar mu{0.5};                  // Tradeoff: 1.0 = interpolation, near 0 = max smoothness
};
```

## The mu Parameter

The `mu` parameter maps to an internal regularisation weight lambda = 2(1-mu) / (3*mu):

| mu | Behaviour |
|----|-----------|
| 1.0 | Exact interpolation (lambda = 0, passes through all waypoints) |
| 0.5 | Balanced smoothness and data fidelity (lambda = 2/3) |
| ~0 | Maximum smoothness (nearly straight line, ignores data) |

The parameter is clamped to (epsilon, 1.0] internally to avoid degenerate lambda values.

## Constructor

```cpp
explicit smoothing_spline(config const& cfg);
```

Constructs a smoothing spline from waypoints and mu parameter. Solves the regularised system (R + lambda * Q^T * Q) * d = Q^T * q for interior second derivatives using dense QR factorisation. Requires at least 2 waypoints. For 2 waypoints, degenerates to a linear segment.

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

#include <ctrlpp/traj/smoothing_spline.h>

#include <iostream>

int main()
{
    // Noisy waypoints -- smoothing removes noise while preserving shape
    ctrlpp::smoothing_spline<double> spline({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.1, 0.4, 1.6, 2.0},  // noisy measurements
        .mu = 0.7,  // moderate smoothing
    });

    double T = spline.duration();
    constexpr double dt = 0.01;
    for (double t = 0.0; t <= T; t += dt) {
        auto pt = spline.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [cubic-spline](cubic-spline.md) -- Exact interpolation with natural, clamped, or periodic BCs
- [Trajectory Generation Theory](../../background/trajectory-generation.md) -- Smoothing spline regularisation formulation
