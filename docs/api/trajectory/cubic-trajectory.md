# cubic_trajectory

Cubic polynomial trajectory segment with arbitrary velocity boundary conditions. Interpolates between two points with specified endpoint velocities using a degree-3 polynomial in normalized time `tau = t/T`.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/cubic_trajectory.h` |
| **Factory** | `ctrlpp::make_cubic_trajectory` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | none (`typename`) | Floating-point type; the factory tests it with `std::isfinite` |
| `ND` | `std::size_t` | Number of dimensions, deduced by the factory from the rows of its fixed-size vector arguments |

## Struct Members

| Field | Type | Description |
|-------|------|-------------|
| `coeffs` | `std::array<Vector<Scalar, ND>, 4>` | Normalized polynomial coefficients `[c0, c1, c2, c3]` |
| `dur` | `Scalar` | Segment duration, in the caller's time unit |

`cubic_trajectory` is a public aggregate. Built directly rather than through the factory, nothing checks it: `evaluate` clamps `t` to `[0, dur]` and divides by `dur`, so `dur` must be finite and strictly positive.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, ND> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Segment duration `T` |

Evaluation uses Horner's method for numerical stability. Time is clamped to `[0, T]`.

## Factory Function

```cpp
template <typename Scalar, int Rows>
ctrlpp::expected<cubic_trajectory<Scalar, static_cast<std::size_t>(Rows)>, trajectory_error>
make_cubic_trajectory(
    Eigen::Matrix<Scalar, Rows, 1> const& q0,   // start position
    Eigen::Matrix<Scalar, Rows, 1> const& q1,   // end position
    Eigen::Matrix<Scalar, Rows, 1> const& v0,   // start velocity
    Eigen::Matrix<Scalar, Rows, 1> const& v1,   // end velocity
    Scalar duration
);
```

The coefficients are those of the polynomial in normalized time `tau = t/T`: with `h = q1 - q0`, `c0 = q0`, `c1 = v0*T`, `c2 = 3*h - T*(2*v0 + v1)`, `c3 = -2*h + T*(v0 + v1)`.

The duration is a divisor (`evaluate` scales derivatives by `1/T`), so its exact domain is finite and strictly positive. Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| NaN/Inf duration or any non-finite boundary entry | `trajectory_error::non_finite_input` |
| `duration <= 0` | `trajectory_error::non_positive_duration` |

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include "ctrlpp/trajectory/cubic_trajectory.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    auto result = ctrlpp::make_cubic_trajectory(
        Vec1{0.0}, Vec1{1.0},   // q0, q1
        Vec1{0.0}, Vec1{0.0},   // v0, v1 (rest-to-rest)
        2.0);                    // duration
    if (!result.has_value())
        return 1;
    auto const& traj = *result;
    for (double t = 0; t <= traj.duration(); t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [cubic-path](cubic-path.md)<br/> normalized cubic (zero-velocity BCs only)
- [quintic-trajectory](quintic-trajectory.md)<br/> adds acceleration BCs
- [piecewise-trajectory](piecewise-trajectory.md)<br/> multi-segment composition
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
