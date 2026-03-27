# cubic_trajectory

Cubic polynomial trajectory segment with arbitrary velocity boundary conditions. Interpolates between two points with specified endpoint velocities using a degree-3 polynomial in normalized time `tau = t/T`.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/cubic_trajectory.h` |
| **Factory** | `ctrlpp::make_cubic_trajectory` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |
| `ND` | Number of dimensions (`std::size_t`) |

## Struct Members

| Field | Type | Description |
|-------|------|-------------|
| `coeffs` | `std::array<Vector<Scalar, ND>, 4>` | Normalized polynomial coefficients `[c0, c1, c2, c3]` |
| `dur` | `Scalar` | Segment duration in seconds |

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, ND> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Segment duration `T` |

Evaluation uses Horner's method for numerical stability. Time is clamped to `[0, T]`.

## Factory Function

```cpp
template <typename Scalar, int Rows>
auto make_cubic_trajectory(
    Eigen::Matrix<Scalar, Rows, 1> const& q0,   // start position
    Eigen::Matrix<Scalar, Rows, 1> const& q1,   // end position
    Eigen::Matrix<Scalar, Rows, 1> const& v0,   // start velocity
    Eigen::Matrix<Scalar, Rows, 1> const& v1,   // end velocity
    Scalar duration
) -> cubic_trajectory<Scalar, ND>;
```

Coefficients are computed in normalized time to prevent ill-conditioning.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include "ctrlpp/traj/cubic_trajectory.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    auto traj = ctrlpp::make_cubic_trajectory(
        Vec1{0.0}, Vec1{1.0},   // q0, q1
        Vec1{0.0}, Vec1{0.0},   // v0, v1 (rest-to-rest)
        2.0);                    // duration
    for (double t = 0; t <= traj.duration(); t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [cubic-path](cubic-path.md) -- normalized cubic (zero-velocity BCs only)
- [quintic-trajectory](quintic-trajectory.md) -- adds acceleration BCs
- [piecewise-trajectory](piecewise-trajectory.md) -- multi-segment composition
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
