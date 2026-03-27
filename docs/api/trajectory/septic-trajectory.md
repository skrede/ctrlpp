# septic_trajectory

Septic polynomial trajectory segment with arbitrary boundary conditions up to jerk. Interpolates between two points using a degree-7 polynomial in normalized time `tau = t/T`.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/septic_trajectory.h` |
| **Factory** | `ctrlpp::make_septic_trajectory` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |
| `ND` | Number of dimensions (`std::size_t`) |

## Struct Members

| Field | Type | Description |
|-------|------|-------------|
| `coeffs` | `std::array<Vector<Scalar, ND>, 8>` | Normalized polynomial coefficients `[c0..c7]` |
| `dur` | `Scalar` | Segment duration in seconds |

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, ND> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Segment duration `T` |

## Factory Function

```cpp
template <typename Scalar, int Rows>
auto make_septic_trajectory(
    Eigen::Matrix<Scalar, Rows, 1> const& q0,   // start position
    Eigen::Matrix<Scalar, Rows, 1> const& q1,   // end position
    Eigen::Matrix<Scalar, Rows, 1> const& v0,   // start velocity
    Eigen::Matrix<Scalar, Rows, 1> const& v1,   // end velocity
    Eigen::Matrix<Scalar, Rows, 1> const& a0,   // start acceleration
    Eigen::Matrix<Scalar, Rows, 1> const& a1,   // end acceleration
    Eigen::Matrix<Scalar, Rows, 1> const& j0,   // start jerk
    Eigen::Matrix<Scalar, Rows, 1> const& j1,   // end jerk
    Scalar duration
) -> septic_trajectory<Scalar, ND>;
```

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include <ctrlpp/trajectory/septic_trajectory.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    Vec1 zero{0.0};
    double T = 2.0;
    auto traj = ctrlpp::make_septic_trajectory(
        Vec1{0.0}, Vec1{1.0},   // q0, q1
        zero, zero,              // v0, v1
        zero, zero,              // a0, a1
        zero, zero,              // j0, j1
        T);                      // duration

    constexpr double dt = 0.01;
    for (double t = 0.0; t <= T; t += dt) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0)
                  << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [septic-path](septic-path.md) -- normalized septic (zero BCs only)
- [quintic-trajectory](quintic-trajectory.md) -- lower degree (no jerk BCs)
- [piecewise-trajectory](piecewise-trajectory.md) -- multi-segment composition
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
