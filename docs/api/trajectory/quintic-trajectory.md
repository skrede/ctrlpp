# quintic_trajectory

Quintic polynomial trajectory segment with arbitrary velocity and acceleration boundary conditions. Interpolates between two points using a degree-5 polynomial in normalized time `tau = t/T`.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/quintic_trajectory.h` |
| **Factory** | `ctrlpp::make_quintic_trajectory` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |
| `ND` | Number of dimensions (`std::size_t`) |

## Struct Members

| Field | Type | Description |
|-------|------|-------------|
| `coeffs` | `std::array<Vector<Scalar, ND>, 6>` | Normalized polynomial coefficients `[c0..c5]` |
| `dur` | `Scalar` | Segment duration in seconds |

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, ND> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Segment duration `T` |

## Factory Function

```cpp
template <typename Scalar, int Rows>
auto make_quintic_trajectory(
    Eigen::Matrix<Scalar, Rows, 1> const& q0,   // start position
    Eigen::Matrix<Scalar, Rows, 1> const& q1,   // end position
    Eigen::Matrix<Scalar, Rows, 1> const& v0,   // start velocity
    Eigen::Matrix<Scalar, Rows, 1> const& v1,   // end velocity
    Eigen::Matrix<Scalar, Rows, 1> const& a0,   // start acceleration
    Eigen::Matrix<Scalar, Rows, 1> const& a1,   // end acceleration
    Scalar duration
) -> quintic_trajectory<Scalar, ND>;
```

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include <ctrlpp/traj/quintic_trajectory.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    double T = 2.0;
    auto traj = ctrlpp::make_quintic_trajectory(
        Vec1{0.0}, Vec1{1.0},   // q0, q1
        Vec1{0.0}, Vec1{0.0},   // v0, v1
        Vec1{0.0}, Vec1{0.0},   // a0, a1
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

- [quintic-path](quintic-path.md) -- normalized quintic (zero BCs only)
- [cubic-trajectory](cubic-trajectory.md) -- lower degree (velocity BCs only)
- [septic-trajectory](septic-trajectory.md) -- higher degree (adds jerk BCs)
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
