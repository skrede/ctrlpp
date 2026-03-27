# time_scaling

Kinematic time scaling for normalized paths. Computes the minimum trajectory duration T from displacement magnitude, peak normalized path derivatives, and kinematic limits (v_max, a_max, j_max). For use with elementary paths only -- composite profiles (trapezoidal, double-S) solve timing internally.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/time_scaling.h` |

## Functions

### compute_min_duration (scalar)

```cpp
template <typename Scalar>
auto compute_min_duration(
    Scalar h,
    std::array<Scalar, 3> const& peak_derivs,
    Scalar v_max, Scalar a_max, Scalar j_max) -> Scalar;
```

Computes minimum duration for a single DOF:

```
T = max(h*dq_max/v_max, sqrt(h*ddq_max/a_max), cbrt(h*dddq_max/j_max))
```

| Parameter | Description |
|-----------|-------------|
| `h` | Displacement magnitude (\|q1 - q0\|) |
| `peak_derivs` | `{dq_max, ddq_max, dddq_max}` from path's `*_peak_derivatives()` |
| `v_max` | Maximum velocity limit |
| `a_max` | Maximum acceleration limit |
| `j_max` | Maximum jerk limit |

### compute_min_duration (vector)

```cpp
template <typename Scalar, std::size_t N>
auto compute_min_duration(
    std::array<Scalar, N> const& h,
    std::array<Scalar, 3> const& peak_derivs,
    Scalar v_max, Scalar a_max, Scalar j_max) -> std::array<Scalar, N>;
```

Returns per-DOF minimum durations. Each DOF is independently computed.

### compute_min_duration_sync

```cpp
template <typename Scalar, std::size_t N>
auto compute_min_duration_sync(
    std::array<Scalar, N> const& h,
    std::array<Scalar, 3> const& peak_derivs,
    Scalar v_max, Scalar a_max, Scalar j_max) -> Scalar;
```

Returns a single T (maximum of per-DOF durations) suitable for synchronizing all joints.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include "ctrlpp/trajectory/time_scaling.h"
#include "ctrlpp/trajectory/quintic_path.h"
#include "ctrlpp/trajectory/trajectory.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    auto peaks = ctrlpp::quintic_path_peak_derivatives<double>();
    double T = ctrlpp::compute_min_duration(5.0, peaks, 2.0, 10.0, 100.0);

    Eigen::Matrix<double, 1, 1> q0{0.0}, q1{5.0};
    auto traj = ctrlpp::make_trajectory(ctrlpp::quintic_path<double>, q0, q1, T);
    for (double t = 0; t <= traj.duration(); t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [cubic-path](cubic-path.md), [quintic-path](quintic-path.md), [septic-path](septic-path.md) -- paths providing peak derivatives
- [harmonic-path](harmonic-path.md), [cycloidal-path](cycloidal-path.md) -- trigonometric paths
- [trajectory](trajectory.md) -- adapter using computed duration
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
