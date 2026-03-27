# cycloidal_path

Cycloidal motion law: `q(tau) = tau - sin(2*pi*tau) / (2*pi)`. Zero acceleration at endpoints with a smooth sinusoidal acceleration profile.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/cycloidal_path.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

## Functions

### cycloidal_path

```cpp
template <typename Scalar>
auto cycloidal_path(Scalar tau) -> path_point<Scalar>;
```

Evaluate the cycloidal motion law at normalized time `tau` in [0,1].

### cycloidal_path_peak_derivatives

```cpp
template <typename Scalar>
auto cycloidal_path_peak_derivatives() -> std::array<Scalar, 3>;
```

Returns `{dq_max, ddq_max, dddq_max}` = `{2, 2*pi, 4*pi^2}`.

## Peak Derivatives

| Derivative | Value | Description |
|------------|-------|-------------|
| dq_max | 2.0 | Peak normalized velocity |
| ddq_max | 2*pi ~ 6.2832 | Peak normalized acceleration |
| dddq_max | 4*pi^2 ~ 39.478 | Peak normalized jerk |

## Endpoint Behavior

| Boundary | q | dq | ddq |
|----------|---|-----|------|
| tau = 0 | 0 | 0 | 0 |
| tau = 1 | 1 | 0 | 0 |

Zero acceleration at both endpoints makes this path C2-continuous at segment boundaries.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include <ctrlpp/traj/cycloidal_path.h>
#include <ctrlpp/traj/time_scaling.h>
#include <ctrlpp/traj/trajectory.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    auto peaks = ctrlpp::cycloidal_path_peak_derivatives<double>();
    double T = ctrlpp::compute_min_duration(2.0, peaks, 1.0, 5.0, 50.0);

    Eigen::Vector2d q0{0.0, 0.0}, q1{2.0, 0.0};
    auto traj = ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q1, T);

    constexpr double dt = 0.01;
    for (double t = 0.0; t <= T; t += dt) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0)
                  << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [harmonic-path](harmonic-path.md) -- sinusoidal velocity (non-zero endpoint acceleration)
- [quintic-path](quintic-path.md) -- polynomial alternative (C2)
- [time-scaling](time-scaling.md) -- uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
