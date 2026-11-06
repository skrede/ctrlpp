# harmonic_path

Harmonic motion law: `q(tau) = (1 - cos(pi*tau)) / 2`. Produces a sinusoidal velocity profile. Non-zero acceleration at endpoints (not C2 at boundaries).

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/harmonic_path.h` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | none (`typename`) | Floating-point type |

## Functions

### harmonic_path

```cpp
template <typename Scalar>
path_point<Scalar> harmonic_path(Scalar tau);
```

Evaluate the harmonic motion law at normalized time `tau`, meant for `[0, 1]`. Returns position, velocity, acceleration, and jerk, as derivatives with respect to `tau`. `tau` is not clamped: outside `[0, 1]` the law is evaluated as written, so a caller clamps `tau` itself when it can leave the interval.

### harmonic_path_peak_derivatives

```cpp
template <typename Scalar>
std::array<Scalar, 3> harmonic_path_peak_derivatives();
```

Returns the peak magnitudes `{dq_max, ddq_max, dddq_max}` = `{pi/2, pi^2/2, pi^3/2}` over `[0, 1]`, the argument `compute_min_duration` takes.

## Peak Derivatives

| Derivative | Value | Description |
|------------|-------|-------------|
| dq_max | pi/2 ~ 1.5708 | Peak normalized velocity |
| ddq_max | pi^2/2 ~ 4.9348 | Peak normalized acceleration |
| dddq_max | pi^3/2 ~ 15.503 | Peak normalized jerk |

## Endpoint Behavior

| Boundary | q | dq | ddq |
|----------|---|-----|------|
| tau = 0 | 0 | 0 | pi^2/2 |
| tau = 1 | 1 | 0 | -pi^2/2 |

Non-zero acceleration at endpoints means this path is C1 but not C2 at segment boundaries.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include <ctrlpp/trajectory/harmonic_path.h>
#include <ctrlpp/trajectory/trajectory.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    Eigen::Vector2d q0{0.0, 0.0}, q1{1.0, 2.0};
    double T = 3.0;
    auto const traj = ctrlpp::make_trajectory(ctrlpp::harmonic_path<double>, q0, q1, T);
    if (!traj.has_value())
        return 1;

    constexpr double dt = 0.01;
    for (double t = 0.0; t <= T; t += dt) {
        auto pt = traj->evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0)
                  << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [cycloidal-path](cycloidal-path.md)<br/> zero acceleration at endpoints (C2)
- [cubic-path](cubic-path.md)<br/> polynomial alternative (C1)
- [time-scaling](time-scaling.md)<br/> uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
