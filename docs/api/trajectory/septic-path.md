# septic_path

Normalized septic polynomial path: `q(tau) = 35*tau^4 - 84*tau^5 + 70*tau^6 - 20*tau^7`. C3-continuous with zero velocity, acceleration, and jerk at endpoints.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/septic_path.h` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | none (`typename`) | Floating-point type |

## Functions

### septic_path

```cpp
template <typename Scalar>
path_point<Scalar> septic_path(Scalar tau);
```

Evaluate the septic motion law at normalized time `tau`, meant for `[0, 1]`. Returns position, velocity, acceleration, and jerk, as derivatives with respect to `tau`. `tau` is not clamped: outside `[0, 1]` the polynomial is evaluated as written, so a caller clamps `tau` itself when it can leave the interval.

### septic_path_peak_derivatives

```cpp
template <typename Scalar>
std::array<Scalar, 3> septic_path_peak_derivatives();
```

Returns the peak magnitudes `{dq_max, ddq_max, dddq_max}` = `{35/16, 7.5132, 52.5}` over `[0, 1]`, the argument `compute_min_duration` takes. `ddq_max` is the literal `7.5132`, the peak rounded to five significant figures; its exact value has no simple closed form.

## Peak Derivatives

| Derivative | Value | Description |
|------------|-------|-------------|
| dq_max | 2.1875 | Peak normalized velocity (35/16) |
| ddq_max | 7.5132 | Peak normalized acceleration (rounded literal) |
| dddq_max | 52.5 | Peak normalized jerk |

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include <ctrlpp/trajectory/septic_path.h>
#include <ctrlpp/trajectory/trajectory.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    Eigen::Vector2d q0{0.0, 0.0}, q1{5.0, 3.0};
    double T = 4.0;
    auto const traj = ctrlpp::make_trajectory(ctrlpp::septic_path<double>, q0, q1, T);
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

- [septic-trajectory](septic-trajectory.md)<br/> polynomial with arbitrary BCs up to jerk
- [quintic-path](quintic-path.md)<br/> lower continuity (C2)
- [time-scaling](time-scaling.md)<br/> uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
