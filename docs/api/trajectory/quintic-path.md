# quintic_path

Normalized quintic polynomial path: `q(tau) = 10*tau^3 - 15*tau^4 + 6*tau^5`. C2-continuous with zero velocity and acceleration at endpoints.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/quintic_path.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

## Functions

### quintic_path

```cpp
template <typename Scalar>
auto quintic_path(Scalar tau) -> path_point<Scalar>;
```

Evaluate the quintic motion law at normalized time `tau` in [0,1].

### quintic_path_peak_derivatives

```cpp
template <typename Scalar>
auto quintic_path_peak_derivatives() -> std::array<Scalar, 3>;
```

Returns `{dq_max, ddq_max, dddq_max}` = `{15/8, 10*sqrt(3)/3, 60}`. Not constexpr due to `sqrt`.

## Peak Derivatives

| Derivative | Value | Description |
|------------|-------|-------------|
| dq_max | 1.875 | Peak normalized velocity |
| ddq_max | 5.7735 | Peak normalized acceleration (10*sqrt(3)/3) |
| dddq_max | 60.0 | Peak normalized jerk |

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include <ctrlpp/trajectory/quintic_path.h>
#include <ctrlpp/trajectory/trajectory.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    Eigen::Vector2d q0{0.0, 0.0}, q1{2.0, 1.0};
    double T = 3.0;
    auto traj = ctrlpp::make_trajectory(ctrlpp::quintic_path<double>, q0, q1, T);

    constexpr double dt = 0.01;
    for (double t = 0.0; t <= T; t += dt) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0)
                  << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [quintic-trajectory](quintic-trajectory.md) -- polynomial with arbitrary velocity+acceleration BCs
- [cubic-path](cubic-path.md) -- lower continuity (C1)
- [septic-path](septic-path.md) -- higher continuity (C3)
- [time-scaling](time-scaling.md) -- uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
