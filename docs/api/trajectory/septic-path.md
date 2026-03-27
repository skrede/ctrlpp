# septic_path

Normalized septic polynomial path: `q(tau) = 35*tau^4 - 84*tau^5 + 70*tau^6 - 20*tau^7`. C3-continuous with zero velocity, acceleration, and jerk at endpoints.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/septic_path.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

## Functions

### septic_path

```cpp
template <typename Scalar>
auto septic_path(Scalar tau) -> path_point<Scalar>;
```

Evaluate the septic motion law at normalized time `tau` in [0,1].

### septic_path_peak_derivatives

```cpp
template <typename Scalar>
auto septic_path_peak_derivatives() -> std::array<Scalar, 3>;
```

Returns `{dq_max, ddq_max, dddq_max}` = `{35/16, ~7.5132, 52.5}`. The ddq_max value involves nested radicals without a simple closed form.

## Peak Derivatives

| Derivative | Value | Description |
|------------|-------|-------------|
| dq_max | 2.1875 | Peak normalized velocity (35/16) |
| ddq_max | 7.5132 | Peak normalized acceleration (approximate) |
| dddq_max | 52.5 | Peak normalized jerk |

## Usage Example

```cpp
#include "ctrlpp/traj/septic_path.h"
#include "ctrlpp/traj/trajectory.h"

auto pp = ctrlpp::septic_path(0.5);  // pp.q == 0.5

Eigen::Vector2d q0{0, 0}, q1{5, 3};
auto traj = ctrlpp::make_trajectory(ctrlpp::septic_path<double>, q0, q1, 4.0);
```

## See Also

- [septic-trajectory](septic-trajectory.md) -- polynomial with arbitrary BCs up to jerk
- [quintic-path](quintic-path.md) -- lower continuity (C2)
- [time-scaling](time-scaling.md) -- uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
