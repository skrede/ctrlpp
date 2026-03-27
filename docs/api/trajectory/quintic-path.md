# quintic_path

Normalized quintic polynomial path: `q(tau) = 10*tau^3 - 15*tau^4 + 6*tau^5`. C2-continuous with zero velocity and acceleration at endpoints.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/quintic_path.h` |

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
#include "ctrlpp/traj/quintic_path.h"
#include "ctrlpp/traj/trajectory.h"

auto pp = ctrlpp::quintic_path(0.5);  // pp.q == 0.5, pp.dq == 1.875

Eigen::Vector3d q0{0, 0, 0}, q1{1, 1, 1};
auto traj = ctrlpp::make_trajectory(ctrlpp::quintic_path<double>, q0, q1, 3.0);
```

## See Also

- [quintic-trajectory](quintic-trajectory.md) -- polynomial with arbitrary velocity+acceleration BCs
- [cubic-path](cubic-path.md) -- lower continuity (C1)
- [septic-path](septic-path.md) -- higher continuity (C3)
- [time-scaling](time-scaling.md) -- uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
