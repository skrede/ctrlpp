# cubic_path

Normalized cubic polynomial path: `q(tau) = 3*tau^2 - 2*tau^3`. C1-continuous with zero velocity at endpoints but non-zero acceleration.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/cubic_path.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

## Functions

### cubic_path

```cpp
template <typename Scalar>
auto cubic_path(Scalar tau) -> path_point<Scalar>;
```

Evaluate the cubic motion law at normalized time `tau` in [0,1]. Returns position, velocity, acceleration, and jerk.

### cubic_path_peak_derivatives

```cpp
template <typename Scalar>
auto cubic_path_peak_derivatives() -> std::array<Scalar, 3>;
```

Returns `{dq_max, ddq_max, dddq_max}` = `{1.5, 6.0, 12.0}`. Used by `compute_min_duration` for time scaling.

## Peak Derivatives

| Derivative | Value | Description |
|------------|-------|-------------|
| dq_max | 1.5 | Peak normalized velocity |
| ddq_max | 6.0 | Peak normalized acceleration |
| dddq_max | 12.0 | Peak normalized jerk (constant) |

## Usage Example

```cpp
#include "ctrlpp/traj/cubic_path.h"
#include "ctrlpp/traj/trajectory.h"

// As standalone path evaluation
auto pp = ctrlpp::cubic_path(0.5);  // pp.q == 0.5, pp.dq == 1.5

// Lifted to physical trajectory via adapter
Eigen::Vector2d q0{0, 0}, q1{1, 2};
auto traj = ctrlpp::make_trajectory(ctrlpp::cubic_path<double>, q0, q1, 2.0);
```

## See Also

- [cubic-trajectory](cubic-trajectory.md) -- polynomial with arbitrary velocity BCs
- [quintic-path](quintic-path.md) -- higher continuity (C2)
- [time-scaling](time-scaling.md) -- uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
