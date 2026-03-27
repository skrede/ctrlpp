# harmonic_path

Harmonic motion law: `q(tau) = (1 - cos(pi*tau)) / 2`. Produces a sinusoidal velocity profile. Non-zero acceleration at endpoints (not C2 at boundaries).

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/harmonic_path.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

## Functions

### harmonic_path

```cpp
template <typename Scalar>
auto harmonic_path(Scalar tau) -> path_point<Scalar>;
```

Evaluate the harmonic motion law at normalized time `tau` in [0,1].

### harmonic_path_peak_derivatives

```cpp
template <typename Scalar>
auto harmonic_path_peak_derivatives() -> std::array<Scalar, 3>;
```

Returns `{dq_max, ddq_max, dddq_max}` = `{pi/2, pi^2/2, pi^3/2}`.

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
#include "ctrlpp/traj/harmonic_path.h"
#include "ctrlpp/traj/trajectory.h"

auto pp = ctrlpp::harmonic_path(0.5);  // pp.q == 0.5, pp.dq == pi/2

Eigen::Vector2d q0{0, 0}, q1{1, 2};
auto traj = ctrlpp::make_trajectory(ctrlpp::harmonic_path<double>, q0, q1, 3.0);
```

## See Also

- [cycloidal-path](cycloidal-path.md) -- zero acceleration at endpoints (C2)
- [cubic-path](cubic-path.md) -- polynomial alternative (C1)
- [time-scaling](time-scaling.md) -- uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../reference/trajectory-generation.md)
