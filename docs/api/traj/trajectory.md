# trajectory

Adapter that lifts a normalized path (evaluating in [0,1]) into a physical trajectory segment with position, velocity, and acceleration over a specified duration and displacement.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/trajectory.h` |
| **Factory** | `ctrlpp::make_trajectory` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Law` | Callable `(Scalar tau) -> path_point<Scalar>` |
| `Scalar` | Floating-point type |
| `ND` | Number of dimensions (`std::size_t`) |

## Kinematic Scaling

The adapter applies displacement and duration scaling to the normalized path output:

- **position:** `q0 + h * law(tau).q`
- **velocity:** `h/T * law(tau).dq`
- **acceleration:** `h/T^2 * law(tau).ddq`

where `h = q1 - q0`, `T = duration`, and `tau = clamp(t, 0, T) / T`.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, ND> evaluate(Scalar t) const` | Evaluate at time `t`, clamped to `[0, T]` |
| `duration` | `Scalar duration() const` | Total duration `T` |

## Factory Function

```cpp
template <typename Law, typename Scalar, int Rows>
auto make_trajectory(
    Law law,
    Eigen::Matrix<Scalar, Rows, 1> const& q0,
    Eigen::Matrix<Scalar, Rows, 1> const& q1,
    Scalar duration) -> trajectory<Law, Scalar, ND>;
```

Computes displacement `h = q1 - q0` internally.

## Usage Example

```cpp
#include "ctrlpp/traj/trajectory.h"
#include "ctrlpp/traj/cycloidal_path.h"

Eigen::Vector2d q0{0.0, 0.0}, q1{1.0, 2.0};
auto traj = ctrlpp::make_trajectory(ctrlpp::cycloidal_path<double>, q0, q1, 2.0);
auto pt = traj.evaluate(1.0);  // midpoint
// pt.position, pt.velocity, pt.acceleration are Eigen::Vector2d
```

## See Also

- [trajectory-segment](trajectory-segment.md) -- concept this type satisfies
- [path-segment](path-segment.md) -- concept for the `Law` parameter
- [Trajectory Generation Theory](../../reference/trajectory-generation.md)
