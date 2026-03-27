# trajectory_types

Core output types for trajectory generation. `trajectory_point` holds ND-dimensional position, velocity, and acceleration vectors. `path_point` holds scalar normalized values for paths evaluated over [0,1].

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/trajectory_types.h` |

## Type: `trajectory_point<Scalar, ND>`

Output of trajectory segment evaluation.

| Field | Type | Description |
|-------|------|-------------|
| `position` | `Vector<Scalar, ND>` | Position vector |
| `velocity` | `Vector<Scalar, ND>` | Velocity vector |
| `acceleration` | `Vector<Scalar, ND>` | Acceleration vector |

## Type: `path_point<Scalar>`

Output of normalized path evaluation. Fields represent derivatives of the normalized position with respect to normalized time tau in [0,1].

| Field | Type | Description |
|-------|------|-------------|
| `q` | `Scalar` | Normalized position [0,1] |
| `dq` | `Scalar` | First derivative dq/dtau |
| `ddq` | `Scalar` | Second derivative d2q/dtau2 |
| `dddq` | `Scalar` | Third derivative d3q/dtau3 |

Physical values are obtained via kinematic scaling: `vel = h/T * dq`, `acc = h/T^2 * ddq`, `jerk = h/T^3 * dddq`.

## Usage Example

```cpp
#include "ctrlpp/traj/trajectory_types.h"

ctrlpp::trajectory_point<double, 3> pt;
pt.position = Eigen::Vector3d{1.0, 2.0, 3.0};
pt.velocity = Eigen::Vector3d::Zero();
pt.acceleration = Eigen::Vector3d::Zero();

ctrlpp::path_point<double> pp{.q = 0.5, .dq = 1.0, .ddq = 0.0, .dddq = 0.0};
```

## See Also

- [trajectory-segment](trajectory-segment.md) -- concept using `trajectory_point`
- [path-segment](path-segment.md) -- concept using `path_point`
- [Trajectory Generation Theory](../../reference/trajectory-generation.md)
