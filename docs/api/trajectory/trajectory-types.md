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
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include "ctrlpp/traj/cubic_trajectory.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    auto traj = ctrlpp::make_cubic_trajectory(Vec1{0.0}, Vec1{1.0}, Vec1{0.0}, Vec1{0.0}, 2.0);
    for (double t = 0; t <= 2.0; t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [trajectory-segment](trajectory-segment.md) -- concept using `trajectory_point`
- [path-segment](path-segment.md) -- concept using `path_point`
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
