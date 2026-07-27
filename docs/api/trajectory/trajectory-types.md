# trajectory_types

Core output and error types for trajectory generation. `trajectory_point` holds ND-dimensional position, velocity, and acceleration vectors. `path_point` holds scalar normalized values for paths evaluated over [0,1]. `spline_error` and `trajectory_error` enumerate the structured failure modes of the spline and trajectory factories; each forms the error channel of a `ctrlpp::expected<T, E>` contract.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/trajectory_types.h` |

## Enum: `spline_error`

Structured failure modes for the spline factories: `cubic_spline`, `smoothing_spline`, `bspline_trajectory`, and `make_bspline_interpolation`. Returned through `ctrlpp::expected<T, spline_error>` from each `try_create`.

| Enumerator | Meaning |
|------------|---------|
| `too_few_points` | Fewer waypoints than the factory minimum (2 for cubic and smoothing splines, `Degree + 1` for B-spline interpolation) |
| `size_mismatch` | `times` and `positions` differ in length |
| `non_increasing_times` | Knot times are not strictly increasing |
| `periodic_endpoint_mismatch` | Periodic boundary conditions require the first and last positions to match within the endpoint rounding budget |
| `periodic_too_few_points` | Periodic boundary conditions require at least 3 waypoints |
| `too_few_control_points` | A B-spline of degree p requires at least p + 1 control points |
| `bad_knot_count` | Knot vector size differs from `control_points.size() + Degree + 1` |
| `non_monotonic_knots` | Knot vector is not non-decreasing |
| `mu_out_of_range` | Smoothing parameter mu lies outside (0, 1] |

## Enum: `trajectory_error`

Structured failure modes for the point-to-point trajectory factories and the online trajectory planners.

| Enumerator | Meaning |
|------------|---------|
| `non_positive_velocity_limit` | The velocity limit must be positive |
| `non_positive_acceleration_limit` | The acceleration limit must be positive |
| `non_positive_jerk_limit` | The jerk limit must be positive |
| `non_positive_duration` | The requested duration must be positive |
| `non_finite_input` | A boundary value or limit is NaN/Inf |
| `boundary_velocity_exceeds_limit` | A commanded boundary velocity is larger in magnitude than the velocity limit it is commanded under. The limit is a precondition of the point-to-point profiles, not a value they raise to fit: raising it would violate a bound the caller asked for, and honoring it would require a ramp that runs backwards in time |
| `unrepresentable_duration` | A duration of the constructed profile is not representable in the scalar type, either because it left the finite range or because the total underflowed to zero on a command with a nonzero displacement, which would report an instantaneous traversal |
| `unreachable_boundary_velocity` | The commanded displacement is shorter than the distance the fastest admissible transition between the two boundary velocities already sweeps, so no profile of the requested shape realizes it |
| `duration_shorter_than_current` | Time rescaling only slows a profile down. The profile already runs at the fastest shape its limits allow, so a duration below the current one is not realizable |
| `unreachable_duration` | The requested duration lies outside the set the commanded displacement, the kinematic limits, and the boundary velocities can realize together |

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

#include "ctrlpp/trajectory/cubic_trajectory.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec1 = Eigen::Matrix<double, 1, 1>;
    auto traj = ctrlpp::make_cubic_trajectory(Vec1{0.0}, Vec1{1.0}, Vec1{0.0}, Vec1{0.0}, 2.0).value();
    for (double t = 0; t <= 2.0; t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [trajectory-segment](trajectory-segment.md)<br/> concept using `trajectory_point`
- [path-segment](path-segment.md)<br/> concept using `path_point`
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
