# cubic_spline

Cubic spline interpolation through waypoints with C2 continuity. Supports natural (zero endpoint acceleration), clamped (assigned endpoint velocities), and periodic (cyclic) boundary conditions. Internally solves a tridiagonal system for spline velocities using the velocity-based formulation.

## Header

```cpp
#include "ctrlpp/traj/cubic_spline.h"
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (`float`, `double`, `long double`) |

## Boundary Conditions

```cpp
enum class boundary_condition { natural, clamped, periodic };
```

| Value | Meaning |
|-------|---------|
| `natural` | Zero acceleration at endpoints (M_0 = M_n = 0) |
| `clamped` | Endpoint velocities v_0 and v_n specified by user |
| `periodic` | Cyclic: v_0 = v_n and a_0 = a_n (requires q_0 = q_n) |

## Config

```cpp
struct config {
    std::vector<Scalar> times;       // Knot times t_0 ... t_n (n+1 entries)
    std::vector<Scalar> positions;   // Waypoint positions q_0 ... q_n (n+1 entries)
    boundary_condition bc{boundary_condition::natural};
    Scalar v0{};                     // Endpoint velocity for clamped BC
    Scalar vn{};                     // Endpoint velocity for clamped BC
};
```

## Constructor

```cpp
explicit cubic_spline(config const& cfg);
```

Constructs a cubic spline from waypoints and boundary conditions. Requires at least 2 waypoints and strictly increasing time values.

## Methods

### evaluate

```cpp
auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>;
```

Evaluate position, velocity, and acceleration at time `t`. Time is clamped to [t_0, t_n]. Uses binary search to find the active span followed by Horner polynomial evaluation.

### duration

```cpp
auto duration() const -> Scalar;
```

Returns total spline duration: t_n - t_0.

## Concept Satisfaction

`cubic_spline<Scalar>` satisfies `trajectory_segment<cubic_spline<Scalar>, Scalar, 1>`.

## Usage Example

```cpp
#include "ctrlpp/traj/cubic_spline.h"

ctrlpp::cubic_spline<double> spline({
    .times = {0.0, 1.0, 2.0, 3.0, 4.0},
    .positions = {0.0, 1.0, 0.5, 1.5, 2.0},
    .bc = ctrlpp::boundary_condition::natural,
});

auto pt = spline.evaluate(1.5);
// pt.position, pt.velocity, pt.acceleration
```

## See Also

- [smoothing-spline](smoothing-spline.md) -- Smoothing spline approximation with data/smoothness tradeoff
- [bspline-trajectory](bspline-trajectory.md) -- B-spline trajectory with configurable degree
- [Trajectory Generation Theory](../../background/trajectory-generation.md) -- Mathematical background for spline interpolation
