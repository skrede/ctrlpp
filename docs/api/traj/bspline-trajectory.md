# bspline_trajectory

B-spline trajectory with compile-time degree and de Boor evaluation. Supports configurable degree (cubic = 3, quintic = 5, etc.), auto-generated uniform clamped knot vectors, and user-provided custom knot vectors. The factory function `make_bspline_interpolation()` solves for control points that pass through specified waypoints.

## Header

```cpp
#include "ctrlpp/traj/bspline_trajectory.h"
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (`float`, `double`, `long double`) |
| `Degree` | `int` | B-spline degree (e.g. 3 for cubic, 5 for quintic) |

## Config

```cpp
struct config {
    std::vector<Scalar> control_points;  // n+1 control points
    std::vector<Scalar> knot_vector{};   // If empty, auto-generate uniform clamped
};
```

If `knot_vector` is left empty, a uniform clamped knot vector is generated automatically with `m + 1 = n + Degree + 2` total knots.

## Constructor

```cpp
explicit bspline_trajectory(config const& cfg);
```

Constructs a B-spline trajectory from control points and an optional knot vector. Requires at least `Degree + 1` control points. Throws `std::invalid_argument` if the knot vector has wrong size or is non-monotonic.

## Methods

### evaluate

```cpp
auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>;
```

Evaluate position, velocity, and acceleration at parameter `t` using de Boor's algorithm. The parameter is clamped to the active range [U[p], U[n+1]]. Velocity and acceleration are computed from derivative control points.

### duration

```cpp
auto duration() const -> Scalar;
```

Returns the active parameter range: U[n+1] - U[p].

## Factory Function

### make_bspline_interpolation

```cpp
template <typename Scalar, int Degree>
auto make_bspline_interpolation(
    std::vector<Scalar> const& times,
    std::vector<Scalar> const& positions) -> bspline_trajectory<Scalar, Degree>;
```

Constructs a B-spline that passes through all waypoints at the given parameter values. Generates a clamped knot vector using de Boor's averaging method and solves the interpolation matrix N * P = Q for control points.

## Free Function

### basis_function

```cpp
template <typename Scalar>
auto basis_function(int i, int p, Scalar t, std::vector<Scalar> const& U) -> Scalar;
```

Evaluate the B-spline basis function B_{i,p}(t) using Cox-de Boor recursion. Primarily used internally by the interpolation factory.

## Concept Satisfaction

`bspline_trajectory<Scalar, Degree>` satisfies `trajectory_segment<bspline_trajectory<Scalar, Degree>, Scalar, 1>`.

## Usage Example

```cpp
#include "ctrlpp/traj/bspline_trajectory.h"

// Cubic B-spline through waypoints via interpolation factory
auto bspline = ctrlpp::make_bspline_interpolation<double, 3>(
    {0.0, 1.0, 2.0, 3.0, 4.0},  // parameter values
    {0.0, 1.0, 0.5, 1.5, 2.0}   // positions
);

auto pt = bspline.evaluate(1.5);
// pt.position, pt.velocity, pt.acceleration

// Or construct directly with known control points
ctrlpp::bspline_trajectory<double, 3> manual({
    .control_points = {0.0, 0.5, 1.0, 1.5, 2.0},
});
```

## See Also

- [cubic-spline](cubic-spline.md) -- Simpler cubic interpolation for moderate waypoint counts
- [smoothing-spline](smoothing-spline.md) -- Spline approximation with noise filtering
- [Trajectory Generation Theory](../../background/trajectory-generation.md) -- B-spline basis functions and de Boor's algorithm
