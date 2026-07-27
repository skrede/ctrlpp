# bspline_trajectory

B-spline trajectory with compile-time degree and de Boor evaluation. Supports configurable degree (cubic = 3, quintic = 5, etc.), auto-generated uniform clamped knot vectors, and user-provided custom knot vectors. The factory function `make_bspline_interpolation()` solves for control points that pass through specified waypoints.

## Header

```cpp
#include "ctrlpp/trajectory/bspline_trajectory.h"
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

## Factory

```cpp
static auto create(config const& cfg)
    -> ctrlpp::expected<bspline_trajectory, spline_error>;
```

`create` is the only construction path. There is no non-fallible constructor: a rejected configuration is a value the caller has to inspect, never an object that quietly stands in for one. It validates the configuration and constructs a B-spline trajectory from control points and an optional knot vector. Requires at least `Degree + 1` control points. An empty knot vector skips the knot checks; a uniform clamped knot vector is generated instead, which is valid by construction.

Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| Fewer than `Degree + 1` control points | `spline_error::too_few_control_points` |
| Knot vector size differs from `control_points.size() + Degree + 1` | `spline_error::bad_knot_count` |
| Knot vector not non-decreasing | `spline_error::non_monotonic_knots` |

## Methods

### evaluate

```cpp
auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>;
```

Evaluate position, velocity, and acceleration at parameter `t` using de Boor's algorithm. The parameter is clamped to the active range [U[p], U[n+1]]. Velocity and acceleration are computed from derivative control points. Near-zero knot spans are handled with a relative epsilon threshold scaled by knot magnitude, preventing overflow when adjacent knots are very close but not identical.

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
    std::vector<Scalar> const& positions)
    -> ctrlpp::expected<bspline_trajectory<Scalar, Degree>, spline_error>;
```

Constructs a B-spline that passes through all waypoints at the given parameter values. Generates a clamped knot vector using de Boor's averaging method and solves the interpolation matrix N * P = Q for control points.

Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| `times` and `positions` differ in length | `spline_error::size_mismatch` |
| Fewer than `Degree + 1` waypoints | `spline_error::too_few_points` |

Any downstream `bspline_trajectory::create` failure is propagated.

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
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include <ctrlpp/trajectory/bspline_trajectory.h>

#include <iostream>

int main()
{
    // Cubic B-spline through waypoints via interpolation factory
    auto bspline = ctrlpp::make_bspline_interpolation<double, 3>(
        {0.0, 1.0, 2.0, 3.0, 4.0},  // parameter values
        {0.0, 1.0, 0.5, 1.5, 2.0}   // positions
    );
    if (!bspline.has_value()) {
        return 1;
    }

    double T = bspline->duration();
    constexpr double dt = 0.01;
    for (double t = 0.0; t <= T; t += dt) {
        auto pt = bspline->evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [cubic-spline](cubic-spline.md)<br/> Simpler cubic interpolation for moderate waypoint counts
- [smoothing-spline](smoothing-spline.md)<br/> Spline approximation with noise filtering
- [trajectory-types](trajectory-types.md)<br/> `spline_error` enumerators returned by `create`
- [Trajectory Generation Theory](../../background/trajectory-generation.md)<br/> B-spline basis functions and de Boor's algorithm
