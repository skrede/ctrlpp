# bspline_trajectory

Scalar B-spline trajectory with compile-time degree and de Boor evaluation. Supports configurable degree (cubic = 3, quintic = 5, etc.), auto-generated uniform clamped knot vectors, and user-provided custom knot vectors. The factory function `make_bspline_interpolation()` solves for control points that pass through specified waypoints.

## Header

```cpp
#include "ctrlpp/trajectory/bspline_trajectory.h"
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type (`float`, `double`, `long double`) |
| `Degree` | `int`, at least 1 (a `static_assert`) | B-spline degree (e.g. 3 for cubic, 5 for quintic) |

## Config

```cpp
struct config {
    std::vector<Scalar> control_points;  // n+1 control points
    std::vector<Scalar> knot_vector{};   // If empty, auto-generate uniform clamped
};
```

If `knot_vector` is left empty, a uniform clamped knot vector is generated automatically with `m + 1 = n + Degree + 2` total knots over `[0, 1]`: `Degree + 1` zeros, equally spaced interior knots, `Degree + 1` ones. The parameter then runs over `[0, 1]` and `duration()` is 1, whatever the number of control points; a time scale comes from a knot vector the caller supplies or from `make_bspline_interpolation`. The control points are not passed through, apart from the first and last under a clamped knot vector.

## Factory

```cpp
static ctrlpp::expected<bspline_trajectory, spline_error> create(config const& cfg);
```

`create` is the only construction path. There is no non-fallible constructor: a rejected configuration is a value the caller has to inspect, never an object that quietly stands in for one. It validates the configuration and constructs a B-spline trajectory from control points and an optional knot vector. Requires at least `Degree + 1` control points. An empty knot vector skips the knot checks; a uniform clamped knot vector is generated instead, which is valid by construction.

Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| Fewer than `Degree + 1` control points | `spline_error::too_few_control_points` |
| A NaN or infinite control point or knot | `spline_error::non_finite_input` |
| Knot vector size differs from `control_points.size() + Degree + 1` | `spline_error::bad_knot_count` |
| Knot vector not non-decreasing | `spline_error::non_monotonic_knots` |
| Knot range `knot_vector.back() - knot_vector.front()` overflows `Scalar` | `spline_error::unrepresentable_spline` |
| Active domain `knot_vector[Degree]` .. `knot_vector[n + 1]` of zero length | `spline_error::non_increasing_times` |
| A derivative control point leaves the representable range of `Scalar` | `spline_error::unrepresentable_spline` |

The last check reads the derivative control points that `evaluate` uses for velocity and acceleration. Finite control points far enough apart that their difference quotient overflows are refused at construction rather than evaluating to NaN.

## Methods

### evaluate

```cpp
trajectory_point<Scalar, 1> evaluate(Scalar t) const;
```

Evaluate position, velocity, and acceleration at parameter `t` using de Boor's algorithm. The parameter is clamped to the active range [U[p], U[n+1]]. Velocity and acceleration are computed from derivative control points formed once in `create`, so `evaluate` does not allocate. For `Degree = 1` the acceleration is zero. Repeated knots are handled by the 0/0 := 0 convention of the B-spline recurrences, and any strictly positive knot span is evaluated as given, however short.

### duration

```cpp
Scalar duration() const;
```

Returns the active parameter range: U[n+1] - U[p].

## Factory Function

### make_bspline_interpolation

```cpp
template <typename Scalar, int Degree>
ctrlpp::expected<bspline_trajectory<Scalar, Degree>, spline_error> make_bspline_interpolation(
    std::vector<Scalar> const& times,
    std::vector<Scalar> const& positions);
```

Constructs a B-spline that passes through all waypoints at the given parameter values. Generates a clamped knot vector over `[times.front(), times.back()]` using de Boor's averaging method and solves the collocation system N * P = Q for the control points with a column-pivoting Householder QR. `Degree` is not deducible from the arguments, so both template arguments are written: `make_bspline_interpolation<double, 3>(times, positions)`.

Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| `times` and `positions` differ in length | `spline_error::size_mismatch` |
| Fewer than `Degree + 1` waypoints | `spline_error::too_few_points` |
| A NaN or infinite time or position | `spline_error::non_finite_input` |
| `times` not strictly increasing | `spline_error::non_increasing_times` |
| A span `times[i + 1] - times[i]` overflows `Scalar` | `spline_error::unrepresentable_spline` |
| A solved control point leaves the representable range of `Scalar` | `spline_error::unrepresentable_spline` |

The solved control points are checked directly because finite waypoints can still drive the solve out of range: a single waypoint of `1e308` among unit-scale neighbors yields infinite control points.

Any downstream `bspline_trajectory::create` failure is propagated.

## Free Function

### basis_function

```cpp
template <typename Scalar>
Scalar basis_function(std::int32_t i, std::int32_t p, Scalar t, std::vector<Scalar> const& U);
```

Evaluate the B-spline basis function B_{i,p}(t) using Cox-de Boor recursion over the knot vector `U`. It is a public entry to the recursion the interpolation factory uses to build its collocation matrix. Nothing is bounds-checked: `U` must hold at least `i + p + 2` knots.

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
