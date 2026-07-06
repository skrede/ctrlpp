# cubic_spline

Cubic spline interpolation through waypoints with C2 continuity. Supports natural (zero endpoint acceleration), clamped (assigned endpoint velocities), and periodic (cyclic) boundary conditions. Internally solves a tridiagonal system for spline velocities using the velocity-based formulation.

## Header

```cpp
#include "ctrlpp/trajectory/cubic_spline.h"
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
| `periodic` | Cyclic: v_0 = v_n and a_0 = a_n (requires q_0 = q_n and at least 3 waypoints) |

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

## Factory

```cpp
[[nodiscard]] static auto try_create(config const& cfg)
    -> ctrlpp::expected<cubic_spline, spline_error>;
```

Validates the configuration and constructs a cubic spline. Requires at least 2 waypoints and strictly increasing time values. Periodic boundary conditions additionally require at least 3 waypoints and matching first/last positions: with only 2 waypoints the cyclic system for the interior velocities is empty and the closed curve degenerates, so that configuration is rejected.

Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| Fewer than 2 waypoints | `spline_error::too_few_points` |
| `times` and `positions` differ in length | `spline_error::size_mismatch` |
| Knot times not strictly increasing | `spline_error::non_increasing_times` |
| Periodic BC with fewer than 3 waypoints | `spline_error::periodic_too_few_points` |
| Periodic BC with q_0 != q_n beyond the rounding budget | `spline_error::periodic_endpoint_mismatch` |

## Constructor

```cpp
explicit cubic_spline(config const& cfg);  // requires CTRLPP_HAS_EXCEPTIONS
```

Throwing convenience wrapper over `try_create`: delegates to `try_create(cfg).value()`, so an invalid configuration throws the `value()` exception of `ctrlpp::expected`. Compiled out when `CTRLPP_HAS_EXCEPTIONS` is 0.

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
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include <ctrlpp/trajectory/cubic_spline.h>

#include <iostream>

int main()
{
    auto spline = ctrlpp::cubic_spline<double>::try_create({
        .times = {0.0, 1.0, 2.0, 3.0, 4.0},
        .positions = {0.0, 1.0, 0.5, 1.5, 2.0},
        .bc = ctrlpp::boundary_condition::natural,
    });
    if (!spline.has_value()) {
        return 1;
    }

    double T = spline->duration();
    constexpr double dt = 0.01;
    for (double t = 0.0; t <= T; t += dt) {
        auto pt = spline->evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [smoothing-spline](smoothing-spline.md)<br/> Smoothing spline approximation with data/smoothness tradeoff
- [bspline-trajectory](bspline-trajectory.md)<br/> B-spline trajectory with configurable degree
- [trajectory-types](trajectory-types.md)<br/> `spline_error` enumerators returned by `try_create`
- [Trajectory Generation Theory](../../background/trajectory-generation.md)<br/> Mathematical background for spline interpolation
