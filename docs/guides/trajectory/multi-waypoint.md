# Multi-Waypoint Paths

Interpolating sequences of waypoints with splines and B-splines.

## When to Use What

| Spline Type | Endpoints | Best For |
|-------------|----------|----------|
| [Natural cubic spline](../../api/trajectory/cubic-spline.md) | Zero-acceleration BCs | General-purpose waypoint interpolation |
| [Clamped cubic spline](../../api/trajectory/cubic-spline.md) | Specified endpoint velocities | Known start/end velocity requirements |
| [Periodic cubic spline](../../api/trajectory/cubic-spline.md) | Cyclic (wraps around) | Repeating motions and closed paths |
| [Smoothing spline](../../api/trajectory/smoothing-spline.md) | Configurable smoothing | Noisy waypoints needing approximation |
| [B-spline](../../api/trajectory/bspline-trajectory.md) | Configurable degree | Local control, derivative constraints |

## Quick Start: Cubic Spline

A natural cubic spline interpolates through all waypoints with C2
continuity. Provide time stamps and positions:

```cpp
#include <ctrlpp/trajectory/cubic_spline.h>

#include <Eigen/Core>
#include <iostream>
#include <vector>

std::vector<double> times = {0.0, 1.0, 2.5, 4.0};

std::vector<Eigen::Vector2d> positions = {
    {0.0, 0.0}, {1.0, 0.5}, {2.0, 1.5}, {3.0, 0.0}};

auto spline = ctrlpp::cubic_spline(
    times, positions, ctrlpp::boundary::natural);

auto [pos, vel, acc] = spline.evaluate(1.5);
std::cout << "position at t=1.5: " << pos.transpose() << "\n";
```

See [example 06](../../../examples/trajectory/ctrlpp_trajectory_06_cubic_spline.cpp)
for a runnable version.

## Quick Start: B-spline

B-splines offer local control -- moving one control point only affects
nearby segments. The degree parameter controls smoothness:

```cpp
#include <ctrlpp/trajectory/bspline_trajectory.h>

#include <Eigen/Core>
#include <iostream>
#include <vector>

std::vector<Eigen::Vector3d> control_points = {
    {0.0, 0.0, 0.0}, {1.0, 2.0, 0.0},
    {3.0, 2.0, 1.0}, {4.0, 0.0, 0.0}};

auto bspline = ctrlpp::bspline_trajectory<double, 3>(
    control_points, 3);  // degree 3

auto [pos, vel, acc] = bspline.evaluate(0.5);
std::cout << "position at t=0.5: " << pos.transpose() << "\n";
```

See [example 07](../../../examples/trajectory/ctrlpp_trajectory_07_bspline.cpp)
for a runnable version.

## Smoothing Noisy Data

The [smoothing spline](../../api/trajectory/smoothing-spline.md) balances
interpolation fidelity against curve smoothness via a single parameter `mu`:

- `mu = 0` -- exact interpolation (passes through every waypoint)
- `mu = 1` -- maximum smoothing (least-squares straight line)

For noisy sensor data, values around `mu = 0.5` to `mu = 0.9` typically
produce good results. The optimal value depends on the noise level relative
to the signal.

## Piecewise Composition

When different path segments need different trajectory types (e.g., a
linear ramp followed by a spline followed by a polynomial), use
[piecewise_trajectory](../../api/trajectory/piecewise-trajectory.md) to stitch
them into a single evaluable trajectory with automatic time management.

## Links

**API reference:**
[cubic spline](../../api/trajectory/cubic-spline.md) |
[smoothing spline](../../api/trajectory/smoothing-spline.md) |
[B-spline](../../api/trajectory/bspline-trajectory.md) |
[piecewise trajectory](../../api/trajectory/piecewise-trajectory.md)

**Examples:**
[06 cubic spline](../../../examples/trajectory/ctrlpp_trajectory_06_cubic_spline.cpp) |
[07 B-spline](../../../examples/trajectory/ctrlpp_trajectory_07_bspline.cpp)
