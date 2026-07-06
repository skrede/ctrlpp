# online_planner_3rd

3rd-order online trajectory planner that generates double-S (jerk-limited) velocity profiles in real time. On each `update(target)`, the planner computes a time-optimal profile from the current state to the target position, respecting `v_max`, `a_max`, and `j_max` constraints. Produces smoother motion than the 2nd-order planner at the cost of slightly longer move times.

Unlike pre-computed trajectory segments, online planners are stateful filters with no fixed duration and do not satisfy `trajectory_segment`.

## Header

```cpp
#include "ctrlpp/trajectory/online_planner_3rd.h"
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (`float`, `double`, `long double`) |

## Config

```cpp
struct config {
    Scalar v_max;   // Maximum velocity magnitude
    Scalar a_max;   // Maximum acceleration magnitude
    Scalar j_max;   // Maximum jerk magnitude
};
```

All three limits divide in the planner math (cruise duration `h / v_max`, jerk-phase durations `a_max / j_max` and `|a| / j_max`, the `a_max`-reached threshold `a_max^2 / j_max`), so the domain of each is finite and strictly positive.

## Construction

```cpp
[[nodiscard]] static auto try_create(config const& cfg)
    -> ctrlpp::expected<online_planner_3rd, trajectory_error>;
```

Validates the kinematic limits and constructs a planner with initial state at rest at q = 0 with zero acceleration. Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| NaN/Inf or non-positive `v_max` | `trajectory_error::non_positive_velocity_limit` |
| NaN/Inf or non-positive `a_max` | `trajectory_error::non_positive_acceleration_limit` |
| NaN/Inf or non-positive `j_max` | `trajectory_error::non_positive_jerk_limit` |

```cpp
explicit online_planner_3rd(config const& cfg);
```

Throwing convenience wrapper over `try_create`; delegates to `try_create(cfg).value()`. Only available when `CTRLPP_HAS_EXCEPTIONS` is 1; prefer `try_create` on exception-free builds.

## Methods

### update

```cpp
void update(Scalar target);
```

Set a new target position and replan from the current state. Computes a time-optimal double-S profile from (q, v, a) to (target, 0, 0) respecting `v_max`, `a_max`, and `j_max`. A same-direction move carries the current velocity through the profile (no full-stop dip); a velocity pointing away from the target, or too large to stop in the available distance, is braked to rest first and then replanned.

### sample

```cpp
auto sample(Scalar t) const -> trajectory_point<Scalar, 1>;
```

Evaluate the planned trajectory at absolute time `t`. Returns position, velocity, and acceleration. Updates internal state for future `update()` calls.

### is_settled

```cpp
auto is_settled() const -> bool;
```

Returns true when the planner has reached the target with zero velocity and zero acceleration.

### reset

```cpp
void reset(Scalar q0);
```

Reset state to position `q0` with zero velocity and zero acceleration.

## Profile Behavior

The planner produces double-S velocity profiles composed of constant-jerk phases (up to 11 phases total):

1. **Acceleration ramp-up**<br/>apply +j_max until a_max reached
2. **Constant acceleration**<br/>hold at a_max (may be zero duration)
3. **Acceleration ramp-down**<br/>apply -j_max to bring acceleration to zero
4. **Cruise**<br/>hold at v_max with zero acceleration (may be zero duration)
5. **Deceleration ramp-up**<br/>apply -j_max to build deceleration
6. **Constant deceleration**<br/>hold at -a_max (may be zero duration)
7. **Deceleration ramp-down**<br/>apply +j_max to bring acceleration and velocity to zero

Degenerate cases (v_max or a_max not reached) automatically reduce the number of phases. Mid-motion replanning uses a brake-to-zero-then-replan strategy for robustness.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include <ctrlpp/trajectory/online_planner_3rd.h>

#include <iostream>

int main()
{
    auto result = ctrlpp::online_planner_3rd<double>::try_create({
        .v_max = 1.0,
        .a_max = 5.0,
        .j_max = 50.0,
    });
    if (!result.has_value())
        return 1;
    auto& planner = *result;
    planner.update(10.0);  // move to position 10

    constexpr double dt = 0.001;
    for (double t = 0.0; !planner.is_settled(); t += dt) {
        auto pt = planner.sample(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [online-planner-2nd](online-planner-2nd.md)<br/> 2nd-order variant (faster, acceleration-limited only)
- [double-s-trajectory](double-s-trajectory.md)<br/> Pre-computed double-S profile
- [Trajectory Generation Theory](../../background/trajectory-generation.md)<br/> Double-S profile mathematics and jerk limitation
