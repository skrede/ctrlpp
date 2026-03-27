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

## Constructor

```cpp
explicit online_planner_3rd(config const& cfg);
```

Constructs a planner with kinematic limits. Initial state is at rest at q = 0 with zero acceleration.

## Methods

### update

```cpp
void update(Scalar target);
```

Set a new target position and replan from the current state. Computes a time-optimal double-S profile from (q, v, a) to (target, 0, 0) respecting `v_max`, `a_max`, and `j_max`. Multi-phase deceleration is used when current velocity or acceleration require braking before replanning.

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

## Profile Behaviour

The planner produces double-S velocity profiles composed of constant-jerk phases (up to 11 phases total):

1. **Acceleration ramp-up** -- apply +j_max until a_max reached
2. **Constant acceleration** -- hold at a_max (may be zero duration)
3. **Acceleration ramp-down** -- apply -j_max to bring acceleration to zero
4. **Cruise** -- hold at v_max with zero acceleration (may be zero duration)
5. **Deceleration ramp-up** -- apply -j_max to build deceleration
6. **Constant deceleration** -- hold at -a_max (may be zero duration)
7. **Deceleration ramp-down** -- apply +j_max to bring acceleration and velocity to zero

Degenerate cases (v_max or a_max not reached) automatically reduce the number of phases. Mid-motion replanning uses a brake-to-zero-then-replan strategy for robustness.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include <ctrlpp/trajectory/online_planner_3rd.h>

#include <iostream>

int main()
{
    ctrlpp::online_planner_3rd<double> planner({
        .v_max = 1.0,
        .a_max = 5.0,
        .j_max = 50.0,
    });
    planner.update(10.0);  // move to position 10

    constexpr double dt = 0.001;
    for (double t = 0.0; !planner.is_settled(); t += dt) {
        auto pt = planner.sample(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [online-planner-2nd](online-planner-2nd.md) -- 2nd-order variant (faster, acceleration-limited only)
- [double-s-trajectory](double-s-trajectory.md) -- Pre-computed double-S profile
- [Trajectory Generation Theory](../../background/trajectory-generation.md) -- Double-S profile mathematics and jerk limitation
