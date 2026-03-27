# synchronize

Multi-axis trajectory synchronization. Rescales multiple independent trajectory profiles to finish simultaneously at the duration of the slowest axis. Works with any profile satisfying the `syncable_profile` concept.

## Header

```cpp
#include "ctrlpp/traj/synchronize.h"
```

## Concept

### syncable_profile

```cpp
template <typename T>
concept syncable_profile = requires(T& p, typename T::scalar_type dur) {
    { p.duration() } -> std::convertible_to<typename T::scalar_type>;
    { p.rescale_to(dur) };
};
```

A profile is syncable if it exposes a `scalar_type` alias, a `duration()` method returning the current profile duration, and a `rescale_to(T)` method that stretches the profile to a new target duration.

The following trajectory types satisfy `syncable_profile`:

- `trapezoidal_trajectory`
- `double_s_trajectory`

## Free Functions

### synchronize (variadic)

```cpp
template <syncable_profile... Profiles>
void synchronize(Profiles&... profiles);
```

Synchronize a heterogeneous set of axis profiles. Finds the maximum duration across all profiles and calls `rescale_to(max_duration)` on each. Each profile's `rescale_to()` is responsible for maintaining constraint satisfaction (v_max, a_max, j_max).

### synchronize (vector)

```cpp
template <syncable_profile Profile>
void synchronize(std::vector<Profile>& profiles);
```

Runtime-sized variant for collections of identical profile types. Same behaviour as the variadic version.

## Multi-Axis Coordination

When multiple axes must move together (e.g. a 3-axis Cartesian robot), each axis computes its own time-optimal profile independently. `synchronize()` then stretches the faster axes to match the slowest, ensuring all axes start and stop together.

The time scaling preserves the profile shape -- the trajectory maintains its velocity/acceleration constraint structure, just stretched in time. Each profile's `rescale_to()` implementation handles the constraint-respecting time scaling internally.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'x', '' using 1:3 with lines title 'y', '' using 1:4 with lines title 'z'"

#include <ctrlpp/traj/synchronize.h>
#include <ctrlpp/traj/trapezoidal_trajectory.h>

#include <iostream>

int main()
{
    // Three-axis motion: each axis has different displacement
    ctrlpp::trapezoidal_trajectory<double> x_axis({.q0 = 0, .q1 = 10, .v_max = 2, .a_max = 5});
    ctrlpp::trapezoidal_trajectory<double> y_axis({.q0 = 0, .q1 = 3,  .v_max = 2, .a_max = 5});
    ctrlpp::trapezoidal_trajectory<double> z_axis({.q0 = 0, .q1 = 7,  .v_max = 2, .a_max = 5});

    // Before: each axis has different duration
    // After: all axes finish at the same time
    ctrlpp::synchronize(x_axis, y_axis, z_axis);

    double T = x_axis.duration();
    constexpr double dt = 0.01;
    for (double t = 0.0; t <= T; t += dt) {
        auto px = x_axis.evaluate(t);
        auto py = y_axis.evaluate(t);
        auto pz = z_axis.evaluate(t);
        std::cout << t << "," << px.position(0) << "," << py.position(0)
                  << "," << pz.position(0) << "\n";
    }
}
```

## See Also

- [trapezoidal-trajectory](trapezoidal-trajectory.md) -- Trapezoidal profile with `rescale_to()` support
- [double-s-trajectory](double-s-trajectory.md) -- Double-S profile with `rescale_to()` support
- [Trajectory Generation Theory](../../background/trajectory-generation.md) -- Multi-axis synchronization algorithms
