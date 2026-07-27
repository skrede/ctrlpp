# synchronize

Multi-axis trajectory synchronization. Retimes multiple independent trajectory profiles to finish simultaneously at the duration of the slowest axis. Works with any profile satisfying the `syncable_profile` concept.

## Header

```cpp
#include "ctrlpp/trajectory/synchronize.h"
```

## Concept

### syncable_profile

```cpp
template <typename T>
concept syncable_profile = requires(T& p, T const& cp, typename T::scalar_type dur) {
    { p.duration() } -> std::convertible_to<typename T::scalar_type>;
    { p.rescale_to(dur) } -> std::same_as<ctrlpp::expected<void, trajectory_error>>;
    { cp.can_rescale_to(dur) } -> std::same_as<ctrlpp::expected<void, trajectory_error>>;
};
```

A profile is syncable if it exposes a `scalar_type` alias, a `duration()` method returning the current profile duration, a `rescale_to(T)` method that retimes the profile and reports whether the request is realizable, and a const `can_rescale_to(T)` query with the same return that answers the same question without mutating. The rescaling return is constrained exactly, so a profile whose rescaling cannot report failure fails at the concept rather than at the call site.

The following trajectory types satisfy `syncable_profile`:

- `trapezoidal_trajectory`
- `double_s_trajectory`

## Free Functions

### synchronize (variadic)

```cpp
template <syncable_profile... Profiles>
[[nodiscard]] auto synchronize(Profiles&... profiles)
    -> ctrlpp::expected<void, trajectory_error>;
```

Synchronize a heterogeneous set of axis profiles. Finds the maximum duration across all profiles and retimes each to it.

### synchronize (contiguous view)

```cpp
template <syncable_profile Profile>
[[nodiscard]] auto synchronize(std::span<Profile> profiles)
    -> ctrlpp::expected<void, trajectory_error>;
```

Runtime-sized variant for a contiguous run of identical profile types. The view carries no ownership, so the same call serves a `std::vector`, a plain array, or a statically allocated buffer, and nothing on this path allocates. Construct the view at the call site:

```cpp
std::vector<ctrlpp::trapezoidal_trajectory<double>> axes = /* ... */;
auto const synced = ctrlpp::synchronize(std::span{axes});
```

## All or nothing

Both overloads run in two passes. Every axis is first asked, through `can_rescale_to()`, whether it can reach the target duration; only once all of them have answered yes is any axis retimed. A rejection on any axis returns that axis's `trajectory_error` with nothing mutated, so a partially synchronized set is never handed back.

The two passes go through the same solve inside each profile, so the checking pass cannot pass an axis that the committing pass then fails. Neither pass builds an owning copy of anything, which is why the runtime-sized overload takes a view rather than a container.

The slowest axis never trips the shortening rejection. The target is a bit-exact copy of that axis's own reported duration, so it compares equal and takes the success no-op path; the value compared is the axis's own stored duration rather than a recomputed one, so there is no float-equality fragility in it.

## Multi-Axis Coordination

When multiple axes must move together (for example a 3-axis Cartesian robot), each axis computes its own time-optimal profile independently. `synchronize()` then slows the faster axes to match the slowest, so all axes start and stop together.

Each profile is rebuilt at the longer duration rather than patched, so the commanded displacement and the commanded boundary velocities still hold afterwards. Retiming is only ever a slowdown, and an axis whose displacement, limits, and boundary velocities cannot reach the target duration is reported instead of being retimed to something else. See the rescaling entries on both profile pages for the exact rejection lists.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'x', '' using 1:3 with lines title 'y', '' using 1:4 with lines title 'z'"

#include <ctrlpp/trajectory/synchronize.h>
#include <ctrlpp/trajectory/trapezoidal_trajectory.h>

#include <iostream>

int main()
{
    // Three-axis motion: each axis has different displacement
    ctrlpp::trapezoidal_trajectory<double> x_axis({.q0 = 0, .q1 = 10, .v_max = 2, .a_max = 5});
    ctrlpp::trapezoidal_trajectory<double> y_axis({.q0 = 0, .q1 = 3,  .v_max = 2, .a_max = 5});
    ctrlpp::trapezoidal_trajectory<double> z_axis({.q0 = 0, .q1 = 7,  .v_max = 2, .a_max = 5});

    // Before: each axis has different duration
    // After: all axes finish at the same time
    auto const synced = ctrlpp::synchronize(x_axis, y_axis, z_axis);
    if (!synced) {
        std::cerr << "Synchronization rejected: the slowest axis duration is not "
                     "reachable for every axis\n";
        return 1;
    }

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

- [trapezoidal-trajectory](trapezoidal-trajectory.md)<br/> Trapezoidal profile with `rescale_to()` support
- [double-s-trajectory](double-s-trajectory.md)<br/> Double-S profile with `rescale_to()` support
- [trajectory-types](trajectory-types.md)<br/> `trajectory_error`, the rejection channel of both
- [Trajectory Generation Theory](../../background/trajectory-generation.md)<br/> Multi-axis synchronization algorithms
