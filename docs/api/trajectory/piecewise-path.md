# piecewise_path

Variadic composition of heterogeneous path segments into a single path. Segments are evaluated in sequence with automatic parameter offsets at the breakpoints. The composed path itself satisfies `path_segment`, enabling recursive composition.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/piecewise_path.h` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | none (`typename`) | The parameter and value type |
| `Segments...` | `path_segment<Segment, Scalar>` for each | The segment types, in evaluation order |

## Construction

```cpp
explicit piecewise_path(Segments... segs);
```

Segments are taken by value and moved into internal storage. Breakpoints are the cumulative sums of the segments' `duration()`, computed once here. Nothing is validated: each duration must be finite and nonnegative.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `path_point<Scalar> evaluate(Scalar t) const` | Evaluate at parameter `t`, dispatched to the owning segment |
| `duration` | `Scalar duration() const` | Sum of all segment durations |

## Static Members

| Member | Type | Description |
|--------|------|-------------|
| `segment_count` | `std::size_t` | Number of segments (compile-time) |

## Dispatch Logic

`t` is clamped to `[0, duration()]`, so the composed parameter runs over the sum of the segment durations rather than over `[0, 1]`. Segment `i` owns `[b_(i-1), b_i)`, where `b_i` is the sum of the first `i + 1` durations and `b_(-1) = 0`, and is evaluated at the local parameter `t - b_(i-1)`. A parameter exactly on a breakpoint belongs to the later segment, at its local parameter zero; the last segment also owns the end, `t = duration()`. A zero-duration segment other than the last is never evaluated.

Only the parameter is offset. The values are each segment's own: `piecewise_path` does not shift positions or match derivatives at a breakpoint, so a continuous composition is one whose segments end where the next begins.

## Usage Example

The elementary path functions are free functions, so each segment is a small struct. A rise, a dwell and a return make a rise-dwell-return cam law:

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'q', '' using 1:3 with lines title 'dq'"

#include "ctrlpp/trajectory/quintic_path.h"
#include "ctrlpp/trajectory/piecewise_path.h"

#include <iostream>

struct rise
{
    ctrlpp::path_point<double> evaluate(double tau) const { return ctrlpp::quintic_path(tau); }
    double duration() const { return 1.0; }
};

struct dwell
{
    ctrlpp::path_point<double> evaluate(double) const { return {1.0, 0.0, 0.0, 0.0}; }
    double duration() const { return 0.5; }
};

struct fall
{
    ctrlpp::path_point<double> evaluate(double tau) const
    {
        auto const p = ctrlpp::quintic_path(tau);
        return {1.0 - p.q, -p.dq, -p.ddq, -p.dddq};
    }
    double duration() const { return 1.0; }
};

int main()
{
    ctrlpp::piecewise_path<double, rise, dwell, fall> cam{rise{}, dwell{}, fall{}};
    for (double s = 0.0; s <= cam.duration(); s += 0.005) {
        auto const pp = cam.evaluate(s);
        std::cout << s << "," << pp.q << "," << pp.dq << "\n";
    }
}
```

## See Also

- [piecewise-trajectory](piecewise-trajectory.md)<br/> analogous composition for trajectories
- [path-segment](path-segment.md)<br/> concept each segment must satisfy
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
