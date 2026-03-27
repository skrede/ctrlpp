# path_segment

Concept constraining types that represent a normalized path segment evaluated over [0,1]. A path segment provides `evaluate(tau)` returning a `path_point` and `duration()` returning the normalized duration (typically 1.0).

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/path_segment.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `S` | The type being constrained |
| `Scalar` | Floating-point type |

## Concept Definition

```cpp
template <typename S, typename Scalar>
concept path_segment = requires(const S& seg, Scalar tau) {
    { seg.evaluate(tau) } -> std::convertible_to<path_point<Scalar>>;
    { seg.duration() } -> std::convertible_to<Scalar>;
};
```

## Required Expressions

| Expression | Return type | Description |
|------------|-------------|-------------|
| `seg.evaluate(tau)` | convertible to `path_point<Scalar>` | Normalized position and derivatives at `tau` |
| `seg.duration()` | convertible to `Scalar` | Normalized duration (typically 1.0) |

## Relationship to Free Functions

The elementary path functions (`cubic_path`, `quintic_path`, etc.) are free functions, not types. They produce `path_point` values but do not satisfy this concept directly. Use `piecewise_path` for composing path segments that satisfy this concept.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'q', '' using 1:3 with lines title 'dq'"

#include "ctrlpp/trajectory/cubic_path.h"

#include <iostream>

int main()
{
    for (double tau = 0; tau <= 1.0; tau += 0.005) {
        auto pp = ctrlpp::cubic_path(tau);
        std::cout << tau << "," << pp.q << "," << pp.dq << "\n";
    }
}
```

## See Also

- [trajectory-segment](trajectory-segment.md) -- analogous concept for physical trajectories
- [trajectory-types](trajectory-types.md) -- `path_point` struct
- [piecewise-path](piecewise-path.md) -- variadic composition of path segments
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
