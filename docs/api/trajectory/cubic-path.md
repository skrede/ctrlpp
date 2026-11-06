# cubic_path

Normalized cubic polynomial path: `q(tau) = 3*tau^2 - 2*tau^3`. C1-continuous with zero velocity at endpoints but non-zero acceleration.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/trajectory/cubic_path.h` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | none (`typename`) | Floating-point type |

## Functions

### cubic_path

```cpp
template <typename Scalar>
path_point<Scalar> cubic_path(Scalar tau);
```

Evaluate the cubic motion law at normalized time `tau`, meant for `[0, 1]`. Returns position, velocity, acceleration, and jerk, as derivatives with respect to `tau`. `tau` is not clamped: outside `[0, 1]` the polynomial is evaluated as written, so a caller clamps `tau` itself when it can leave the interval.

### cubic_path_peak_derivatives

```cpp
template <typename Scalar>
std::array<Scalar, 3> cubic_path_peak_derivatives();
```

Returns the peak magnitudes `{dq_max, ddq_max, dddq_max}` = `{1.5, 6.0, 12.0}` over `[0, 1]`, the argument `compute_min_duration` takes.

## Peak Derivatives

| Derivative | Value | Description |
|------------|-------|-------------|
| dq_max | 1.5 | Peak normalized velocity |
| ddq_max | 6.0 | Peak normalized acceleration |
| dddq_max | 12.0 | Peak normalized jerk (constant) |

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'q', '' using 1:3 with lines title 'dq', '' using 1:4 with lines title 'ddq'"

#include "ctrlpp/trajectory/cubic_path.h"

#include <iostream>

int main()
{
    for (double tau = 0; tau <= 1.0; tau += 0.005) {
        auto pp = ctrlpp::cubic_path(tau);
        std::cout << tau << "," << pp.q << "," << pp.dq << "," << pp.ddq << "\n";
    }
}
```

## See Also

- [cubic-trajectory](cubic-trajectory.md)<br/> polynomial with arbitrary velocity BCs
- [quintic-path](quintic-path.md)<br/> higher continuity (C2)
- [time-scaling](time-scaling.md)<br/> uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
