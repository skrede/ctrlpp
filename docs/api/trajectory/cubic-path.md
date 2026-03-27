# cubic_path

Normalized cubic polynomial path: `q(tau) = 3*tau^2 - 2*tau^3`. C1-continuous with zero velocity at endpoints but non-zero acceleration.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/cubic_path.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

## Functions

### cubic_path

```cpp
template <typename Scalar>
auto cubic_path(Scalar tau) -> path_point<Scalar>;
```

Evaluate the cubic motion law at normalized time `tau` in [0,1]. Returns position, velocity, acceleration, and jerk.

### cubic_path_peak_derivatives

```cpp
template <typename Scalar>
auto cubic_path_peak_derivatives() -> std::array<Scalar, 3>;
```

Returns `{dq_max, ddq_max, dddq_max}` = `{1.5, 6.0, 12.0}`. Used by `compute_min_duration` for time scaling.

## Peak Derivatives

| Derivative | Value | Description |
|------------|-------|-------------|
| dq_max | 1.5 | Peak normalized velocity |
| ddq_max | 6.0 | Peak normalized acceleration |
| dddq_max | 12.0 | Peak normalized jerk (constant) |

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'q', '' using 1:3 with lines title 'dq', '' using 1:4 with lines title 'ddq'"

#include "ctrlpp/traj/cubic_path.h"

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

- [cubic-trajectory](cubic-trajectory.md) -- polynomial with arbitrary velocity BCs
- [quintic-path](quintic-path.md) -- higher continuity (C2)
- [time-scaling](time-scaling.md) -- uses peak derivatives for duration computation
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
