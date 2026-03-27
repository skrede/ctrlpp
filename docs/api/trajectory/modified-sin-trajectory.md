# modified_sin_trajectory

Modified sinusoidal velocity profile with harmonic+cycloidal blend. Uses a sinusoidal central region with cycloidal transition regions at start and end. Produces very smooth acceleration with reduced vibration excitation compared to standard profiles. Duration is user-specified.

| Property | Value |
|----------|-------|
| **Header** | `ctrlpp/traj/modified_sin_trajectory.h` |

## Template Parameters

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type |

Output dimension is fixed at ND=1 (scalar trajectory).

## Config

```cpp
struct config {
    Scalar q0;  // start position
    Scalar q1;  // end position
    Scalar T;   // total duration (user-specified)
};
```

## Constructor

```cpp
explicit modified_sin_trajectory(config const& cfg);
```

## Profile Structure

Three regions spanning the full profile, symmetric about T/2:

| Region | Interval | Description |
|--------|----------|-------------|
| A-B | [0, T/8) | Cycloidal start |
| B-C | [T/8, 7T/8) | Sinusoidal middle |
| C-D | [7T/8, T] | Cycloidal end (symmetric to A-B) |

Normalized velocity profile (tau = t/T) with `A = pi/(pi+4)`:

- **A-B:** `v(tau) = A*(1 - cos(4*pi*tau))`
- **B-C:** `v(tau) = A*(1 + 3*sin(4*pi*tau/3 - pi/6))`
- **C-D:** `v(tau) = A*(1 - cos(4*pi*(1-tau)))` (symmetry)

## Key Constants

- **Peak velocity:** `4*pi*h / ((pi+4)*T)`
- **Peak acceleration:** `4*pi^2*h / ((pi+4)*T^2)`

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `evaluate` | `trajectory_point<Scalar, 1> evaluate(Scalar t) const` | Position, velocity, acceleration at time `t` |
| `duration` | `Scalar duration() const` | Total duration `T` |
| `peak_velocity` | `Scalar peak_velocity() const` | `4*pi*h / ((pi+4)*T)` |

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel', '' using 1:4 with lines title 'acc'"

#include "ctrlpp/traj/modified_sin_trajectory.h"

#include <iostream>

int main()
{
    ctrlpp::modified_sin_trajectory<double> traj({
        .q0 = 0.0, .q1 = 5.0, .T = 2.0
    });
    for (double t = 0; t <= traj.duration(); t += 0.01) {
        auto pt = traj.evaluate(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "," << pt.acceleration(0) << "\n";
    }
}
```

## See Also

- [modified-trap-trajectory](modified-trap-trajectory.md) -- alternative smooth profile
- [harmonic-path](harmonic-path.md) -- elementary harmonic motion law
- [cycloidal-path](cycloidal-path.md) -- elementary cycloidal motion law
- [Trajectory Generation Theory](../../background/trajectory-generation.md)
