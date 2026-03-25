# deriv_filter

First-order low-pass filter on the derivative term to attenuate high-frequency noise amplification. Without filtering, the derivative action amplifies measurement noise proportionally to frequency. The filter coefficient N controls the bandwidth: larger N allows more high-frequency content through (less filtering), while smaller N provides heavier smoothing at the cost of slower derivative response.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::deriv_filter` | `#include <ctrlpp/control/pid_policies.h>` |

## Config Fields

Fields added to `pid_config` when this policy is active:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `n` | `std::array<Scalar, N>` | `{}` (zeros) | Derivative filter coefficient per channel. Typical range: 3--20. Higher values mean less filtering. |

## Behavior

Replaces the pure derivative s*Kd with the filtered form Kd*N*s / (1 + N/s), implemented in discrete time. This limits the derivative gain at high frequencies to Kd*N rather than letting it grow unbounded.

## Usage Example

```cpp
#include <ctrlpp/pid.h>

#include <cmath>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::deriv_filter>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(0.5);
    cfg.kd = Vec::Constant(0.1);
    cfg.template policy<ctrlpp::deriv_filter>().n = {10.0};

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 5.0; t += dt) {
        // Simulate noisy measurement
        double noise = 0.05 * std::sin(100.0 * t);
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y + noise);
        auto u = ctrl.compute(sp, meas, dt);
        y = 0.9 * y + 0.1 * u[0];
        std::cout << t << "," << y << "," << u[0] << "\n";
    }
}
```

## See Also

- [PID overview](README.md) -- parent PID documentation
- [setpoint-filter](setpoint-filter.md) -- reference signal filtering
- [guides/pid/composition](../../guides/pid/composition.md) -- composing policies
- [reference/pid-theory](../../reference/pid-theory.md) -- derivative filtering theory
