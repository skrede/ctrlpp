# rate_limit

Limits the rate of change of the PID output between consecutive time steps. This protects actuators from excessive slew rates and reduces mechanical stress. The output change per step is clamped to +/- rate_max, regardless of what the PID computation requests.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::rate_limit` | `#include <ctrlpp/control/pid_policies.h>` |

## Config Fields

Fields added to `pid_config` when this policy is active:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `rate_max` | `std::array<Scalar, N>` | `{}` (zeros) | Maximum allowable output change per time step, per channel. Units are output-units per sample period. |

## Behavior

After the PID output is computed, the change from the previous output is clamped: delta_u = clamp(u_new - u_prev, -rate_max, +rate_max). The actual output becomes u_prev + delta_u. This is applied independently per channel.

Rate limiting interacts with anti-windup: if both are active, anti-windup sees the rate-limited (actual) output as the saturated value.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'output', '' using 1:3 with lines title 'control'"

#include <ctrlpp/control/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::rate_limit>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(10.0);   // aggressive gain
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.0);
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {0.1};  // max 0.1 per step

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 5.0; t += dt) {
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y);
        auto u = ctrl.compute(sp, meas, dt);
        y = 0.9 * y + 0.1 * u[0];
        std::cout << t << "," << y << "," << u[0] << "\n";
    }
}
```

## See Also

- [PID overview](README.md) -- parent PID documentation
- [anti-windup](anti-windup.md) -- integrator saturation handling
- [guides/pid/composition](../../guides/pid/composition.md) -- composing policies
- [reference/pid-theory](../../reference/pid-theory.md) -- rate limiting theory
