# anti_windup

Prevents integral windup when the PID output saturates against its limits. Three strategies are available, selectable at compile time via a strategy tag: back-calculation, clamping, and conditional integration. Each modifies how the integrator behaves when the controller output is clipped.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::anti_windup<Strategy>` | `#include <ctrlpp/control/pid_policies.h>` |

## Template Variants

```cpp
ctrlpp::anti_windup<ctrlpp::back_calc>                // back-calculation (default)
ctrlpp::anti_windup<ctrlpp::clamping>                 // integrator clamping
ctrlpp::anti_windup<ctrlpp::conditional_integration>  // conditional integration
```

## Config Fields

Fields added to `pid_config` when this policy is active:

### back_calc

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `kb` | `std::array<Scalar, N>` | `{}` (zeros) | Back-calculation gain per channel. Controls how fast the integrator is driven back when saturated. Typical value: 1/Ti. |

### clamping

No additional config fields. The integrator is simply frozen when the output is at its limit and the integrator would push it further.

### conditional_integration

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `error_threshold` | `std::array<Scalar, N>` | `{}` (zeros) | Integration is suspended when the absolute error exceeds this threshold per channel. |

## Behavior

- **back_calc**: When the PID output is clipped, the difference between the clipped and unclipped output is fed back through gain `kb` to adjust the integrator state. This provides smooth recovery from saturation.
- **clamping**: The integrator is held at its current value whenever the output is saturated and the error would increase the saturation. Simple and robust.
- **conditional_integration**: Integration is paused entirely when the error magnitude exceeds `error_threshold`. Useful for large setpoint changes where integration during the transient would cause overshoot.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'output', '' using 1:3 with lines title 'control'"

#include <ctrlpp/control/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::anti_windup<ctrlpp::back_calc>>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.0);
    cfg.output_min = Vec::Constant(-1.0);
    cfg.output_max = Vec::Constant(1.0);
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb = {1.0};

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 10.0; t += dt) {
        auto sp = Vec::Constant(5.0);   // large setpoint to trigger saturation
        auto meas = Vec::Constant(y);
        auto u = ctrl.compute(sp, meas, dt);
        y = 0.95 * y + 0.05 * u[0];
        std::cout << t << "," << y << "," << u[0] << "\n";
    }
}
```

## See Also

- [PID overview](README.md) -- parent PID documentation
- [velocity-form](velocity-form.md) -- bumpless transfer via incremental output
- [rate-limit](rate-limit.md) -- output rate limiting
- [reference/pid-theory](../../reference/pid-theory.md) -- anti-windup theory
