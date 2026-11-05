# setpoint_filter / pv_filter

Two first-order low-pass filters on the PID inputs. `setpoint_filter` filters the setpoint, which softens the response to setpoint steps without affecting disturbance rejection. `pv_filter` filters the process variable before it enters the PID computation. Setpoint weighting, the `b` and `c` fields, is part of every `pid_config` and needs neither policy; it is described here because it shapes the same reference path.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::setpoint_filter` | `#include <ctrlpp/control/pid_policies.h>` |
| `ctrlpp::pv_filter` | `#include <ctrlpp/control/pid_policies.h>` |

## Config Fields

### Setpoint weights (base fields)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `b` | `Vector<Scalar, NY>` | `1.0` (all channels) | Proportional setpoint weight. The proportional term uses b*r - y instead of r - y. `create` refuses a `NaN` or infinite component with `non_finite_setpoint_weight`. |
| `c` | `Vector<Scalar, NY>` | `1.0` (all channels) | Derivative setpoint weight. With `derivative_on_error` set, the derivative term uses c*r - y instead of r - y; without it the derivative acts on -y alone and `c` is not read. `create` refuses a `NaN` or infinite component with `non_finite_setpoint_weight`. |

### setpoint_filter

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `tf` | `std::array<Scalar, N>` | `{}` (zeros) | Setpoint prefilter time constant per channel, in the time unit of `dt`. Zero leaves the channel unfiltered. `create` refuses a `NaN`, negative or `+inf` `tf` with `invalid_filter_time_constant`. |

### pv_filter

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `tf` | `std::array<Scalar, N>` | `{}` (zeros) | Process variable filter time constant per channel, in the time unit of `dt`. Larger values mean heavier smoothing. Zero leaves the channel unfiltered. `create` refuses a `NaN`, negative or `+inf` `tf` with `invalid_filter_time_constant`. |

## Behavior

- **Setpoint weights**: b and c scale the setpoint contribution to the proportional and derivative terms respectively. Setting b < 1 reduces the proportional kick on setpoint changes. The derivative acts on the measurement by default, which already avoids a derivative kick; with `derivative_on_error` set, c = 0 restores that. The integral always acts on the unweighted error r - y.
- **setpoint_filter**: each channel with a positive `tf` is filtered as `r_f = alpha * r_f + (1 - alpha) * r` with `alpha = tf / (tf + dt)`. The filtered setpoint replaces the setpoint in the error, the proportional and derivative terms and the performance metrics; the feed-forward callable still receives the unfiltered one.
- **pv_filter**: the same first-order filter on the measurement, which then replaces it in the error, the proportional and derivative terms and the metrics. Useful for reducing measurement noise without modifying the derivative filter coefficient.

Both filter states start at zero after `create` and `reset`, so the first cycles see a setpoint or measurement rising from zero, not the signal's current value.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'reference', '' using 1:3 with lines title 'output', '' using 1:4 with lines title 'control'"

#include <ctrlpp/control/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, ctrlpp::setpoint_filter>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(4.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.5);
    cfg.b = Vec::Constant(0.7);   // reduce proportional kick
    cfg.template policy<ctrlpp::setpoint_filter>().tf = {0.2};  // soften the step

    auto created = Pid::create(cfg);
    if (!created.has_value()) return 1;
    auto& ctrl = *created;

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 5.0; t += dt) {
        double r = (t >= 1.0) ? 1.0 : 0.0;  // step at t=1
        auto sp = Vec::Constant(r);
        auto meas = Vec::Constant(y);
        // A refused cycle produced no command; the caller decides what the
        // actuator does. This sample stops.
        auto step = ctrl.compute(sp, meas, dt);
        if (!step.has_value()) return 1;
        const auto& u = *step;
        y = 0.9 * y + 0.1 * u[0];
        std::cout << t << "," << r << "," << y << "," << u[0] << "\n";
    }
}
```

## See Also

- [PID overview](README.md)<br/> parent PID documentation
- [derivative-filter](derivative-filter.md)<br/>noise attenuation on the derivative term
- [guides/pid/composition](../../../guides/pid/composition.md)<br/>composing policies
- [background/pid](../../../background/pid.md)<br/>setpoint weighting theory
