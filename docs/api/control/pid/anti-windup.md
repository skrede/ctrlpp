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
| `kb` | `std::array<Scalar, N>` | `{}` (zeros) | Back-calculation gain per channel (units of 1/time). Controls how fast the integrator is driven back when saturated. Leaving a channel at zero selects the automatic default described below. `create` refuses a `NaN` or infinite `kb` with `invalid_back_calculation_gain`, and a default it cannot form with `undefined_tracking_time` or `unrepresentable_gain`. |

The back-calculation gain multiplies the saturation error (the clipped minus unclipped output, in output units) into the integrator, whose rate is in output-per-time, so `kb` must carry units of 1/time. When a channel's `kb` is left at zero, it is auto-computed using the Astrom tracking-time-constant form `kb = 1/Tt` with `Tt = sqrt(Ti*Td)`:

- with derivative action (`kd != 0`): `kb = sqrt(ki/kd)`,
- without derivative action (`kd == 0`): the tracking time collapses to the integral time, giving the fallback `kb = ki/kp`,
- for a pure integral controller (`kp == 0` and `kd == 0`): there is no proportional or derivative reference for a tracking time, so back-calculation is disabled (`kb = 0`).

The gains above are the internal parallel gains (in ISA form `ki = Kp/Ti` and `kd = Kp*Td`). The tracking time exists only when `Ti = kp/ki` and `Td = kd/kp` are positive, which holds exactly when the nonzero gains share a sign. A channel with integral action (`ki != 0`) whose nonzero gains do not share a sign has no default: `sqrt(ki/kd)` is not a number and `ki/kp` is negative, so `create` and `set_params` refuse it with `undefined_tracking_time` while its `kb` is left at zero. An explicit nonzero `kb` is accepted for such a channel. A reverse-acting controller, with every gain negative, shares a sign and gets the positive default. A default that overflows `Scalar` (for example `ki/kd` beyond the largest finite value) is refused with `unrepresentable_gain`.

### clamping

No additional config fields. The integrator is simply frozen when the output is at its limit and the integrator would push it further.

### conditional_integration

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `error_threshold` | `std::array<Scalar, N>` | `{}` (zeros) | Integration is suspended when the absolute error exceeds this threshold per channel. `create` refuses a `NaN` threshold with `invalid_error_threshold`. |

## Behavior

- **back_calc**: The difference between the applied output and the unconstrained command is fed back through gain `kb` to adjust the integrator state. The unconstrained command is taken before both the rate limiter and the output clamp, so a rate-limited ramp also engages the anti-windup, not only output saturation. This provides smooth recovery from saturation.
- **clamping**: The integrator is held whenever the applied output is constrained and the error would increase the constraint. The decision is made per channel, so saturation in one channel of a MIMO controller does not freeze the integrators of unsaturated channels. Simple and robust.
- **conditional_integration**: Integration is paused entirely when the error magnitude exceeds `error_threshold`. Useful for large setpoint changes where integration during the transient would cause overshoot.

A channel counts as constrained when the rate limiter clipped it or the output clamp saturated it. The verdict comes from the limiter and the clamp themselves, never from comparing the applied output with the unconstrained command, so a composed rate limiter that does not clip leaves anti-windup exactly as it is without one.

`anti_windup` does not compose with `velocity_form`: the combination is a compile error whose message states that clamping the accumulated output already prevents windup. See [velocity-form](velocity-form.md).

While `freeze_integral` holds the integrator, no strategy changes it, back-calculation included: a frozen integrator cannot wind up, so the explicit hold takes precedence. Back-calculation resumes on the first cycle after the integrator is unfrozen.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'output', '' using 1:3 with lines title 'control'"

#include <ctrlpp/control/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, ctrlpp::anti_windup<ctrlpp::back_calc>>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.0);
    cfg.output_min = Vec::Constant(-1.0);
    cfg.output_max = Vec::Constant(1.0);
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb = {1.0};

    auto created = Pid::create(cfg);
    if (!created.has_value()) return 1;
    auto& ctrl = *created;

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 10.0; t += dt) {
        auto sp = Vec::Constant(5.0);   // large setpoint to trigger saturation
        auto meas = Vec::Constant(y);
        // A refused cycle produced no command; the caller decides what the
        // actuator does. This sample stops.
        auto step = ctrl.compute(sp, meas, dt);
        if (!step.has_value()) return 1;
        const auto& u = *step;
        y = 0.95 * y + 0.05 * u[0];
        std::cout << t << "," << y << "," << u[0] << "\n";
    }
}
```

## See Also

- [PID overview](README.md)<br/> parent PID documentation
- [velocity-form](velocity-form.md)<br/> incremental output, which refuses anti_windup
- [rate-limit](rate-limit.md)<br/>output rate limiting
- [guides/pid/composition](../../../guides/pid/composition.md)<br/>composing policies
- [background/pid](../../../background/pid.md)<br/>anti-windup theory
