# velocity_form

Switches the PID from position form (output = absolute value) to velocity form (output = incremental change). In velocity form, the controller outputs delta-u at each step, and the actuator integrates. The controller keeps no integrator of its own: it accumulates the increments it emits so it can clamp that running sum to the output limits, which avoids windup by construction. A switch from manual to automatic is bumpless when the tracking overload of `compute` re-seeds that sum, as described below.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::velocity_form` | `#include <ctrlpp/control/pid_policies.h>` |

## Config Fields

No additional config fields. The velocity form policy modifies the internal computation structure without requiring extra parameters.

## Behavior

Instead of computing u(k) directly, the controller computes the change delta_u(k) = u(k) - u(k-1). The actual output is accumulated externally (or by the actuator). This means:

- There is no integrator separate from the accumulated output, the running sum of the emitted increments that the output clamp acts on. `integral()` stays zero, and `set_integral` is not available: calling it on a velocity-form controller is a compile error.
- After `create` or `reset` the previous proportional term is formed from a zero setpoint and a zero measurement, and the accumulated output is zero. The first increment therefore carries the full proportional term `kp*(b*setpoint - measurement)` plus the integral increment, and the output limits act on a sum that starts at zero, not at the actuator's position. Switching from manual to automatic is bumpless only through tracking: while in manual, run the tracking overload of `compute` with the manual output on every cycle. Each such cycle re-seeds the accumulated output to the manual output clamped to the limits and keeps the previous setpoint and measurement current, so the first automatic cycle emits only the change in the proportional and derivative terms plus the integral increment, accumulated from the manual output. This replaces `set_integral` as the way to seed a velocity-form controller.
- The increment is the difference of the position form's own proportional and derivative terms plus the integral increment: `delta_u(k) = (P(k) - P(k-1)) + delta_I(k) + (D(k) - D(k-1))`. `P(k-1)` is recomputed from the previous filtered setpoint and measurement with the current gains, and `D(k-1)` is the derivative term of the previous cycle. When no limit binds, the running sum of the increments equals the position form's output.
- The first step emits only the proportional and integral increments. The derivative term is zero on the first step, as in the position form, so step 2 emits a first difference and later steps a second difference. A setpoint step does not kick the output and leaves no parked offset.
- The derivative acts on the measurement unless `derivative_on_error` is set, as in the position form. With `derivative_on_error` set it differentiates `c*setpoint - measurement`.
- The setpoint weights `b` and `c`, `deriv_filter`, the integral discretization (the default `ki*e(k)*dt`, or `forward_euler` and `tustin` when composed), `setpoint_filter`, `pv_filter` and `freeze_integral` are honored. While the integral is frozen the integral increment is zero.
- Non-uniform sampling is exact: each cycle's derivative divides its own measurement difference by its own `dt`, so a measurement ramp sampled at a varying rate produces no spurious derivative increment.
- `output_min` and `output_max` bound the accumulated output, not the individual increment. The emitted increment is the change needed to bring the accumulated output to that clamped value, so the output can both rise and fall against asymmetric limits (for example a valve with `output_min = 0`). `saturated()` is true exactly on the cycles where this clamp changed the accumulated output.
- Feed-forward is injected incrementally: the emitted increment carries the change in the feed-forward level, so a constant feed-forward adds nothing to the increment and the actuator does not drift.
- With `rate_limit` composed, each channel's increment, feed-forward increment included, is clipped to `+-rate_max * dt` before it is accumulated. A channel within its bound, or with `rate_max` at zero or `+inf`, passes its increment through bit for bit. The clipped part of an increment is dropped, not carried into later cycles: the controller does not catch up on a change the limiter refused.
- The tracking overload of `compute` re-seeds the accumulated output to the tracking signal clamped to `output_min` and `output_max`, after the cycle succeeds. It returns that cycle's computed increment, and the next cycle accumulates from the re-seeded value. `integral()` stays zero.
- Composing `velocity_form` with any `anti_windup` strategy is a compile error: clamping the accumulated output already prevents windup. The accumulated output never moves past a limit, so inhibiting integration while the output saturates, which is equivalent to back-calculation (Astrom, *Control System Design*, 2002, Sec. 6.5), is built into the form.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'output', '' using 1:3 with lines title 'delta_u', '' using 1:4 with lines title 'u_accum'"

#include <ctrlpp/control/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, ctrlpp::velocity_form>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.1);

    auto created = Pid::create(cfg);
    if (!created.has_value()) return 1;
    auto& ctrl = *created;

    double y = 0.0;
    double u_accum = 0.0;  // actuator integrates the incremental output
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 5.0; t += dt) {
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y);
        // A refused cycle produced no command; the caller decides what the
        // actuator does. This sample stops.
        auto step = ctrl.compute(sp, meas, dt);
        if (!step.has_value()) return 1;
        const auto& delta_u = *step;
        u_accum += delta_u[0];  // actuator integration
        y = 0.9 * y + 0.1 * u_accum;
        std::cout << t << "," << y << "," << delta_u[0] << "," << u_accum << "\n";
    }
}
```

## See Also

- [PID overview](README.md)<br/> parent PID documentation
- [anti-windup](anti-windup.md)<br/> anti-windup strategies for the position form; they do not compose with velocity_form
- [isa-form](isa-form.md)<br/>ISA standard PID parametrization
- [guides/pid/composition](../../../guides/pid/composition.md)<br/>composing policies
- [background/pid](../../../background/pid.md)<br/>velocity form derivation
