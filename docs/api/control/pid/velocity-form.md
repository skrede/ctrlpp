# velocity_form

Switches the PID from position form (output = absolute value) to velocity form (output = incremental change). In velocity form, the controller outputs delta-u at each step, and the actuator integrates. This provides automatic bumpless transfer when switching between manual and automatic modes, and avoids integrator windup by construction since the integrator is in the actuator, not the controller.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::velocity_form` | `#include <ctrlpp/control/pid_policies.h>` |

## Config Fields

No additional config fields. The velocity form policy modifies the internal computation structure without requiring extra parameters.

## Behavior

Instead of computing u(k) directly, the controller computes the change delta_u(k) = u(k) - u(k-1). The actual output is accumulated externally (or by the actuator). This means:

- No explicit integrator state is maintained internally, eliminating integrator windup.
- Switching between manual and automatic control is bumpless: the incremental output starts from zero regardless of the previous manual setting.
- The derivative term uses backward differences of the error rather than the filtered derivative.
- `output_min` and `output_max` bound the accumulated output, not the individual increment. The emitted increment is the change needed to bring the accumulated output to that clamped value, so the output can both rise and fall against asymmetric limits (for example a valve with `output_min = 0`).
- Feed-forward is injected incrementally: the emitted increment carries the change in the feed-forward level, so a constant feed-forward adds nothing to the increment and the actuator does not drift.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'output', '' using 1:3 with lines title 'delta_u', '' using 1:4 with lines title 'u_accum'"

#include <ctrlpp/control/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::velocity_form>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.1);

    Pid ctrl(cfg);

    double y = 0.0;
    double u_accum = 0.0;  // actuator integrates the incremental output
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 5.0; t += dt) {
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y);
        auto delta_u = ctrl.compute(sp, meas, dt);
        u_accum += delta_u[0];  // actuator integration
        y = 0.9 * y + 0.1 * u_accum;
        std::cout << t << "," << y << "," << delta_u[0] << "," << u_accum << "\n";
    }
}
```

## See Also

- [PID overview](README.md)<br/> parent PID documentation
- [anti-windup](anti-windup.md)<br/> explicit anti-windup for position form
- [isa-form](isa-form.md)<br/>ISA standard PID parametrization
- [guides/pid/composition](../../../guides/pid/composition.md)<br/>composing policies
- [background/pid](../../../background/pid.md)<br/>velocity form derivation
