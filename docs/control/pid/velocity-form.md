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

## Usage Example

```cpp
#include <ctrlpp/pid.h>

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

- [PID overview](README.md) -- parent PID documentation
- [anti-windup](anti-windup.md) -- explicit anti-windup for position form
- [isa-form](isa-form.md) -- ISA standard PID parameterisation
- [reference/pid-theory](../../reference/pid-theory.md) -- velocity form derivation
