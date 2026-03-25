# isa_form

Switches the PID to ISA (Instrument Society of America) standard form, where the controller is parameterised with proportional gain Kp, integral time Ti, and derivative time Td instead of the independent gains Kp, Ki, Kd. The ISA form uses u = Kp * (e + (1/Ti) * integral(e) + Td * de/dt), tying integral and derivative action to the proportional gain.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::isa_form` | `#include <ctrlpp/control/pid_policies.h>` |

## Config Fields

When `isa_form` is active, the `ki` and `kd` fields in `pid_config` are reinterpreted:

| Field | Interpretation | Description |
|-------|---------------|-------------|
| `kp` | Kp | Proportional gain (same as standard form) |
| `ki` | Ti | Integral time constant. Actual integral gain = Kp / Ti. |
| `kd` | Td | Derivative time constant. Actual derivative gain = Kp * Td. |

No additional config struct fields are added. The ISA form modifies how the existing `kp`, `ki`, `kd` fields are used internally.

## Behavior

The ISA form couples the integral and derivative actions to the proportional gain. When Kp is changed (e.g., for gain scheduling), the integral and derivative actions scale proportionally. This is the standard industrial tuning convention and makes tuning rules (Ziegler-Nichols, Cohen-Coon, etc.) directly applicable.

## Usage Example

```cpp
#include <ctrlpp/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::isa_form>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);   // proportional gain
    cfg.ki = Vec::Constant(0.5);   // Ti = 0.5s (integral time)
    cfg.kd = Vec::Constant(0.1);   // Td = 0.1s (derivative time)

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
- [velocity-form](velocity-form.md) -- incremental output form
- [performance](performance.md) -- controller performance metrics
- [guides/pid/composition](../../guides/pid/composition.md) -- composing policies
- [reference/pid-theory](../../reference/pid-theory.md) -- ISA standard form derivation
