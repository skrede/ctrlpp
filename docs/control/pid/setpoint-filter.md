# setpoint_filter / pv_filter

Two complementary filtering policies for shaping the PID reference path. `setpoint_filter` applies setpoint weighting (b, c parameters) to reduce overshoot on setpoint changes without affecting disturbance rejection. `pv_filter` applies a first-order low-pass filter to the process variable before it enters the PID computation.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::setpoint_filter` | `#include <ctrlpp/control/pid_policies.h>` |
| `ctrlpp::pv_filter` | `#include <ctrlpp/control/pid_policies.h>` |

## Config Fields

### setpoint_filter

Uses the base `pid_config` fields `b` and `c` (already present in all PID configurations):

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `b` | `Vector<Scalar, N>` | `1.0` (all channels) | Proportional setpoint weight. The proportional term uses b*r - y instead of r - y. |
| `c` | `Vector<Scalar, N>` | `1.0` (all channels) | Derivative setpoint weight. The derivative term uses c*r - y instead of r - y. |

Additionally, the policy adds a filter time constant:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `tf` | `std::array<Scalar, N>` | `{}` (zeros) | Setpoint prefilter time constant per channel. |

### pv_filter

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `tf` | `std::array<Scalar, N>` | `{}` (zeros) | Process variable filter time constant per channel. Larger values mean heavier smoothing. |

## Behavior

- **setpoint_filter**: Weights b and c scale the setpoint contribution to the proportional and derivative terms respectively. Setting b < 1 reduces the proportional kick on setpoint changes. Setting c = 0 eliminates derivative kick entirely. The filter time constant `tf` applies a first-order prefilter to the setpoint signal.
- **pv_filter**: Applies a first-order low-pass filter to the process variable measurement before it enters the PID error computation. Useful for reducing measurement noise without modifying the derivative filter coefficient.

## Usage Example

```cpp
#include <ctrlpp/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::setpoint_filter>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(4.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.5);
    cfg.b = Vec::Constant(0.7);   // reduce proportional kick
    cfg.c = Vec::Constant(0.0);   // no derivative kick

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 5.0; t += dt) {
        double r = (t >= 1.0) ? 1.0 : 0.0;  // step at t=1
        auto sp = Vec::Constant(r);
        auto meas = Vec::Constant(y);
        auto u = ctrl.compute(sp, meas, dt);
        y = 0.9 * y + 0.1 * u[0];
        std::cout << t << "," << r << "," << y << "," << u[0] << "\n";
    }
}
```

## See Also

- [PID overview](README.md) -- parent PID documentation
- [derivative-filter](derivative-filter.md) -- noise attenuation on the derivative term
- [guides/pid/composition](../../guides/pid/composition.md) -- composing policies
- [reference/pid-theory](../../reference/pid-theory.md) -- setpoint weighting theory
