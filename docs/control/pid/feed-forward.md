# feed_forward

Adds a feed-forward term to the PID output. Feed-forward injects a control action based on the setpoint (or a known disturbance model) directly into the output, bypassing the feedback loop. This improves disturbance rejection and setpoint tracking speed by not relying solely on the error signal. The feed-forward function is provided as a callable template parameter.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::feed_forward<Callable>` | `#include <ctrlpp/control/pid_policies.h>` |
| `ctrlpp::feed_forward<void>` | `#include <ctrlpp/control/pid_policies.h>` (no-op variant) |

## Config Fields

### feed_forward\<Callable\>

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `ff_func` | `Callable` | -- | A callable that computes the feed-forward contribution. Signature depends on the user's model. |

### feed_forward\<void\>

No additional config fields. The void specialisation is a no-op placeholder.

## Behavior

The feed-forward callable is evaluated at each control step and its output is added to the PID feedback output. This allows the controller to anticipate the required control action for known reference trajectories or measurable disturbances, reducing the error signal that the PID terms must correct.

## Usage Example

```cpp
#include <ctrlpp/pid.h>

#include <iostream>

int main()
{
    // Feed-forward function: scales setpoint by inverse plant gain
    auto ff = [](double setpoint) { return 0.5 * setpoint; };
    using FF = decltype(ff);

    using Pid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::feed_forward<FF>>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(0.5);
    cfg.kd = Vec::Constant(0.0);
    cfg.template policy<ctrlpp::feed_forward<FF>>().ff_func = ff;

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 5.0; t += dt) {
        double r = (t >= 1.0) ? 1.0 : 0.0;
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
- [setpoint-filter](setpoint-filter.md) -- reference shaping alternative
- [guides/pid/composition](../../guides/pid/composition.md) -- composing policies
- [reference/pid-theory](../../reference/pid-theory.md) -- feed-forward design theory
