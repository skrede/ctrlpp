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
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'reference', '' using 1:3 with lines title 'output', '' using 1:4 with lines title 'control'"

#include <ctrlpp/control/pid.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    using Vec = Eigen::Matrix<double, 1, 1>;

    // Feed-forward function: scales setpoint by inverse plant gain.
    // Signature must be (const vector_t&, Scalar) -> vector_t.
    auto ff = [](const Vec& sp, double /*dt*/) -> Vec {
        return Vec::Constant(0.5 * sp[0]);
    };
    using FF = decltype(ff);

    using Pid = ctrlpp::pid<double, 1, 1, 1, ctrlpp::feed_forward<FF>>;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(0.5);
    cfg.kd = Vec::Constant(0.0);
    cfg.template policy<ctrlpp::feed_forward<FF>>().ff_func = ff;

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double dt = 0.01;

    for (int k = 0; k < 500; ++k)
    {
        double t = k * dt;
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
- [reference/pid-theory](../../reference/pid-theory.md) -- feed-forward design theory
