# feed_forward

Adds a feed-forward term to the PID output. Feed-forward injects a control action based on the setpoint (or a known disturbance model) directly into the output, bypassing the feedback loop. This improves disturbance rejection and setpoint tracking speed by not relying solely on the error signal. The feed-forward function is provided as a callable type in the template argument and stored in the configuration.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::feed_forward<Callable>` | `#include <ctrlpp/control/pid_policies.h>` |
| `ctrlpp::feed_forward<void>` | `#include <ctrlpp/control/pid_policies.h>` (no-op variant) |

## Config Fields

### feed_forward\<Callable\>

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `feed_forward_cb` | `Callable` | value-initialized | Called once per cycle as `feed_forward_cb(sp, dt)`, with the setpoint passed to `compute` (before `setpoint_filter`) and the cycle's step; it returns the feed-forward contribution as a `vector_t` or an Eigen expression that converts to one. Any other input it needs, such as a measured disturbance, it holds itself. |

### feed_forward\<void\>

No additional config fields. The void specialization, also the default argument `feed_forward<>`, is a no-op placeholder: nothing is added to the output.

## Behavior

The feed-forward callable is evaluated once per cycle that passes the input checks. In position form its output is added to the PID feedback output before the rate limiter and the output clamp, so both act on the sum, and it counts with the proportional and derivative terms in the integrator the tracking overload of `compute` re-seeds. In velocity form the increment carries the change in the feed-forward value since the previous cycle, so a constant feed-forward adds nothing after the first cycle (see [velocity-form](velocity-form.md)). This allows the controller to anticipate the required control action for known reference trajectories or measurable disturbances, reducing the error signal that the PID terms must correct.

A `NaN` feed-forward value makes the cycle fail with `non_finite_result`, and the controller is rolled back. An infinite one fails the cycle in velocity form, where the feed-forward value is carried to the next cycle, and with `anti_windup<back_calc>`, whose feedback it makes non-finite; otherwise the position form clamps it to the output limit like any other command and the cycle succeeds. A rollback restores the controller, not the callable: a stateful callable keeps whatever it changed during a cycle that was then refused.

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

    using Pid = ctrlpp::pid<double, 1, ctrlpp::feed_forward<FF>>;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(0.5);
    cfg.kd = Vec::Constant(0.0);
    cfg.template policy<ctrlpp::feed_forward<FF>>().feed_forward_cb = ff;

    auto created = Pid::create(cfg);
    if (!created.has_value()) return 1;
    auto& ctrl = *created;

    double y = 0.0;
    constexpr double dt = 0.01;

    for (int k = 0; k < 500; ++k)
    {
        double t = k * dt;
        double r = (t >= 1.0) ? 1.0 : 0.0;
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

## Limitations

### set_params needs a copy-assignable callable

`set_params` copy-assigns the whole configuration, so it compiles only when `Callable` is copy-assignable. A lambda without captures is, and so is a function object with copyable members; a lambda with captures is not, and with it `set_params` is a compile error while `create` and `compute` work. `Pid::config_type cfg{}` also needs a default-constructible `Callable`, which a lambda with captures is not either; such a lambda is passed in when the configuration is built, as `Pid::config_type cfg{.policies = {{ff}}}`.

## See Also

- [PID overview](README.md)<br/> parent PID documentation
- [setpoint-filter](setpoint-filter.md)<br/> reference shaping alternative
- [guides/pid/composition](../../../guides/pid/composition.md)<br/>composing policies
- [background/pid](../../../background/pid.md)<br/>feed-forward design theory
