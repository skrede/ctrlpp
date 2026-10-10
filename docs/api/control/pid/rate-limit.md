# rate_limit

Limits the rate of change of the PID output between consecutive time steps. This protects actuators from excessive slew rates and reduces mechanical stress. The output change over a step of length dt is clamped to +/- rate_max * dt, regardless of what the PID computation requests.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::rate_limit` | `#include <ctrlpp/control/pid_policies.h>` |

## Config Fields

Fields added to `pid_config` when this policy is active:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `rate_max` | `std::array<Scalar, N>` | `{}` (zeros) | Maximum rate of change of the output per channel, in output units per second; a cycle of step `dt` moves the output by at most `rate_max * dt`. Zero and `+inf` leave the channel unlimited. `create` refuses a `NaN` or negative `rate_max` with `invalid_rate_limit`. |

## Behavior

The limiter runs per channel, after the PID command is formed and before the output clamp. `rate_max` is in output units per second, so the bound for one cycle is `rate_max * dt`. When the requested change `u_new - u_prev` exceeds `+rate_max * dt`, the output becomes `u_prev + rate_max * dt`; when it falls below `-rate_max * dt`, it becomes `u_prev - rate_max * dt`; either way the channel is reported as clipped. A channel within its bound, or with `rate_max` at zero or `+inf`, passes the command through bit for bit, so composing a limiter that never binds changes neither the output nor the integrator.

Rate limiting interacts with anti-windup: the constraint verdict is the limiter's own clip report together with the output clamp's, per channel, so a channel the limiter did not clip never engages clamping anti-windup. Back-calculation feeds back against the command taken before the limiter, so a rate-limited ramp drives the integrator back just as output saturation does.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'output', '' using 1:3 with lines title 'control'"

#include <ctrlpp/control/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, ctrlpp::rate_limit>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(10.0);   // aggressive gain
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.0);
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {0.1};  // max 0.1 output units per second

    auto created = Pid::create(cfg);
    if (!created.has_value()) return 1;
    auto& ctrl = *created;

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 5.0; t += dt) {
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y);
        // A refused cycle produced no command; the caller decides what the
        // actuator does. This sample stops.
        auto step = ctrl.compute(sp, meas, dt);
        if (!step.has_value()) return 1;
        const auto& u = *step;
        y = 0.9 * y + 0.1 * u[0];
        std::cout << t << "," << y << "," << u[0] << "\n";
    }
}
```

## See Also

- [PID overview](README.md)<br/> parent PID documentation
- [anti-windup](anti-windup.md)<br/>integrator saturation handling
- [guides/pid/composition](../../../guides/pid/composition.md)<br/>composing policies
- [background/pid](../../../background/pid.md)<br/>rate limiting theory
