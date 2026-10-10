# deriv_filter

First-order low-pass filter on the derivative term to attenuate high-frequency noise amplification. Without filtering, the derivative action amplifies measurement noise proportionally to frequency. The filter time constant Tf sets the bandwidth. By default it is derived from the divisor N as Tf = Td/N: larger N allows more high-frequency content through (less filtering), while smaller N provides heavier smoothing at the cost of slower derivative response. A channel can instead state Tf directly in seconds.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::deriv_filter` | `#include <ctrlpp/control/pid_policies.h>` |

## Config Fields

Fields added to `pid_config` when this policy is active:

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `n` | `std::array<Scalar, N>` | `{}` (zeros) | Derivative filter divisor N per channel, giving the default time constant Tf = kd/(kp*N). Typical values are 8 to 20. Higher values mean less filtering. `create` refuses a `NaN` or infinite `n` with `invalid_derivative_filter`. |
| `tf` | `std::array<std::optional<Scalar>, N>` | `{}` (absent) | Explicit filter time constant Tf per channel, in seconds. When present it replaces the default kd/(kp*N) and `n` is not used for that channel. It must be finite and positive; `create` refuses a `NaN`, infinite, zero or negative `tf` with `invalid_derivative_filter`. |

The gains in this page are the internal parallel gains (in ISA form `kd = Kp*Td`). A channel without derivative action (`kd == 0`) is not filtered, and any finite `n` is accepted for it. A channel with derivative action and no `tf` must have a default time constant that is a positive time, so `create` and `set_params` refuse it with `undefined_derivative_filter` when:

- `kp == 0` (a D-only or I-D channel), since Tf = kd/(kp*N) has no value,
- `n <= 0`, which includes the default `n` of zero,
- `kd` and `kp` differ in sign, which would make Tf negative and the filter unstable.

A default Tf that is not representable as a positive finite `Scalar` (kd/(kp*N) overflowing, or underflowing to zero) is refused with `unrepresentable_gain`. For a D-only or I-D channel, state the time constant directly:

```cpp
using Pid = ctrlpp::pid<double, 1, ctrlpp::deriv_filter>;
using Vec = Pid::vector_t;

Pid::config_type cfg{};
cfg.ki = Vec::Constant(0.5);
cfg.kd = Vec::Constant(0.1);
cfg.template policy<ctrlpp::deriv_filter>().tf = {0.02};  // seconds; kp is zero, so there is no default
```

## Behavior

The ideal derivative s*Td is filtered by a first-order system with time constant Tf = Td/N, so the derivative term becomes Kp*Td*s / (1 + s*Td/N), or kd*s / (1 + s*Tf) in parallel gains. Its gain at high frequencies is limited to kd/Tf, which is Kp*N with the default time constant, rather than growing without bound (Astrom and Hagglund, *Advanced PID Control*, 2006). In discrete time the filter is the backward-difference form D_k = alpha*D_(k-1) + (1 - alpha)*D_raw with alpha = Tf/(Tf + dt), which lies between 0 and 1 for every positive Tf and dt, so the filter is stable for all accepted parameters.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'output', '' using 1:3 with lines title 'control'"

#include <ctrlpp/control/pid.h>

#include <cmath>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, ctrlpp::deriv_filter>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(0.5);
    cfg.kd = Vec::Constant(0.1);
    cfg.template policy<ctrlpp::deriv_filter>().n = {10.0};

    auto created = Pid::create(cfg);
    if (!created.has_value()) return 1;
    auto& ctrl = *created;

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 5.0; t += dt) {
        // Simulate noisy measurement
        double noise = 0.05 * std::sin(100.0 * t);
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y + noise);
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
- [setpoint-filter](setpoint-filter.md)<br/>reference signal filtering
- [guides/pid/composition](../../../guides/pid/composition.md)<br/>composing policies
- [background/pid](../../../background/pid.md)<br/>derivative filtering theory
