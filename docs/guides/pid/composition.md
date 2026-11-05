# Policy-Based PID Composition

ctrlpp's PID controller uses compile-time policy composition. Each policy
adds a specific behavior &mdash; anti-windup, derivative filtering, rate limiting
-- as a template parameter. The code of a policy you do not compose is
compiled out.

## How It Works

The `pid` class template accepts a variadic policy pack:

```cpp
ctrlpp::pid<Scalar, NY, Policies...>
```

Each policy type in `Policies...` that declares a configuration struct adds
it to the `policies` tuple of the controller's `config_type`, and
`cfg.template policy<Policy>()` returns a reference to it. The base fields
(`kp`, `ki`, `kd`, the limits and the setpoint weights) are the same for every
composition.

## Composing Three Policies

This example composes `anti_windup<back_calc>`, `deriv_filter`, and
`rate_limit`. Together they give:

- **Anti-windup**: prevents integrator saturation using back-calculation
- **Derivative filtering**: low-pass filters the derivative term to reduce
  noise amplification
- **Rate limiting**: bounds the rate of change of the control signal

```cpp
// Usage: ./pid_composition | gnuplot -p -e "set datafile separator ','; plot '-' skip 1 using 1:3 with lines title 'measurement', '' using 1:4 with lines title 'control'"
#include <ctrlpp/control/pid.h>

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1,
                             ctrlpp::anti_windup<ctrlpp::back_calc>,
                             ctrlpp::deriv_filter,
                             ctrlpp::rate_limit>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(3.0);
    cfg.ki = Vec::Constant(1.5);
    cfg.kd = Vec::Constant(0.5);
    cfg.output_min = Vec::Constant(-5.0);
    cfg.output_max = Vec::Constant(5.0);

    // Anti-windup: back-calculation gain in 1/time (left at zero it auto-defaults to
    // 1/sqrt(Ti*Td), or 1/Ti without derivative action)
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb = {1.0};

    // Derivative filter: N sets the bandwidth (higher = less filtering)
    cfg.template policy<ctrlpp::deriv_filter>().n = {10.0};

    // Rate limit: maximum change per second
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {2.0};

    auto created = Pid::create(cfg);
    if (!created.has_value()) return 1;
    auto& ctrl = *created;

    // Simulated first-order plant
    double y = 0.0;
    constexpr double setpoint = 1.0;
    constexpr double dt = 0.01;
    constexpr double a = 0.9;
    constexpr double duration = 10.0;

    std::cout << "time,setpoint,measurement,control\n";

    for (double t = 0.0; t < duration; t += dt)
    {
        double disturbance = (t >= 5.0) ? -0.5 : 0.0;

        auto sp = Vec::Constant(setpoint);
        auto meas = Vec::Constant(y + disturbance);
        // A refused cycle produced no command; the caller decides what the
        // actuator does. This sample stops.
        auto step = ctrl.compute(sp, meas, dt);
        if (!step.has_value()) return 1;
        const auto& u = *step;

        y = a * y + (1.0 - a) * u[0];

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << setpoint << "," << y << ","
                  << u[0] << "\n";
    }
}
```

## The Config Struct

When you add policies, the config struct grows. For the three-policy
composition above, `cfg` includes:

| Source                    | Field        | Purpose                              |
| ------------------------ | ------------ | ------------------------------------ |
| Base PID                 | `kp`         | Proportional gain                    |
| Base PID                 | `ki`         | Integral gain                        |
| Base PID                 | `kd`         | Derivative gain                      |
| Base PID                 | `output_min` | Lower output clamp                   |
| Base PID                 | `output_max` | Upper output clamp                   |
| Base PID                 | `derivative_on_error` | Differentiate `c*r - y` instead of `-y` |
| Base PID                 | `b`, `c`     | Proportional and derivative setpoint weights |
| `anti_windup<back_calc>` | `kb`         | Back-calculation gain                |
| `deriv_filter`           | `n`          | Derivative filter divisor N (default Tf = kd/(kp*N)) |
| `deriv_filter`           | `tf`         | Explicit derivative filter time      |
| `rate_limit`             | `rate_max`   | Maximum output rate of change        |

Access policy fields through `cfg.template policy<PolicyType>()`. Composing
`anti_windup<conditional_integration>` instead would add its
`error_threshold`, which must be set positive in every channel.

## Available Policies

| Policy                                    | Purpose                                  |
| ----------------------------------------- | ---------------------------------------- |
| `anti_windup<back_calc>`                  | Back-calculation anti-windup             |
| `anti_windup<clamping>`                   | Clamping anti-windup                     |
| `anti_windup<conditional_integration>`    | Integration suspended above an error threshold |
| `deriv_filter`                            | Low-pass filter on derivative term       |
| `setpoint_filter`                         | Setpoint pre-filter                      |
| `pv_filter`                               | Process variable filter                  |
| `feed_forward<Callable>`                  | Feed-forward from a callable of setpoint and step |
| `rate_limit`                              | Rate-of-change output limiter            |
| `velocity_form`                           | Incremental (velocity) PID form          |
| `isa_form`                                | ISA standard form (Ti, Td parameters)    |
| `perf_assessment<Metrics...>`             | Online performance metrics (IAE, ISE, ITAE, oscillation detection) |
| `forward_euler`, `tustin`                 | Integral discretization (default backward Euler) |

`velocity_form` does not compose with `anti_windup`: the combination is a compile error, because clamping
the accumulated output already prevents windup.

## Compile-Time Composition

Policies are resolved entirely at compile time with `if constexpr`: adding
`deriv_filter` adds only the filter arithmetic, with no virtual dispatch and
no run-time branch on policy presence. Every composition still carries the
same set of state members, checks its inputs and the finiteness of its state,
and snapshots its carried state at the start of each cycle so that a refused
cycle can be rolled back.

## Next Steps

- [Anti-Windup](../../api/control/pid/anti-windup.md)<br/> detailed anti-windup API
- [Derivative Filter](../../api/control/pid/derivative-filter.md)<br/> filter bandwidth selection
- [Rate Limit](../../api/control/pid/rate-limit.md)<br/> rate limiter API
- [PID API Reference](../../api/control/pid/README.md)<br/> full method signatures
- [PID Theory](../../background/pid.md)<br/> mathematical background
