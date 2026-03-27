# Policy-Based PID Composition

ctrlpp's PID controller uses compile-time policy composition. Each policy
adds a specific behaviour -- anti-windup, derivative filtering, rate limiting
-- as a template parameter. You pay zero runtime cost for policies you do not
use.

## How It Works

The `pid` class template accepts a variadic policy pack:

```cpp
ctrlpp::pid<Scalar, NX, NU, NY, Policies...>
```

Each policy type in `Policies...` extends the controller's `config_type` with
policy-specific fields. When you compose multiple policies, the config struct
merges all their fields automatically.

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
    using Pid = ctrlpp::pid<double, 1, 1, 1,
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

    // Anti-windup: back-calculation gain (typically 1/Ti)
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb = {1.0};

    // Derivative filter: N sets the bandwidth (higher = less filtering)
    cfg.template policy<ctrlpp::deriv_filter>().n = {10.0};

    // Rate limit: maximum change per second
    cfg.template policy<ctrlpp::rate_limit>().rate_max = {2.0};

    Pid ctrl(cfg);

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
        auto u = ctrl.compute(sp, meas, dt);

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
| `anti_windup<back_calc>` | `kb`         | Back-calculation gain                |
| `deriv_filter`           | `n`          | Derivative filter bandwidth          |
| `rate_limit`             | `rate_max`   | Maximum output rate of change        |

Access policy fields through `cfg.template policy<PolicyType>()`.

## Available Policies

| Policy                                    | Purpose                                  |
| ----------------------------------------- | ---------------------------------------- |
| `anti_windup<back_calc>`                  | Back-calculation anti-windup             |
| `anti_windup<clamping>`                   | Clamping (conditional) anti-windup       |
| `anti_windup<conditional_integration>`    | Conditional integration anti-windup      |
| `deriv_filter`                            | Low-pass filter on derivative term       |
| `setpoint_filter`                         | Setpoint pre-filter                      |
| `pv_filter`                               | Process variable filter                  |
| `feed_forward<Callable>`                  | Feed-forward from setpoint               |
| `rate_limit`                              | Rate-of-change output limiter            |
| `velocity_form`                           | Incremental (velocity) PID form          |
| `isa_form`                                | ISA standard form (Ti, Td parameters)    |
| `perf_assessment<Metrics...>`             | Online performance metrics (IAE, ISE, ITAE) |

## Zero-Cost Abstraction

Policies are resolved entirely at compile time. A `pid<double, 1, 1, 1>` with
no policies compiles to a minimal P+I+D computation. Adding `deriv_filter`
adds only the filter arithmetic -- no virtual dispatch, no branch on policy
presence.

The optimiser sees through the template instantiation and inlines everything
into a single function body.

## Next Steps

- [Anti-Windup](../../api/control/pid/anti-windup.md) -- detailed anti-windup API
- [Derivative Filter](../../api/control/pid/derivative-filter.md) -- filter bandwidth selection
- [Rate Limit](../../api/control/pid/rate-limit.md) -- rate limiter API
- [PID API Reference](../../api/control/pid/README.md) -- full method signatures
- [PID Theory](../../background/pid.md) -- mathematical background
