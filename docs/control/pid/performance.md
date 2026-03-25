# perf_assessment

Online performance assessment policy that computes control loop quality metrics during operation. Four metric types are available as compile-time template parameters: IAE (Integral of Absolute Error), ISE (Integral of Squared Error), ITAE (Integral of Time-weighted Absolute Error), and oscillation_detect (zero-crossing rate monitoring). Multiple metrics can be active simultaneously.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::perf_assessment<Metrics...>` | `#include <ctrlpp/control/pid_policies.h>` |

## Template Variants

```cpp
ctrlpp::perf_assessment<ctrlpp::IAE>                          // integral of absolute error
ctrlpp::perf_assessment<ctrlpp::ISE>                          // integral of squared error
ctrlpp::perf_assessment<ctrlpp::ITAE>                         // integral of time-weighted absolute error
ctrlpp::perf_assessment<ctrlpp::oscillation_detect>           // zero-crossing rate
ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE>             // multiple metrics
ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::oscillation_detect>  // mixed
```

## Config Fields

### perf_assessment (all metric variants)

No additional config fields for IAE, ISE, or ITAE.

### oscillation_detect

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `crossing_rate_threshold` | `double` | `5.0` | Zero-crossing rate (crossings/second) above which the controller is considered oscillating. |

## Metric Tags

```cpp
struct IAE {};                 // sum of |e| * dt
struct ISE {};                 // sum of e^2 * dt
struct ITAE {};                // sum of t * |e| * dt
struct oscillation_detect {};  // counts error sign changes per unit time
```

## Behavior

The performance assessment policy accumulates metrics passively during each `compute()` call. It does not modify the PID output -- it is read-only observation. The accumulated metrics can be queried at any time for logging, tuning evaluation, or online diagnostics.

## Usage Example

```cpp
#include <ctrlpp/control/pid.h>

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1,
        ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE>>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(3.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.1);

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 10.0; t += dt) {
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y);
        auto u = ctrl.compute(sp, meas, dt);
        y = 0.9 * y + 0.1 * u[0];
    }

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Simulation complete after 10s\n";
    std::cout << "Final output: " << y << "\n";
}
```

## See Also

- [PID overview](README.md) -- parent PID documentation
- [anti-windup](anti-windup.md) -- performance depends on proper windup handling
- [derivative-filter](derivative-filter.md) -- noise filtering affects performance metrics
- [guides/pid/composition](../../guides/pid/composition.md) -- composing policies
- [reference/pid-theory](../../reference/pid-theory.md) -- performance index definitions
