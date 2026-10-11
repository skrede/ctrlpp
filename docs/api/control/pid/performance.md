# perf_assessment

Online performance assessment policy that computes control loop quality metrics during operation. Four metric types are available as compile-time template parameters: IAE (Integral of Absolute Error), ISE (Integral of Squared Error), ITAE (Integral of Time-weighted Absolute Error), and oscillation_detect (Hagglund's load-disturbance oscillation detector [1, 2]). Multiple metrics can be active simultaneously.

## Header

| Form | Header |
|------|--------|
| `ctrlpp::perf_assessment<Metrics...>` | `#include <ctrlpp/control/pid_policies.h>` |

## Template Variants

```cpp
ctrlpp::perf_assessment<ctrlpp::IAE>                          // integral of absolute error
ctrlpp::perf_assessment<ctrlpp::ISE>                          // integral of squared error
ctrlpp::perf_assessment<ctrlpp::ITAE>                         // integral of time-weighted absolute error
ctrlpp::perf_assessment<ctrlpp::oscillation_detect>           // load-disturbance oscillation detector
ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE>             // multiple metrics
ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::oscillation_detect>  // mixed
```

## Config Fields

### perf_assessment

`perf_assessment<Metrics...>::config<Scalar, NY>` carries the oscillation detector's parameters. IAE, ISE, and ITAE read none of them, and the fields are checked only when `oscillation_detect` is among the metrics. Set them through the policy config, for example `cfg.template policy<PA>().amplitude = {0.01};`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `amplitude` | `std::array<Scalar, NY>` | none (required) | `a`, in error units: the lower limit of an oscillation amplitude worth reporting, at the ultimate frequency. Hagglund suggests 1% of the signal range [2, Sec. 2.1]. Zero is the unset value and is refused. |
| `ultimate_frequency` | `std::array<Scalar, NY>` | none (required) | `omega_u`, in rad/s: the loop's ultimate frequency, for example from a relay autotuning experiment. Zero is the unset value and is refused. |
| `load_count_limit` | `Scalar` | `10` | `n_lim`: the weighted load count above which a channel is oscillating [1], as restated in [2, Sec. 2.1]. |
| `supervision_periods` | `Scalar` | `50` | The supervision window `T_sup` in ultimate periods `T_u = 2*pi/omega_u`. The published window is `T_sup = 5 * n_lim * T_u` [1], [2, Sec. 2.1], so the default 50 belongs to the default `n_lim = 10`: changing `load_count_limit` means changing this field too if the published coupling is to be kept. |

`create()` and `set_params()` refuse an `amplitude` or `ultimate_frequency` component that is not positive and finite with `pid_config_error::missing_oscillation_parameters`, and a `load_count_limit` or `supervision_periods` that is not positive and finite with `pid_config_error::invalid_oscillation_supervision`.

## Metric Tags

```cpp
struct IAE {};                 // sum of |e| * dt
struct ISE {};                 // sum of e^2 * dt
struct ITAE {};                // sum of t * |e| * dt
struct oscillation_detect {};  // weighted count of load disturbances
```

## Behavior

The performance assessment policy accumulates metrics passively during each `compute()` call. It does not modify the PID output &mdash; it is read-only observation. The accumulated metrics can be queried at any time for logging, tuning evaluation, or online diagnostics.

### Oscillation detection

Each channel runs Hagglund's recursive detector [2, Sec. 2.1] on the control error `e` with the cycle's step `h`:

1. The integrated absolute error since the last sign change of `e`, IAE, grows by `|e| * h` every cycle.
2. At a sign change, the interval just closed is a load disturbance when its IAE exceeds `IAE_lim = 2 * a / omega_u`; the IAE then restarts from the current sample. An error of exactly zero is not a sign change and adds nothing.
3. The load events are counted with exponential forgetting, `x = gamma * x + load`, where `gamma = 1 - h / T_sup` [2, Eq. (2)] and `T_sup = supervision_periods * 2 * pi / omega_u`. Eq. (2) assumes `h` much shorter than `T_sup`; a step longer than the whole window sets `gamma` to zero rather than negative.

`metric<oscillation_detect>()` returns `x` per channel, the rate of load detections plotted in [2, Fig. 2]. `oscillating()` returns true while `x > load_count_limit` in any channel. It is a query of the current state: nothing resets `x` on detection, so the verdict clears by itself as `x` decays once the oscillation stops. Hagglund's implemented program instead latches an alarm and resets `x` [2, Secs. 2.1 and 5]; latching and acknowledgement are left to the application.

With the defaults and a load at every half period, `gamma` over half a period is `exp(-1/100)`, about 0.990, whatever `T_u` is. The weighted count after ten loads is about 9.56 and after eleven about 10.47, so `oscillating()` first returns true at the eleventh load event, about 5.5 ultimate periods after the oscillation starts. A sustained oscillation drives `x` toward about 100.

`x` rises in a step of 1 at each load event and decays smoothly between them. When a sustained oscillation holds `x` near `load_count_limit`, `x` saw-tooths across the threshold and `oscillating()` toggles at roughly the half-period cadence. A caller that needs a steady signal should apply hysteresis to `metric<oscillation_detect>()` rather than use the boolean.

The detector assumes the controller has integral action, so that the average control error is zero [2, Sec. 2.1]. Data from manual control, from periods when the controller tracks an external signal, and from periods when the signals are saturated should be excluded [2, Sec. 3.1]; call `reset_metrics()` when such a period ends, which discards what was accumulated during it (the other metrics included).

In `float`, `gamma = 1 - h / T_sup` rounds to exactly 1 when `h / T_sup` is below about 3e-8 (for example a 1 ms step with a window longer than about nine hours), and `x` then stops decaying. Use `double` for such loops.

## Usage Example

```cpp
#include <ctrlpp/control/pid.h>

#include <iomanip>
#include <iostream>
#include <numbers>

int main()
{
    using PA = ctrlpp::perf_assessment<ctrlpp::IAE, ctrlpp::ISE, ctrlpp::oscillation_detect>;
    using Pid = ctrlpp::pid<double, 1, PA>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(3.0);
    cfg.ki = Vec::Constant(1.0);
    // 1% of a unit signal range, and the loop's ultimate frequency from a relay experiment.
    cfg.policy<PA>().amplitude = {0.01};
    cfg.policy<PA>().ultimate_frequency = {2.0 * std::numbers::pi / 4.0};

    auto created = Pid::create(cfg);
    if (!created.has_value()) return 1;
    auto& ctrl = *created;

    double y = 0.0;
    constexpr double dt = 0.01;

    for (double t = 0.0; t < 10.0; t += dt) {
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y);
        // A refused cycle produced no command; the caller decides what the
        // actuator does. This sample stops.
        auto step = ctrl.compute(sp, meas, dt);
        if (!step.has_value()) return 1;
        const auto& u = *step;
        y = 0.9 * y + 0.1 * u[0];
    }

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Final output: " << y << "\n";
    std::cout << "IAE: " << ctrl.metric<ctrlpp::IAE>()[0] << "\n";
    std::cout << "Load detection rate: " << ctrl.metric<ctrlpp::oscillation_detect>()[0] << "\n";
    std::cout << "Oscillating: " << std::boolalpha << ctrl.oscillating() << "\n";
}
```

## See Also

- [PID overview](README.md)<br/> parent PID documentation
- [anti-windup](anti-windup.md)<br/> performance depends on proper windup handling
- [derivative-filter](derivative-filter.md)<br/>noise filtering affects performance metrics
- [guides/pid/composition](../../../guides/pid/composition.md)<br/>composing policies
- [background/pid](../../../background/pid.md)<br/>performance index definitions

## References

[1] T. Hagglund, "A Control-Loop Performance Monitor," Control Engineering Practice, vol. 3, no. 11, pp. 1543-1551, 1995.

[2] T. Hagglund, "Industrial Applications of Automatic Performance Monitoring Tools," Proc. 15th IFAC World Congress, Barcelona, 2002.
