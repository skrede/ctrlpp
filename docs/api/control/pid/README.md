# pid

A policy-based PID controller with compile-time feature composition. The template
signature `pid<Scalar, NY, Policies...>` lets you opt into exactly the
features you need &mdash; anti-windup, derivative filtering, setpoint weighting,
velocity form, ISA form, feed-forward, rate limiting, and performance assessment --
without paying for what you don't use.

The base PID provides proportional, integral, and derivative action with output
clamping and setpoint weighting (b, c parameters). Each policy extends the
controller with one additional behavior, composed at compile time via variadic
template packs.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::pid<Scalar, NY, Policies...>` | `#include <ctrlpp/control/pid.h>` |
| (convenience) | `#include <ctrlpp/pid.h>` |
| `ctrlpp::pid_config<Scalar, NY, Policies...>` | `#include <ctrlpp/control/pid_config.h>` (convenience `<ctrlpp/pid_config.h>`; included by `pid.h`) |
| the policy tags | `#include <ctrlpp/control/pid_policies.h>` (convenience `<ctrlpp/pid_policies.h>`; included by `pid.h`) |
| `pid_config_error`, `pid_step_error`, `pid_health` | `#include <ctrlpp/control/pid_error.h>` (included by `pid.h`) |
| `ctrlpp::pid_performance_tracker` | `#include <ctrlpp/control/pid_performance.h>` (included by `pid.h`; see [performance](performance.md)) |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | A floating-point type (`float`, `double`, `long double`). The constraint sits on `pid_config`, so `pid` with any other type does not compile. |
| `NY` | `std::size_t`, positive | Number of channels &mdash; the size of the gain, setpoint, measurement and output vectors. `pid_config` refuses `NY = 0` with a `static_assert`. |
| `Policies...` | policy types | Zero or more policy tags from `pid_policies.h`, or the integration tags `forward_euler` and `tustin` (see [Integral discretization](#integral-discretization)) |

## Type Aliases

```cpp
using config_type = pid_config<Scalar, NY, Policies...>;
using vector_t    = Vector<Scalar, NY>;
```

## Config (`pid_config`)

The config struct holds all tuning parameters. Policy-specific fields are accessed
via `cfg.template policy<PolicyType>()`.

### Base Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `kp` | `Vector<Scalar, NY>` | Zero | Proportional gain per channel |
| `ki` | `Vector<Scalar, NY>` | Zero | Integral gain per channel (or Ti when `isa_form` active) |
| `kd` | `Vector<Scalar, NY>` | Zero | Derivative gain per channel (or Td when `isa_form` active) |
| `output_min` | `Vector<Scalar, NY>` | `lowest()` | Lower output clamp per channel |
| `output_max` | `Vector<Scalar, NY>` | `max()` | Upper output clamp per channel |
| `derivative_on_error` | `bool` | `false` | If true, the D term differentiates the weighted error `c*r - y`; if false, it differentiates `-y` only, so a setpoint step causes no derivative kick and `c` is not read |
| `b` | `Vector<Scalar, NY>` | `1.0` | Proportional setpoint weight: P term uses `b*r - y` |
| `c` | `Vector<Scalar, NY>` | `1.0` | Derivative setpoint weight: D term uses `c*r - y` when `derivative_on_error` is true |
| `policies` | `std::tuple` of the policy configs | each value-initialized | One configuration struct per composed policy that declares one; read and written through `policy<P>()` |

`r` and `y` are the setpoint and measurement after the `setpoint_filter` and `pv_filter` policies, when those are composed. The derivative term is zero on the first cycle after `create` or `reset`, since there is no previous sample to difference.

### Policy Config Access

```cpp
cfg.template policy<anti_windup<back_calc>>().kb = {1.0};
cfg.template policy<deriv_filter>().n = {10.0};
cfg.template policy<setpoint_filter>().tf = {0.1};
cfg.template policy<pv_filter>().tf = {0.05};
cfg.template policy<rate_limit>().rate_max = {0.5};
cfg.template policy<feed_forward<FF>>().feed_forward_cb = my_feed_forward;
```

## Available Policies

| Policy | Description | Config fields | Page |
|--------|-------------|---------------|------|
| `anti_windup<back_calc>` | Back-calculation anti-windup | `kb` | [anti-windup](anti-windup.md) |
| `anti_windup<clamping>` | Integrator clamping | (none) | [anti-windup](anti-windup.md) |
| `anti_windup<conditional_integration>` | Error-threshold gated integration | `error_threshold` | [anti-windup](anti-windup.md) |
| `deriv_filter` | Low-pass filtered derivative | `n`, `tf` | [derivative-filter](derivative-filter.md) |
| `setpoint_filter` | Setpoint prefilter | `tf` | [setpoint-filter](setpoint-filter.md) |
| `pv_filter` | Process variable filter | `tf` | [setpoint-filter](setpoint-filter.md) |
| `velocity_form` | Incremental output form; does not compose with `anti_windup` | (none) | [velocity-form](velocity-form.md) |
| `isa_form` | ISA standard form (Kp, Ti, Td) | (none) | [isa-form](isa-form.md) |
| `feed_forward<Callable>` | Additive feed-forward | `feed_forward_cb` | [feed-forward](feed-forward.md) |
| `rate_limit` | Output rate limiting | `rate_max` | [rate-limit](rate-limit.md) |
| `perf_assessment<Metrics...>` | Online performance metrics | `amplitude`, `ultimate_frequency`, `load_count_limit`, `supervision_periods` (read only by `oscillation_detect`) | [performance](performance.md) |

## Integral discretization

The integral increment per cycle is chosen by an optional tag from `ctrlpp/model/discretize.h` (included by `pid_policies.h`) in the policy pack, with `e` the error after the input filters and `ki` the internal parallel integral gain:

| Tag | Increment |
|-----|-----------|
| none (default) | `ki * e[k] * dt`, the backward Euler rule |
| `ctrlpp::forward_euler` | `ki * e[k-1] * dt`, so the first cycle after `create` or `reset` adds nothing |
| `ctrlpp::tustin` | `ki * (e[k] + e[k-1]) / 2 * dt`, the trapezoidal rule |

The other tags of that header have no effect on `pid`. The increment is the same in velocity form.

## create

```cpp
static expected<pid, pid_config_error> create(const config_type& cfg);
```

Validates the configuration and constructs the controller, or reports the first fault it finds. A configuration that cannot work is refused where it is given rather than discovered on the first cycle. The checks run in the order of the table, so a configuration with several faults reports the earliest; a check on a policy field runs only when that policy is composed. A condition on a derived gain reads the internal parallel gains: in ISA form `kp/Ti` and `kp*Td`. `create` is the only construction path: the constructor is private, so no configuration reaches the controller without these checks. The table lists every enumerator of `pid_config_error` in declaration order; `create` and `set_params` test the first fourteen, and the last belongs to `set_integral`.

```cpp
auto made = ctrlpp::pid<double, 1>::create(cfg);
if(!made) {
    // made.error() names the first field that cannot work
}
```

| Enumerator | Cause | Repair |
|---|---|---|
| `non_finite_gain` | A `kp`, `ki` or `kd` component is `NaN` or infinite. In ISA form `ki` is the integral time, and a `NaN` or negative-infinite one is refused; `+inf` is accepted and means no integral action. | The tuning that produced the gains. |
| `zero_integral_time` | In ISA form an integral time `ki` is zero, which is an infinite integral gain Kp / Ti, not "off". | A nonzero integral time, or for no integral action an integral time of `+inf` or the parallel form with `ki = 0`. |
| `invalid_output_limits` | An output limit is `NaN`, `output_min` exceeds `output_max` in some channel, `output_min` is `+inf`, or `output_max` is `-inf`. Equal limits are accepted. | The actuator range. |
| `non_finite_setpoint_weight` | A `b` or `c` component is `NaN` or infinite. | The setpoint weights. |
| `invalid_filter_time_constant` | A `setpoint_filter` or `pv_filter` time constant is `NaN`, negative or `+inf`. Zero means no filter. | The filter's `tf`. |
| `invalid_rate_limit` | A `rate_limit` `rate_max` is `NaN` or negative. Zero and `+inf` mean no limit. | The rate limit. |
| `invalid_back_calculation_gain` | An explicit `anti_windup<back_calc>` gain `kb` is `NaN`, infinite or negative. Zero selects the default gain. | A positive tracking gain. |
| `invalid_error_threshold` | An `anti_windup<conditional_integration>` `error_threshold` is not positive: `NaN`, zero (the default) or negative. | A positive integration threshold in every channel. |
| `invalid_derivative_filter` | A `deriv_filter` `n` is `NaN` or infinite, or an explicit `tf` is `NaN`, infinite, zero or negative. | The filter divisor or time constant. |
| `missing_oscillation_parameters` | `perf_assessment` includes `oscillation_detect` and an `amplitude` or `ultimate_frequency` component is not positive and finite. Both default to zero, the unset value. | The acceptable oscillation amplitude and the loop's ultimate frequency. |
| `invalid_oscillation_supervision` | `perf_assessment` includes `oscillation_detect` and `load_count_limit` or `supervision_periods` is not positive and finite. | The detector's load count limit and window. |
| `undefined_tracking_time` | `anti_windup<back_calc>` leaves a channel's `kb` at zero, the channel has integral action (`ki != 0`), and either it has neither proportional nor derivative action or its nonzero gains do not share a sign, so the default tracking time does not exist. A channel with `ki = 0` keeps `kb = 0`, the limit of an infinite tracking time. | An explicit positive `kb`, or gains of one sign. |
| `undefined_derivative_filter` | `deriv_filter` is composed, a channel has derivative action (`kd != 0`) and no `tf`, and its default Tf = kd/(kp*n) is not a positive time: `kp` is zero, `n` is zero or negative, or `kd` and `kp` differ in sign. | An explicit `tf`, or a positive `n`. |
| `unrepresentable_gain` | The fields are well posed, but a quantity derived from them is not representable in `Scalar`: in ISA form `kp` over the integral time, or `kp` times the derivative time, overflows, or rounds to zero while both factors are nonzero (an integral time of `+inf` excepted); the default derivative filter time constant kd/(kp*n) is not a positive finite number; or the default back-calculation gain is not finite, or rounds to zero on a channel with integral action. | Gains of less extreme magnitude, an explicit `tf` or `kb`, or a wider `Scalar`. |
| `non_finite_integral` | Returned only by `set_integral`: a component of the value is `NaN` or infinite. | A finite integrator value. |

## Methods

### compute

```cpp
expected<vector_t, pid_step_error> compute(const vector_t& sp, const vector_t& meas, Scalar dt);
```

Computes the control output given a setpoint, measurement, and time step. On success the result carries the clamped output vector, or in velocity form the increment.

The first four causes in the table are tested before any member is written. A cycle whose result turns out non-finite is rolled back from a snapshot taken when it began. Either way a rejected cycle leaves every piece of carried state -- the integrator, the error and derivative histories, the input filter states, the accumulated output, the performance metrics -- bitwise unchanged, and the caller may retry on the next sample.

**A refused cycle produced no command.** That is a different situation from an estimator declining a measurement, where the estimate simply stands: an actuator will be driven by something regardless of what `compute` returns, and the caller must choose what. The three defensible responses are to hold the command the last successful cycle returned, to command a configured safe value, or to fail over to a redundant channel. The controller does not choose, because the right answer is a property of the plant: holding is correct for a slow thermal loop and dangerous for an unstable attitude loop.

A non-positive step is refused with `invalid_timestep` rather than answered with the previous output, which the caller could not tell from a freshly computed command: a stopped clock would read as a steady loop.

The table lists every enumerator of `pid_step_error` in declaration order. A cycle tests the first four in that order, the most upstream cause first.

| Enumerator | Cause | Where the caller repairs it |
|---|---|---|
| `non_finite_state` | Some piece of carried state is already non-finite when the cycle begins. This refusal latches `health()`. | `reset`, or `set_integral` with a finite value for the integrator alone. Reported first because it is upstream of every argument. |
| `invalid_timestep` | The step is not a positive finite duration. A `NaN` step cannot be caught by `dt <= 0` alone, because every comparison against a `NaN` is false. | The caller's clock. |
| `non_finite_setpoint` | The commanded setpoint has a non-finite component. | The trajectory, outer loop, or operator input that generates the reference. |
| `non_finite_measurement` | The process variable has a non-finite component. | The sensor or estimator. Separate from the setpoint even though the error is their difference, because the two arrive from different subsystems. |
| `non_finite_tracking_signal` | The tracking signal handed to the four-argument overload is non-finite. That overload tests it before it runs the cycle, so it is reported ahead of the four causes above. | The signal being tracked. |
| `non_finite_result` | The cycle ran on finite inputs, but its output, a piece of the state it would carry into the next cycle, or a performance metric came out non-finite: for example an integrator that overflows `Scalar`, or a `NaN` from an infinite term meeting its opposite. A proportional, derivative or feed-forward term that overflows to infinity is clamped to the output limit like any other command, so it fails the cycle only where it reaches carried state. In the tracking overload it is also returned when the re-seeded integrator would be non-finite. The cycle is rolled back. | Gains, limits or signals of less extreme magnitude, or a wider `Scalar`. |

### compute (with tracking signal)

```cpp
expected<vector_t, pid_step_error> compute(const vector_t& sp, const vector_t& meas, Scalar dt,
                                           const vector_t& tracking_signal);
```

Overload for controller output tracking (bumpless transfer in cascade configurations). After the cycle succeeds the controller is re-seeded so its output tracks an external signal, and the returned value is that cycle's own computed command, or increment in velocity form. Position form re-seeds the integrator; velocity form re-seeds the accumulated output.

In position form, after the cycle succeeds the integrator is assigned `I = clamp(w, output_min, output_max) - (P + D + FF)`, where `w` is the tracking signal and `P + D + FF` is the sum of the proportional, derivative and feed-forward terms of the cycle that formed the output, taken before rate limiting and saturation. It is not reconstructed from the returned command, which may have been saturated or rate-limited, and the anti-windup step of that cycle has no influence on it. A target outside the output limits is clamped first, so the integrator never holds an offset that drives the unconstrained command past a limit. For example, with `kp = 1`, `ki = 0.1`, limits of `+-10`, setpoint 20 and measurement 0, tracking 5 leaves the integrator at -15, and tracking 50 leaves it at `10 - 20 = -10`.

The position form re-seeds two more pieces of carried state, so the first automatic cycle after tracking starts from the command the actuator actually received:

- With `rate_limit`, the limiter's anchor becomes `clamp(w, output_min, output_max)`, not the discarded automatic command. For example, with `kp = 1`, `ki = 0.5`, `rate_max = 1`, `dt = 1/64`, setpoint 1 and measurement 0, four automatic cycles ramp the output to `4/64`; a tracking cycle at `w = 0.5` returns `5/64`, and the next automatic cycle emits `0.5 + 1/128` instead of ramping on from `5/64` to `6/64`.
- With `anti_windup<clamping>`, the carried constraint verdict is cleared, so the next automatic cycle does not withhold its increment on the constraint of the command that tracking replaced.

Tracking assigns the integrator even while `freeze_integral` holds it. Like `set_integral` it is an explicit assignment, not anti-windup, and the freeze stays in force afterward: the following automatic cycles leave the tracked integrator unchanged until it is unfrozen.

In velocity form, after the cycle succeeds the accumulated output is assigned `clamp(w, output_min, output_max)`, and the next cycle's increment is accumulated from it. For example, with `kp = 1` and limits of `+-10`, a cycle at setpoint 1 tracking 9.5 re-seeds the accumulated output to 9.5, and a following setpoint step to 3 emits 0.5 and reports `saturated()`.

The tracking signal is checked *before* the cycle is delegated, because it is written straight into carried state once the cycle succeeds: admitting a non-finite one would destroy that state after the cycle had already committed the rest. A rejection from the delegated cycle is forwarded unchanged with its own specific cause, and the tracking assignment does not run.

### health

```cpp
pid_health health() const;
```

Reports whether the controller has found itself carrying non-finite state -- a question a per-cycle result cannot answer, because it outlives the call. It is `pid_health::ok` after `create` and `reset`.

It changes to `non_finite_carried_state` only when a cycle is refused with `non_finite_state`, and then latches until `reset`; neither `set_params` nor `set_integral` clears it. Every other refusal, `non_finite_result` included, leaves it as it was, since those cycles leave no state behind. No known public route leads to `non_finite_carried_state`: `create`, `set_params` and `set_integral` refuse non-finite values, and a cycle whose result is non-finite is rolled back. The latch is kept as defense in depth.

### set_params

```cpp
expected<void, pid_config_error> set_params(const config_type& new_cfg);
```

Refuses the configuration with the same checks, in the same order, as `create`; a refusal leaves the current configuration, the integrator and every other piece of state bitwise unchanged. On success it updates tuning parameters with bumpless transfer. The integral state is stored in output units and enters the output directly, so its contribution is already continuous across a gain change; `set_params` leaves the integral state untouched, which is what makes gain scheduling bumpless. No rescaling is applied (rescaling would step the integral contribution instead of preserving it). The filter states, the histories and the performance metrics are kept as well. `set_params` copy-assigns the configuration, so it compiles only when `config_type` is copy-assignable; see [feed-forward](feed-forward.md#limitations) for the one policy where that can fail.

### error

```cpp
const vector_t& error() const;
```

Returns the most recent error vector (setpoint minus measurement, after filtering).

### integral

```cpp
const vector_t& integral() const;
```

Returns the current integral accumulator state, in output units. In velocity form there is no integrator and it stays zero.

### params

```cpp
const config_type& params() const;
```

Returns a const reference to the current config.

### saturated

```cpp
bool saturated() const;
```

Returns true if the most recent successful `compute()` call clamped the output to `output_min` or `output_max` in some channel; rate limiting alone does not set it. In velocity form it is true when the accumulated output was clamped to the output limits on that call. A refused cycle leaves it as it was.

### reset

```cpp
void reset();
```

Resets all internal state (integral, previous error, setpoint and measurement, derivative and filter states, accumulated output, performance metrics) to zero, so the next cycle is treated as the first. It also releases `freeze_integral`, clears `saturated()` and returns `health()` to `ok`. The configuration is kept.

### freeze_integral

```cpp
void freeze_integral(bool freeze = true);
```

Freezes or unfreezes the integral accumulator. While frozen, the integral term does not change, including through back-calculation anti-windup, which resumes once the accumulator is unfrozen. The freeze holds the integrator against the integral increment and anti-windup only: the explicit assignments, `set_integral` and the tracking overload of `compute`, still write it.

### set_integral

```cpp
expected<void, pid_config_error> set_integral(const vector_t& val);
```

Directly sets the integral accumulator to a given value. A `NaN` or infinite component is refused with `pid_config_error::non_finite_integral`, and the integrator is left bitwise unchanged.

A velocity-form controller has no integrator, so `set_integral` is not available there: calling it is a compile error. Seed a velocity-form controller through the tracking overload of `compute`, which re-seeds the accumulated output (see [velocity-form](velocity-form.md)).

### metric (requires `perf_assessment`)

```cpp
template <typename Metric>
const vector_t& metric() const;
```

Returns the accumulated value of a performance metric. `Metric` must be one of the types in the `perf_assessment` pack (e.g., `IAE`, `ISE`, `ITAE`, `oscillation_detect`). For `oscillation_detect` the value is the detector state `x` per channel: the exponentially weighted count of load disturbances (see [performance](performance.md#oscillation-detection)).

### oscillating (requires `perf_assessment` with `oscillation_detect`)

```cpp
bool oscillating() const;
```

Returns true while the weighted load-disturbance count `x` exceeds `load_count_limit` in any channel. It is a query of the current state, not a latch: it clears by itself as `x` decays. Under a sustained oscillation that holds `x` near the limit it toggles at roughly the half-period cadence; apply hysteresis to `metric<oscillation_detect>()` for a steady signal.

### reset_metrics (requires `perf_assessment`)

```cpp
void reset_metrics();
```

Resets all performance metric accumulators to zero, together with the time base `ITAE` weights by. The controller state is not touched.

## Usage Example

```cpp
// PID step response with anti-windup and derivative filtering.
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'output', '' using 1:3 with lines title 'control'"

#include <ctrlpp/control/pid.h>

#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1,
        ctrlpp::anti_windup<ctrlpp::back_calc>,
        ctrlpp::deriv_filter>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(4.0);
    cfg.ki = Vec::Constant(2.0);
    cfg.kd = Vec::Constant(0.5);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);
    cfg.template policy<ctrlpp::anti_windup<ctrlpp::back_calc>>().kb = {2.0};
    cfg.template policy<ctrlpp::deriv_filter>().n = {10.0};

    auto created = Pid::create(cfg);
    if (!created.has_value()) return 1;
    auto& ctrl = *created;

    double y = 0.0;
    constexpr double dt = 0.01;

    for (int k = 0; k < 500; ++k)
    {
        double t = k * dt;
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y);
        // A refused cycle produced no command; the caller decides what the
        // actuator does. This sample stops.
        auto step = ctrl.compute(sp, meas, dt);
        if (!step.has_value()) return 1;
        const auto& u = *step;
        y = 0.95 * y + 0.05 * u[0];
        std::cout << t << "," << y << "," << u[0] << "\n";
    }
}
```

## Policy Pages

- [anti-windup](anti-windup.md)<br/> back-calculation, clamping, and conditional integration
- [derivative-filter](derivative-filter.md)<br/>low-pass filtered derivative term
- [setpoint-filter](setpoint-filter.md)<br/>setpoint and process variable filtering, and the setpoint weights
- [velocity-form](velocity-form.md)<br/>incremental (velocity) form PID
- [isa-form](isa-form.md)<br/>ISA standard form (Ti, Td parametrization)
- [feed-forward](feed-forward.md)<br/>additive feed-forward policy
- [rate-limit](rate-limit.md)<br/>control signal rate limiting
- [performance](performance.md)<br/>IAE, ISE, ITAE, and oscillation detection

## See Also

- [lqr](../lqr.md)<br/>optimal state-feedback for multi-variable systems
- [PID Composition Guide](../../../guides/pid/composition.md)<br/>combining multiple policies
- [PID Theory](../../../background/pid.md)<br/>mathematical background
