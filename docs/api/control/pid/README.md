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
| `pid_config_error`, `pid_step_error`, `pid_health` | `#include <ctrlpp/control/pid_error.h>` (included by `pid.h`) |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (`double`, `float`, `long double`) |
| `NX` | `std::size_t` | State dimension (unused by PID internals, present for interface uniformity) |
| `NU` | `std::size_t` | Input dimension (unused by PID internals, present for interface uniformity) |
| `NY` | `std::size_t` | Output/measurement dimension &mdash; determines the vector size for gains, setpoints, and measurements |
| `Policies...` | policy types | Zero or more policy types from `pid_policies.h` |

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
| `derivative_on_error` | `bool` | `false` | If true, differentiates error; if false, differentiates measurement (avoids derivative kick) |
| `b` | `Vector<Scalar, NY>` | `1.0` | Proportional setpoint weight: P term uses `b*r - y` |
| `c` | `Vector<Scalar, NY>` | `1.0` | Derivative setpoint weight: D term uses `c*r - y` |

### Policy Config Access

```cpp
cfg.template policy<anti_windup<back_calc>>().kb = {1.0};
cfg.template policy<deriv_filter>().n = {10.0};
cfg.template policy<setpoint_filter>().tf = {0.1};
cfg.template policy<pv_filter>().tf = {0.05};
cfg.template policy<rate_limit>().rate_max = {0.5};
cfg.template policy<feed_forward<FF>>().ff_func = my_ff_func;
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
| `feed_forward<Callable>` | Additive feed-forward | `ff_func` | [feed-forward](feed-forward.md) |
| `rate_limit` | Output rate limiting | `rate_max` | [rate-limit](rate-limit.md) |
| `perf_assessment<Metrics...>` | Online performance metrics | (none) | [performance](performance.md) |

## create

```cpp
static auto create(const config_type& cfg) -> expected<pid, pid_config_error>;
```

Validates the configuration and constructs the controller, or reports the first fault it finds. A configuration that cannot work is refused where it is given rather than discovered on the first cycle. The checks run in the order of the table, so a configuration with several faults reports the earliest; a check on a policy field runs only when that policy is composed. `create` is the only construction path: the constructor is private, so no configuration reaches the controller without these checks.

```cpp
auto made = ctrlpp::pid<double, 1>::create(cfg);
if(!made) {
    // made.error() names the first field that cannot work
}
```

| Enumerator | Cause | Repair |
|---|---|---|
| `non_finite_gain` | A `kp`, `ki` or `kd` component is `NaN` or infinite. In ISA form `ki` is the integral time, and a `NaN` or negative-infinite one is refused. | The tuning that produced the gains. |
| `zero_integral_time` | In ISA form an integral time `ki` is zero, which is an infinite integral gain Kp / Ti, not "off". | A positive integral time. |
| `invalid_output_limits` | An output limit is `NaN`, `output_min` exceeds `output_max` in some channel, `output_min` is `+inf`, or `output_max` is `-inf`. Equal limits are accepted. | The actuator range. |
| `non_finite_setpoint_weight` | A `b` or `c` component is `NaN` or infinite. | The setpoint weights. |
| `invalid_filter_time_constant` | A `setpoint_filter` or `pv_filter` time constant is `NaN`, negative or `+inf`. Zero means no filter. | The filter's `tf`. |
| `invalid_rate_limit` | A `rate_limit` `rate_max` is `NaN` or negative. Zero and `+inf` mean no limit. | The rate limit. |
| `invalid_back_calculation_gain` | An explicit `anti_windup<back_calc>` gain `kb` is `NaN` or infinite. Zero selects the default gain. | The tracking gain. |
| `invalid_error_threshold` | An `anti_windup<conditional_integration>` `error_threshold` is `NaN`. | The integration threshold. |
| `invalid_derivative_filter` | A `deriv_filter` `n` is `NaN` or infinite, or an explicit `tf` is `NaN`, infinite, zero or negative. | The filter divisor or time constant. |
| `undefined_tracking_time` | `anti_windup<back_calc>` leaves a channel's `kb` at zero, the channel has integral action (`ki != 0`), and its nonzero gains do not share a sign, so the default tracking time does not exist. | An explicit `kb`, or gains of one sign. |
| `undefined_derivative_filter` | `deriv_filter` is composed, a channel has derivative action (`kd != 0`) and no `tf`, and its default Tf = kd/(kp*n) is not a positive time: `kp` is zero, `n` is zero or negative, or `kd` and `kp` differ in sign. | An explicit `tf`, or a positive `n`. |
| `unrepresentable_gain` | Each field is finite, but a quantity derived from them is not representable in `Scalar`: in ISA form `kp` over the integral time, or `kp` times the derivative time, overflows; the default derivative filter time constant kd/(kp*n) is not a positive finite number; or the default back-calculation gain is not finite. | Smaller gains, an explicit `tf` or `kb`, or a wider `Scalar`. |

## Methods

### compute

```cpp
auto compute(const vector_t& sp, const vector_t& meas, Scalar dt)
    -> expected<vector_t, pid_step_error>;
```

Computes the control output given a setpoint, measurement, and time step. On success the result carries the clamped output vector.

The cycle is classified before any member is written, so a rejected cycle leaves every piece of carried state -- the integrator, both error histories, the input filter states, the accumulated output -- bitwise unchanged, and the caller may retry on the next sample.

**A refused cycle produced no command.** That is a different situation from an estimator declining a measurement, where the estimate simply stands: an actuator will be driven by something regardless of what `compute` returns, and the caller must choose what. The three defensible responses are to hold the command the last successful cycle returned, to command a configured safe value, or to fail over to a redundant channel. The controller does not choose, because the right answer is a property of the plant: holding is correct for a slow thermal loop and dangerous for an unstable attitude loop.

A non-positive step used to be answered with the previous output on the success path, which the caller could not distinguish from a freshly computed command -- so a stopped clock read as a steady loop. It is now `invalid_timestep`.

| Enumerator | Cause | Where the caller repairs it |
|---|---|---|
| `non_finite_state` | Some piece of carried state is already non-finite when the cycle begins. | `reset`, or `set_integral` with a finite value for the integrator alone. Reported first because it is upstream of every argument. |
| `invalid_timestep` | The step is not a positive finite duration. A `NaN` step cannot be caught by `dt <= 0` alone, because every comparison against a `NaN` is false. | The caller's clock. |
| `non_finite_setpoint` | The commanded setpoint has a non-finite component. | The trajectory, outer loop, or operator input that generates the reference. |
| `non_finite_measurement` | The process variable has a non-finite component. | The sensor or estimator. Separate from the setpoint even though the error is their difference, because the two arrive from different subsystems. |
| `non_finite_tracking_signal` | The tracking signal handed to the four-argument overload is non-finite. | The signal being tracked. |

### compute (with tracking signal)

```cpp
auto compute(const vector_t& sp, const vector_t& meas, Scalar dt,
             const vector_t& tracking_signal)
    -> expected<vector_t, pid_step_error>;
```

Overload for controller output tracking (bumpless transfer in cascade configurations). After the cycle succeeds the controller is re-seeded so its output tracks an external signal, and the returned value is that cycle's own computed command, or increment in velocity form. Position form re-seeds the integrator; velocity form re-seeds the accumulated output.

In position form, after the cycle succeeds the integrator is assigned `I = clamp(w, output_min, output_max) - (P + D + FF)`, where `w` is the tracking signal and `P + D + FF` is the sum of the proportional, derivative and feed-forward terms of the cycle that formed the output, taken before rate limiting and saturation. It is not reconstructed from the returned command, which may have been saturated or rate-limited, and the anti-windup step of that cycle has no influence on it. A target outside the output limits is clamped first, so the integrator never holds an offset that drives the unconstrained command past a limit. For example, with `kp = 1`, `ki = 0.1`, limits of `+-10`, setpoint 20 and measurement 0, tracking 5 leaves the integrator at -15, and tracking 50 leaves it at `10 - 20 = -10`.

In velocity form, after the cycle succeeds the accumulated output is assigned `clamp(w, output_min, output_max)`, and the next cycle's increment is accumulated from it. For example, with `kp = 1` and limits of `+-10`, a cycle at setpoint 1 tracking 9.5 re-seeds the accumulated output to 9.5, and a following setpoint step to 3 emits 0.5 and reports `saturated()`.

The tracking signal is checked *before* the cycle is delegated, because it is written straight into carried state once the cycle succeeds: admitting a non-finite one would destroy that state after the cycle had already committed the rest. A rejection from the delegated cycle is forwarded unchanged with its own specific cause, and the tracking assignment does not run.

### health

```cpp
auto health() const -> pid_health;
```

Reports whether the controller is still carrying non-finite state from an earlier cycle -- a question a per-cycle result cannot answer, because it outlives the call. Latches until `reset`.

A **rejected cycle does not set it**: a rejection mutates nothing, so it leaves the controller exactly as healthy as it was. No known public route leads to `non_finite_carried_state`: `create`, `set_params` and `set_integral` refuse non-finite values, and a cycle whose result is non-finite is rejected and leaves no state behind. The latch is kept as defense in depth.

### set_params

```cpp
auto set_params(const config_type& new_cfg) -> expected<void, pid_config_error>;
```

Refuses the configuration with the same checks, in the same order, as `create`; a refusal leaves the current configuration, the integrator and every other piece of state bitwise unchanged. On success it updates tuning parameters with bumpless transfer. The integral state is stored in output units and enters the output directly, so its contribution is already continuous across a gain change; `set_params` leaves the integral state untouched, which is what makes gain scheduling bumpless. No rescaling is applied (rescaling would step the integral contribution instead of preserving it).

### error

```cpp
const vector_t& error() const;
```

Returns the most recent error vector (setpoint minus measurement, after filtering).

### integral

```cpp
const vector_t& integral() const;
```

Returns the current integral accumulator state.

### params

```cpp
const config_type& params() const;
```

Returns a const reference to the current config.

### saturated

```cpp
bool saturated() const;
```

Returns true if the output was clipped on the most recent `compute()` call. In velocity form it is true when the accumulated output was clamped to the output limits on that call.

### reset

```cpp
void reset();
```

Resets all internal state (integral, previous error, filter states, performance metrics) to zero.

### freeze_integral

```cpp
void freeze_integral(bool freeze = true);
```

Freezes or unfreezes the integral accumulator. While frozen, the integral term does not change, including through back-calculation anti-windup, which resumes once the accumulator is unfrozen.

### set_integral

```cpp
auto set_integral(const vector_t& val) -> expected<void, pid_config_error>;
```

Directly sets the integral accumulator to a given value. A `NaN` or infinite component is refused with `pid_config_error::non_finite_integral`, and the integrator is left bitwise unchanged.

### metric (requires `perf_assessment`)

```cpp
template <typename Metric>
const vector_t& metric() const;
```

Returns the accumulated value of a performance metric. `Metric` must be one of the types in the `perf_assessment` pack (e.g., `IAE`, `ISE`, `ITAE`, `oscillation_detect`).

### oscillating (requires `perf_assessment` with `oscillation_detect`)

```cpp
bool oscillating() const;
```

Returns true if the zero-crossing rate exceeds the configured threshold.

### reset_metrics (requires `perf_assessment`)

```cpp
void reset_metrics();
```

Resets all performance metric accumulators to zero.

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
- [setpoint-filter](setpoint-filter.md)<br/>setpoint weighting and process variable filtering
- [velocity-form](velocity-form.md)<br/>incremental (velocity) form PID
- [isa-form](isa-form.md)<br/>ISA standard form (Ti, Td parametrization)
- [feed-forward](feed-forward.md)<br/>additive feed-forward policy
- [rate-limit](rate-limit.md)<br/>control signal rate limiting
- [performance](performance.md)<br/>IAE, ISE, ITAE, and oscillation detection

## See Also

- [lqr](../lqr.md)<br/>optimal state-feedback for multi-variable systems
- [PID Composition Guide](../../../guides/pid/composition.md)<br/>combining multiple policies
- [PID Theory](../../../background/pid.md)<br/>mathematical background
