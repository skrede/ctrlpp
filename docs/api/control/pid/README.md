# pid

A policy-based PID controller with compile-time feature composition. The template
signature `pid<Scalar, NX, NU, NY, Policies...>` lets you opt into exactly the
features you need -- anti-windup, derivative filtering, setpoint weighting,
velocity form, ISA form, feed-forward, rate limiting, and performance assessment --
without paying for what you don't use.

The base PID provides proportional, integral, and derivative action with output
clamping and setpoint weighting (b, c parameters). Each policy extends the
controller with one additional behaviour, composed at compile time via variadic
template packs.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::pid<Scalar, NX, NU, NY, Policies...>` | `#include <ctrlpp/control/pid.h>` |
| (convenience) | `#include <ctrlpp/pid.h>` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (`double`, `float`, `long double`) |
| `NX` | `std::size_t` | State dimension (unused by PID internals, present for interface uniformity) |
| `NU` | `std::size_t` | Input dimension (unused by PID internals, present for interface uniformity) |
| `NY` | `std::size_t` | Output/measurement dimension -- determines the vector size for gains, setpoints, and measurements |
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
| `deriv_filter` | Low-pass filtered derivative | `n` | [derivative-filter](derivative-filter.md) |
| `setpoint_filter` | Setpoint prefilter | `tf` | [setpoint-filter](setpoint-filter.md) |
| `pv_filter` | Process variable filter | `tf` | [setpoint-filter](setpoint-filter.md) |
| `velocity_form` | Incremental output form | (none) | [velocity-form](velocity-form.md) |
| `isa_form` | ISA standard form (Kp, Ti, Td) | (none) | [isa-form](isa-form.md) |
| `feed_forward<Callable>` | Additive feed-forward | `ff_func` | [feed-forward](feed-forward.md) |
| `rate_limit` | Output rate limiting | `rate_max` | [rate-limit](rate-limit.md) |
| `perf_assessment<Metrics...>` | Online performance metrics | (none) | [performance](performance.md) |

## Constructor

```cpp
explicit pid(const config_type& cfg);
```

Constructs the controller from a config struct. Computes internal gains and initialises policy state.

## Methods

### compute

```cpp
auto compute(const vector_t& sp, const vector_t& meas, Scalar dt) -> vector_t;
```

Computes the control output given a setpoint, measurement, and time step. Returns the clamped output vector. If `dt <= 0`, returns the previous output unchanged.

### compute (with tracking signal)

```cpp
vector_t compute(const vector_t& sp, const vector_t& meas, Scalar dt,
                 const vector_t& tracking_signal);
```

Overload for controller output tracking (bumpless transfer in cascade configurations). The tracking signal adjusts the integral term so the controller output tracks an external signal.

### set_params

```cpp
void set_params(const config_type& new_cfg);
```

Updates tuning parameters with bumpless transfer. Rescales the integral state to avoid output discontinuities when gains change.

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

Returns true if the output was clipped on the most recent `compute()` call.

### reset

```cpp
void reset();
```

Resets all internal state (integral, previous error, filter states, performance metrics) to zero.

### freeze_integral

```cpp
void freeze_integral(bool freeze = true);
```

Freezes or unfreezes the integral accumulator. While frozen, the integral term does not change.

### set_integral

```cpp
void set_integral(const vector_t& val);
```

Directly sets the integral accumulator to a given value.

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
    using Pid = ctrlpp::pid<double, 1, 1, 1,
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

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double dt = 0.01;

    for (int k = 0; k < 500; ++k)
    {
        double t = k * dt;
        auto sp = Vec::Constant(1.0);
        auto meas = Vec::Constant(y);
        auto u = ctrl.compute(sp, meas, dt);
        y = 0.95 * y + 0.05 * u[0];
        std::cout << t << "," << y << "," << u[0] << "\n";
    }
}
```

## Policy Pages

- [anti-windup](anti-windup.md) -- back-calculation, clamping, and conditional integration
- [derivative-filter](derivative-filter.md) -- low-pass filtered derivative term
- [setpoint-filter](setpoint-filter.md) -- setpoint weighting and process variable filtering
- [velocity-form](velocity-form.md) -- incremental (velocity) form PID
- [isa-form](isa-form.md) -- ISA standard form (Ti, Td parameterisation)
- [feed-forward](feed-forward.md) -- additive feed-forward policy
- [rate-limit](rate-limit.md) -- control signal rate limiting
- [performance](performance.md) -- IAE, ISE, ITAE, and oscillation detection

## See Also

- [lqr](../lqr.md) -- optimal state-feedback for multi-variable systems
- [PID Composition Guide](../../guides/pid/composition.md) -- combining multiple policies
- [PID Theory](../../reference/pid-theory.md) -- mathematical background
