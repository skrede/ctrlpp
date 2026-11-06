# mrac

Model Reference Adaptive Controller with gradient adaptation and compile-time robustification policy selection. The controller tracks a user-supplied discrete-time reference model by adapting state-feedback and feedforward gain matrices online. Supports SISO (NX=NU=1) and MIMO configurations with matrix adaptation gains. Robustification modes (dead-zone, sigma-modification, e-modification) prevent parameter drift under noise or unmodeled dynamics.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::mrac_controller<Scalar, NX, NU, Robustification>` | `#include <ctrlpp/control/mrac.h>` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type (e.g. `double`, `float`) |
| `NX` | `std::size_t`, default `1` | State dimension |
| `NU` | `std::size_t`, default `1` | Input dimension |
| `Robustification` | policy tag, default `no_robustification` | Robustification mode (`no_robustification`, `dead_zone`, `sigma_modification`, `e_modification`). Any other type with a nested `config` (or `config<Scalar>`) compiles and adapts nothing. |

## Type Aliases

```cpp
using config_type  = mrac_config<Scalar, NX, NU, Robustification>;
using state_type   = Vector<Scalar, NX>;
using input_type   = Vector<Scalar, NU>;
using theta_x_type = Matrix<Scalar, NU, NX>;
using theta_r_type = Matrix<Scalar, NU, NU>;
```

## Config: mrac_config

```cpp
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU,
          typename Robustification = no_robustification>
struct mrac_config;
```

Header: `#include <ctrlpp/control/mrac_config.h>`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `reference_model` | `discrete_state_space<Scalar, NX, NU, NX>` | zero | Discrete-time reference model (A_m, B_m, C_m, D_m) |
| `gamma_x` | `Matrix<Scalar, NX, NX>` | zero | SPD adaptation gain for state parameters (user must set; SPD not enforced; the zero default disables the gradient term) |
| `gamma_r` | `Matrix<Scalar, NU, NU>` | zero | SPD adaptation gain for reference parameters (user must set; SPD not enforced; the zero default disables the gradient term) |
| `sign_b` | `Matrix<Scalar, NU, NU>` | Identity | Sign-definite part of the plant's input gain, applied on the left of the update |
| `theta_x_0` | `Matrix<Scalar, NU, NX>` | zero | Initial state parameter estimate |
| `theta_r_0` | `Matrix<Scalar, NU, NU>` | zero | Initial reference parameter estimate |
| `x_model_0` | `Vector<Scalar, NX>` | zero | Initial reference model state |
| `W` | `Matrix<Scalar, NX, NX>` | Identity | Weight of the error norm `sqrt(e' W e)` that `dead_zone` and `e_modification` read (identity gives the Euclidean norm) |
| `robustification` | `Robustification::config` (`Robustification::config<Scalar>` for `dead_zone`, `sigma_modification` and `e_modification`) | value-initialized | Robustification options (see below); the `threshold`, `sigma` and `delta` fields default to zero |

The adaptation gains `gamma_x` and `gamma_r` are meant to be symmetric positive definite, as the Lyapunov argument behind the law takes them; neither is checked. Larger values give faster adaptation but may cause oscillation. The `sign_b` matrix captures the sign structure of the unknown plant input matrix; for most applications with positive-definite B_p, the default identity suffices.

With `e` the tracking error after the reference model has been propagated, the update is

```
theta_x -= sign_b * (B_m^T e) * x^T * gamma_x  [+ leakage]
theta_r -= sign_b * (B_m^T e) * r^T * gamma_r  [+ leakage]
```

where `B_m` is `reference_model.B` and the leakage term is the robustification's (below). The assumptions this law carries are listed under [Limitations](#limitations).

## Robustification Modes

Robustification prevents unbounded parameter drift caused by noise or unmodeled dynamics (the Rohrs counterexample). The robustification mode is selected at compile time via the fourth template parameter.

### Dead-Zone (dead_zone)

Freezes adaptation when the weighted tracking error norm `sqrt(e' W e)` is at or below a threshold, preventing drift from measurement noise in steady state.

**Config field:** `robustification.threshold` (Scalar)<br/> error norm at or below which adaptation is frozen.

**When to use:** Known noise bound. Set threshold slightly above the expected noise level. Most aggressive at preventing small-signal drift, but requires a noise bound estimate.

**Theory:** Slotine & Li, *Applied Nonlinear Control*, Sec. 8.6.

```cpp
using config = ctrlpp::mrac_config<double, 2, 2, ctrlpp::dead_zone>;
config cfg{};
// ... set reference_model, gamma_x, gamma_r ...
cfg.robustification.threshold = 0.01;
```

### Sigma-Modification (sigma_modification)

Subtracts `sigma * theta` from each parameter matrix every cycle, decaying the parameters toward zero when the adaptation gradient is small.

**Config field:** `robustification.sigma` (Scalar)<br/> leakage rate (typical range 0.001--0.1).

**When to use:** Unknown noise bound. Provides global boundedness without needing a noise estimate. Simple to tune: one scalar parameter.

**Theory:** Ioannou & Sun, *Robust Adaptive Control*, Ch. 8.

```cpp
using config = ctrlpp::mrac_config<double, 2, 2, ctrlpp::sigma_modification>;
config cfg{};
// ... set reference_model, gamma_x, gamma_r ...
cfg.robustification.sigma = 0.01;
```

### E-Modification (e_modification)

Subtracts `delta * sqrt(e' W e) * theta` from each parameter matrix every cycle, giving stronger leakage when the error is large and vanishing leakage at equilibrium.

**Config field:** `robustification.delta` (Scalar)<br/> error-proportional leakage gain.

**When to use:** Desire error-proportional leakage that vanishes at equilibrium, preserving the ideal adaptation law when tracking is accurate.

**Theory:** Narendra & Annaswamy, *Stable Adaptive Systems*, Ch. 8.

```cpp
using config = ctrlpp::mrac_config<double, 2, 2, ctrlpp::e_modification>;
config cfg{};
// ... set reference_model, gamma_x, gamma_r ...
cfg.robustification.delta = 0.05;
```

## Class Methods

### Constructor

```cpp
explicit mrac_controller(const config_type& cfg);
```

Constructs the controller from a configuration struct. Stores the configuration and initializes the reference-model state and the adapted parameters to `x_model_0`, `theta_x_0` and `theta_r_0`. Nothing is validated; a non-finite initial value is reported by the first `evaluate`.

### evaluate

```cpp
expected<input_type, mrac_step_error> evaluate(const state_type& x, const input_type& r);
```

Computes the control output `u = theta_x * x + theta_r * r` from the parameters held before the cycle, then propagates the reference model and updates the adapted parameter matrices by the law above. Call once per time step. On success the result carries the control vector.

The cycle is classified before any member is written, so a rejected cycle leaves the reference-model state, the tracking error and both parameter matrices bitwise unchanged, and a following valid cycle produces exactly what it would have produced had the rejected one never been attempted.

**This guard matters more here than on an estimator's update.** An estimator that admits one bad sample is re-driven towards the truth by the samples that follow. The parameter matrices are not: the adaptation only ever subtracts an increment into them and nothing re-derives them from data, so a single admitted non-finite sample makes them non-finite forever. Reconstruction or `reset` is the only way back.

**A refused cycle produced no command.** The actuator will be driven by something regardless, and the caller must choose what: hold the command the last successful cycle returned, drive a configured safe value, or fail over. The controller does not choose, because the right answer is a property of the plant.

The enumerators follow the adaptation's own data flow, so the cause named is always the most upstream one.

| Enumerator | Cause | Where the caller repairs it |
|---|---|---|
| `non_finite_reference_state` | The internal reference-model state is already non-finite. It is the root cause whenever it holds: the tracking error is the plant state minus this, so it would poison both parameter matrices. A committed cycle never leaves it non-finite, so this traces to a non-finite `x_model_0`. | The initial reference-model state. |
| `non_finite_parameters` | An adaptive parameter matrix is already non-finite. Permanent until `reset`, because nothing re-derives it. A committed cycle never leaves them non-finite, so this traces to a non-finite `theta_x_0` or `theta_r_0`. | The initial parameters. |
| `non_finite_state` | The supplied plant state has a non-finite component. | The sensor or estimator. |
| `non_finite_reference` | The supplied reference command has a non-finite component. | The command generator -- a different subsystem from the one that measures the plant. |
| `non_finite_result` | The cycle's own arithmetic produced a non-finite value: the propagated reference-model state, the tracking error, the command, the projected error `B_m^T e` or either updated parameter matrix. Nothing is committed. A non-finite reference model and finite gains whose product overflows end here. | The reference model, the adaptation gains or the magnitude of the arguments. |

### health

```cpp
mrac_health health() const;
```

Reports whether the controller is carrying non-finite state. Latches until `reset` restores the configured initial parameters.

It is set by the cycle that is refused with `non_finite_reference_state` or `non_finite_parameters`. A committed cycle never leaves the carried state non-finite, so the only route in is a non-finite `x_model_0`, `theta_x_0` or `theta_r_0` supplied at construction. A cycle refused on its arguments (`non_finite_state`, `non_finite_reference`) or on its own arithmetic (`non_finite_result`) does not set it.

### theta_x

```cpp
const theta_x_type& theta_x() const;
```

Returns a const reference to the current state parameter estimate matrix (NU x NX).

### theta_r

```cpp
const theta_r_type& theta_r() const;
```

Returns a const reference to the current reference parameter estimate matrix (NU x NU).

### tracking_error

```cpp
const state_type& tracking_error() const;
```

Returns a const reference to the current tracking error vector (x - x_model).

### x_model

```cpp
const state_type& x_model() const;
```

Returns a const reference to the current reference model state vector.

### reset

```cpp
void reset();
```

Resets adapted parameters and model state to the initial values from the config, the tracking error to zero and `health()` to `mrac_health::ok`.

## Usage Examples

### SISO: First-Order Plant Tracking

```cpp
// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader;
//   plot '<./ctrlpp_mrac_01_tracking' using 1:2 with lines title 'reference',
//        '' using 1:3 with lines title 'plant',
//        '' using 1:4 with lines title 'model'"

#include "ctrlpp/control/mrac.h"

#include <iomanip>
#include <iostream>

int main()
{
    using controller = ctrlpp::mrac_controller<double>;
    using config = controller::config_type;

    ctrlpp::discrete_state_space<double, 1, 1, 1> ref_model{};
    ref_model.A(0, 0) = 0.9;
    ref_model.B(0, 0) = 0.1;
    ref_model.C(0, 0) = 1.0;
    ref_model.D(0, 0) = 0.0;

    config cfg{};
    cfg.reference_model = ref_model;
    cfg.gamma_x << 0.5;
    cfg.gamma_r << 0.5;

    controller ctrl(cfg);

    constexpr double a_p = 0.8, b_p = 0.5;
    double x_plant = 0.0;

    ctrlpp::Vector<double, 1> r;
    r[0] = 1.0;

    for (int k = 0; k < 500; ++k) {
        ctrlpp::Vector<double, 1> x;
        x[0] = x_plant;
        // A refused cycle produced no command AND left the adaptive
        // parameters untouched; the caller decides what the actuator does.
        auto step = ctrl.evaluate(x, r);
        if (!step.has_value()) return 1;
        const auto& u = *step;

        std::cout << std::fixed << std::setprecision(6)
                  << k << "," << 1.0 << "," << x_plant << ","
                  << ctrl.x_model()[0] << "," << u[0] << ","
                  << ctrl.theta_x()(0, 0) << ","
                  << ctrl.theta_r()(0, 0) << "\n";

        x_plant = a_p * x_plant + b_p * u[0];
    }
}
```

### MIMO: Two-Channel Tracking

```cpp
// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader;
//   plot '<./ctrlpp_mrac_02_mimo_tracking' using 1:2 with lines title 'ref_ch1',
//        '' using 1:3 with lines title 'plant_ch1',
//        '' using 1:4 with lines title 'model_ch1',
//        '' using 1:5 with lines title 'ref_ch2',
//        '' using 1:6 with lines title 'plant_ch2',
//        '' using 1:7 with lines title 'model_ch2'"

#include "ctrlpp/control/mrac.h"

#include <iomanip>
#include <iostream>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 2;

    using controller = ctrlpp::mrac_controller<double, NX, NU>;
    using config = controller::config_type;

    ctrlpp::discrete_state_space<double, NX, NU, NX> ref_model{};
    ref_model.A << 0.9, 0.0,
                   0.0, 0.85;
    ref_model.B << 0.1, 0.0,
                   0.0, 0.15;
    ref_model.C = ctrlpp::Matrix<double, NX, NX>::Identity();
    ref_model.D = ctrlpp::Matrix<double, NX, NU>::Zero();

    config cfg{};
    cfg.reference_model = ref_model;
    cfg.gamma_x = 0.3 * ctrlpp::Matrix<double, NX, NX>::Identity();
    cfg.gamma_r = 0.3 * ctrlpp::Matrix<double, NU, NU>::Identity();

    controller ctrl(cfg);

    ctrlpp::Matrix<double, NX, NX> A_p;
    A_p << 0.8, 0.1,
           0.0, 0.7;
    ctrlpp::Matrix<double, NX, NU> B_p;
    B_p << 0.5, 0.0,
           0.0, 0.4;

    ctrlpp::Vector<double, NX> x_plant = ctrlpp::Vector<double, NX>::Zero();
    ctrlpp::Vector<double, NU> r;
    r << 1.0, 0.5;

    for (int k = 0; k < 500; ++k) {
        // A refused cycle produced no command AND left the adaptive
        // parameters untouched; the caller decides what the actuator does.
        auto step = ctrl.evaluate(x_plant, r);
        if (!step.has_value()) return 1;
        const auto& u = *step;

        std::cout << std::fixed << std::setprecision(6)
                  << k << ","
                  << r[0] << "," << x_plant[0] << "," << ctrl.x_model()[0] << ","
                  << r[1] << "," << x_plant[1] << "," << ctrl.x_model()[1] << ","
                  << u[0] << "," << u[1] << "\n";

        x_plant = A_p * x_plant + B_p * u;
    }
}
```

## Limitations

### The adaptation fixes the Lyapunov weight to the identity

The update projects the tracking error through `B_m^T` only: it is the gradient `B_m^T P e` with the Lyapunov weight `P` fixed to the identity. For the discrete-time reference model the controller takes, `P = I` is a valid Lyapunov weight exactly when `I - A_m^T A_m` is positive definite, that is when the spectral norm of `reference_model.A` is below one. A stable `A_m` outside that set needs a different `P`; no configurable `P` weighting is provided, and the condition is not checked.

### The update carries no sample-time factor

The update has no `dt` factor, so `gamma_x` and `gamma_r` absorb the sample time: rescale them proportionally if the sample rate changes. It is not normalized by the signal magnitude either, so the step it takes grows with the square of the state and reference magnitudes, and the gains are chosen for the signal range expected.

### MIMO adaptation assumes the reference input matrix is the plant's up to a positive scalar

The gradient the decrease argument needs is `B_p^T P e`, with `B_p` the plant's input matrix, but the configuration has no plant input matrix and the update uses `reference_model.B` in its place. Writing `B_m = B_p K_r*`, the update direction differs from the gradient by `K_r*^T`. For a single input `B_m^T e` is `b_m e`, a positive multiple of the error the gains absorb when `b_m > 0`, with `sign_b` carrying the sign of the plant's input gain; a negative `b_m` reverses the update's direction. With several inputs the decrease argument holds only when `reference_model.B` equals the plant input matrix times a positive scalar; for any other `B_m` the argument does not cover the update, though the adaptation may still converge.

## See Also

- [Adaptive Control Theory](../../background/adaptive-control.md)<br/> background theory covering MRAC and L1
- [PID](pid/README.md)<br/> structural parallel (stateful controller with compile-time policy composition)
- [state_space](../model/state-space.md)<br/> reference model representation
- [L1 Adaptive Control](l1.md)<br/> L1 adaptive controller with low-pass filtered control output
