# mrac

Model Reference Adaptive Controller with Lyapunov-based adaptation and compile-time robustification policy selection. The controller tracks a user-supplied discrete-time reference model by adapting state-feedback and feedforward gain matrices online. Supports SISO (NX=NU=1) and MIMO configurations with matrix adaptation gains. Robustification modes (dead-zone, sigma-modification, e-modification) prevent parameter drift under noise or unmodeled dynamics.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::mrac_controller<Scalar, NX, NU, Robustification>` | `#include <ctrlpp/control/mrac.h>` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
| `NX` | `std::size_t`, default `1` | State dimension |
| `NU` | `std::size_t`, default `1` | Input dimension |
| `Robustification` | policy tag, default `no_robustification` | Robustification mode (`no_robustification`, `dead_zone`, `sigma_modification`, `e_modification`) |

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
template <typename Scalar, std::size_t NX, std::size_t NU,
          typename Robustification = no_robustification>
struct mrac_config;
```

Header: `#include <ctrlpp/control/mrac_config.h>`

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `reference_model` | `discrete_state_space<Scalar, NX, NU, NX>` | zero | Discrete-time reference model (A_m, B_m, C_m, D_m) |
| `gamma_x` | `Matrix<Scalar, NX, NX>` | zero | SPD adaptation gain for state parameters (user must set; SPD not enforced) |
| `gamma_r` | `Matrix<Scalar, NU, NU>` | zero | SPD adaptation gain for reference parameters (user must set; SPD not enforced) |
| `sign_b` | `Matrix<Scalar, NU, NU>` | Identity | Sign-definite part of high-frequency gain matrix |
| `theta_x_0` | `Matrix<Scalar, NU, NX>` | zero | Initial state parameter estimate |
| `theta_r_0` | `Matrix<Scalar, NU, NU>` | zero | Initial reference parameter estimate |
| `x_model_0` | `Vector<Scalar, NX>` | zero | Initial reference model state |
| `W` | `Matrix<Scalar, NX, NX>` | Identity | Error weight matrix for weighted norm in robustification (defaults to identity for standard Euclidean norm) |
| `robustification` | (policy-dependent) | &mdash; | Robustification options (see below) |

The adaptation gains `gamma_x` and `gamma_r` must be symmetric positive definite for Lyapunov stability guarantees. Larger values give faster adaptation but may cause oscillation. The `sign_b` matrix captures the sign structure of the unknown plant input matrix; for most applications with positive-definite B_p, the default identity suffices.

**Adaptation law assumptions.** The parameter update projects the tracking error through `B_m^T` only. This is the Lyapunov gradient `B_m^T P e` with `P = I`, which is exact when the reference model is chosen so that `A_m^T + A_m` is negative definite (`P = I` solves `A_m^T P + P A_m = -Q`). For a general stable `A_m`, a different `P` would be required; no configurable `P` weighting is provided. The update also carries no explicit `dt` factor, so `gamma_x` and `gamma_r` absorb the sample time: rescale them proportionally if the sample rate changes.

## Robustification Modes

Robustification prevents unbounded parameter drift caused by noise or unmodeled dynamics (the Rohrs counterexample). The robustification mode is selected at compile time via the fourth template parameter.

### Dead-Zone (dead_zone)

Freezes adaptation when the tracking error norm falls below a threshold, preventing drift from measurement noise in steady state.

**Config field:** `robustification.threshold` (Scalar)<br/> error norm below which adaptation is frozen.

**When to use:** Known noise bound. Set threshold slightly above the expected noise level. Most aggressive at preventing small-signal drift, but requires a noise bound estimate.

**Theory:** Slotine & Li, *Applied Nonlinear Control*, Sec. 8.6.

```cpp
using config = ctrlpp::mrac_config<double, 2, 2, ctrlpp::dead_zone>;
config cfg{};
// ... set reference_model, gamma_x, gamma_r ...
cfg.robustification.threshold = 0.01;
```

### Sigma-Modification (sigma_modification)

Adds a leakage term proportional to the current parameter values, gradually decaying parameters toward zero when the adaptation gradient is small.

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

Adds a leakage term proportional to both the parameter values and the tracking error norm, giving stronger leakage when the error is large and vanishing leakage at equilibrium.

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

Constructs the controller from a configuration struct. Stores the reference model and initializes adapted parameters to their initial values.

### evaluate

```cpp
auto evaluate(const state_type& x, const input_type& r) -> input_type;
```

Computes the control output `u = theta_x * x + theta_r * r` and updates the adapted parameter matrices using the Lyapunov-based adaptation law. Call once per time step. The reference model is propagated internally.

### theta_x

```cpp
auto theta_x() const -> const theta_x_type&;
```

Returns a const reference to the current state parameter estimate matrix (NU x NX).

### theta_r

```cpp
auto theta_r() const -> const theta_r_type&;
```

Returns a const reference to the current reference parameter estimate matrix (NU x NU).

### tracking_error

```cpp
auto tracking_error() const -> const state_type&;
```

Returns a const reference to the current tracking error vector (x - x_model).

### x_model

```cpp
auto x_model() const -> const state_type&;
```

Returns a const reference to the current reference model state vector.

### reset

```cpp
void reset();
```

Resets adapted parameters and model state to the initial values from the config.

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
        auto u = ctrl.evaluate(x, r);

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
        auto u = ctrl.evaluate(x_plant, r);

        std::cout << std::fixed << std::setprecision(6)
                  << k << ","
                  << r[0] << "," << x_plant[0] << "," << ctrl.x_model()[0] << ","
                  << r[1] << "," << x_plant[1] << "," << ctrl.x_model()[1] << ","
                  << u[0] << "," << u[1] << "\n";

        x_plant = A_p * x_plant + B_p * u;
    }
}
```

## See Also

- [Adaptive Control Theory](../../background/adaptive-control.md)<br/> background theory covering MRAC and L1
- [PID](pid/README.md)<br/> structural parallel (stateful controller with compile-time policy composition)
- [state_space](../model/state-space.md)<br/> reference model representation
- [L1 Adaptive Control](l1.md)<br/> L1 adaptive controller with low-pass filtered control output
