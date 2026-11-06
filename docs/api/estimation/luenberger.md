# luenberger_observer

Discrete-time Luenberger state observer with a fixed gain matrix L. Provides deterministic state estimation for linear systems where noise statistics are not available or not needed. The observer splits each step into `predict` and `update`, so its error dynamics are those of `A - LCA`, not `A - LC`; [Gain design](#gain-design) shows how to place them with `place_observer`.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::luenberger_observer<Scalar, NX, NU, NY>` | `#include <ctrlpp/estimation/luenberger.h>` |
| `ctrlpp::luenberger_observer<Scalar, NX, NU, NY>` | `#include <ctrlpp/luenberger.h>` (convenience) |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type (e.g. `double`, `float`) |
| `NX` | `std::size_t` | State dimension |
| `NU` | `std::size_t` | Input dimension |
| `NY` | `std::size_t` | Output dimension |

## Type Aliases

```cpp
using state_vector_t  = Eigen::Matrix<Scalar, nx, 1>;
using input_vector_t  = Eigen::Matrix<Scalar, nu, 1>;
using output_vector_t = Eigen::Matrix<Scalar, ny, 1>;
using gain_matrix_t   = Eigen::Matrix<Scalar, nx, ny>;
using system_t        = discrete_state_space<Scalar, NX, NU, NY>;
```

## Constructor

```cpp
luenberger_observer(system_t sys, gain_matrix_t L, state_vector_t x0);
```

Constructs the observer from a discrete state-space model, observer gain matrix L, and initial state estimate. The constructor is not fallible and validates nothing: a non-finite `x0` is refused by the first `update` as `non_finite_state`, and a non-finite `L` or model surfaces as `non_finite_result` or, through `predict`, as `non_finite_state`.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates the state estimate: x = Ax + Bu.

### update

```cpp
ctrlpp::expected<void, luenberger_update_error> update(const output_vector_t& z);
```

Corrects the state estimate with measurement: x = x + L(z - Cx). The feedthrough `D` is not used.

A refused step leaves the state bitwise unchanged, so the caller may retry with the next sample. The observer carries no covariance, so it has no covariance cause; the first two causes are checked before anything is computed and the third on the computed state, which is then discarded. Each is an exact domain condition rather than a tuning preference.

| Condition | Error | Why it is a separate cause |
|---|---|---|
| the carried state estimate is already non-finite | `luenberger_update_error::non_finite_state` | The fault is upstream of the measurement, so it is reported ahead of it: a caller told the measurement is bad would replace a working sensor while the real fault sits in the prediction that poisoned the state |
| the supplied measurement has a non-finite component | `luenberger_update_error::non_finite_measurement` | The correction is the single fused expression x + L(z - Cx), so every state component whose gain row is nonzero becomes non-finite and the observer carries that state forward with no path back |
| the corrected state is non-finite although the carried state and the measurement were finite | `luenberger_update_error::non_finite_result` | The product overflowed, or the gain or `C` carries a non-finite entry; the fault is in the configuration, not in the estimate or the sample |

`predict` is deliberately not fallible. Its input is a command the caller already owns and the plant already took, so refusing it would leave the filter with no propagation for a step that happened. A prediction that poisons the carried estimate is reported by [`health`](#health) instead, which the next `update` latches.

### state

```cpp
const state_vector_t& state() const;
```

Returns the current state estimate.

### health

```cpp
luenberger_health health() const;
```

Returns the persistent state-health status, one of `luenberger_health::ok` or `luenberger_health::non_finite_estimate`. It answers a question a per-call result cannot, because the question outlives the call: whether the carried estimate is still degraded from a step several samples ago. The status starts at `ok` and latches to `non_finite_estimate` when an `update` is refused with `non_finite_state`, which is how a poisoned `predict` becomes visible; [`reset`](#reset) clears it, because it replaces the very state the status describes. An update refused with `non_finite_measurement` or `non_finite_result` does not set it, since neither leaves a non-finite estimate behind.

### set_gain

```cpp
void set_gain(const gain_matrix_t& L);
```

Replaces the observer gain. Not validated.

### set_model

```cpp
void set_model(system_t sys);
```

Replaces the state-space model. Not validated.

### reset

```cpp
void reset(const state_vector_t& x0);
```

Resets the state estimate to a new initial value, unvalidated, and the [`health`](#health) status to `ok`. It is the one operation that clears a latched status, because it replaces the state that status describes.

## Gain design

A step is `predict(u)` followed by `update(z)` with the measurement of the state just predicted, so the corrected estimate obeys `x_k = (I - LC)(A x_{k-1} + B u_{k-1}) + L z_k`, the current-estimator form. Its error evolves as `e_k = (I - LC) A e_{k-1}`, whose eigenvalues are those of `A - LCA`; calling `update` before `predict` gives the same error poles. `place_observer(A, C, poles)` places the eigenvalues of `A - LC`, the prediction-form observer of the [observer background](../../background/observers.md), which this class does not implement. To place the error poles of this observer, pass `CA` as the output matrix: `place_observer(A, C * A, poles)` returns `L` with `eig(A - LCA)` at the requested poles.


## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'true', '' using 1:3 with lines title 'estimate'"

#include <ctrlpp/estimation/luenberger.h>
#include <ctrlpp/control/place.h>
#include <ctrlpp/model/state_space.h>

#include <Eigen/Dense>

#include <array>
#include <complex>
#include <iostream>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2, NU = 1, NY = 1;

    ctrlpp::discrete_state_space<Scalar, NX, NU, NY> sys{};
    sys.A << 1.0, 0.1, 0.0, 1.0;
    sys.B << 0.005, 0.1;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    // Error poles of the predict-update observer at 0.3 +/- 0.1j: eig(A - L C A).
    std::array<std::complex<Scalar>, NX> poles = {
        std::complex<Scalar>{0.3, 0.1},
        std::complex<Scalar>{0.3, -0.1}
    };
    Eigen::Matrix<Scalar, 1, 2> CA = sys.C * sys.A;
    auto L_opt = ctrlpp::place_observer<Scalar, NX, NY>(sys.A, CA, poles);
    if (!L_opt.has_value()) {
        // L_opt.error() is the place_error naming the refusal.
        std::cerr << "observer pole placement refused the design\n";
        return 1;
    }

    Eigen::Vector2d x0_est = Eigen::Vector2d::Zero();
    ctrlpp::luenberger_observer<Scalar, NX, NU, NY> obs(sys, *L_opt, x0_est);

    Eigen::Vector2d x_true;
    x_true << 1.0, 0.5;

    for (int k = 0; k < 50; ++k) {
        Eigen::Matrix<Scalar, 1, 1> u;
        u << 0.0;

        x_true = sys.A * x_true + sys.B * u;
        Eigen::Matrix<Scalar, 1, 1> z = sys.C * x_true;

        obs.predict(u);
        if(!obs.update(z))
        {
            std::cerr << "observer rejected the measurement at step " << k << "\n";
            return 1;
        }

        auto est = obs.state();
        std::cout << k * 0.1 << "," << x_true(0) << "," << est(0) << "\n";
    }
}
```

## Limitations

### The feedthrough matrix is not used

The correction compares `z` with `Cx`. `update` takes no input, so the `D` of the `discrete_state_space` model never enters it, and a plant with `D != 0` is observed as if `D` were zero. A caller with such a plant subtracts `Du` from the measurement before passing it to `update`.

## See Also

- [kalman](kalman.md)<br/> optimal stochastic observer
- [observer-policy](observer-policy.md)<br/> concept satisfied by this type
- [place](../control/place.md)<br/> pole placement for observer gain design
- [background/observers](../../background/observers.md)<br/> observer theory and comparison with the Kalman filter
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> composing observers with controllers
