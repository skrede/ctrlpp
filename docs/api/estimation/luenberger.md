# luenberger_observer

Discrete-time Luenberger state observer with a fixed gain matrix L. Provides deterministic state estimation for linear systems where noise statistics are not available or not needed. The observer gain L is typically designed via pole placement (`place_observer`) to set the convergence rate of the estimation error.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::luenberger_observer<Scalar, NX, NU, NY>` | `#include <ctrlpp/estimation/luenberger.h>` |
| `ctrlpp::luenberger_observer<Scalar, NX, NU, NY>` | `#include <ctrlpp/luenberger.h>` (convenience) |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
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

Constructs the observer from a discrete state-space model, observer gain matrix L, and initial state estimate.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates the state estimate: x = Ax + Bu.

### update

```cpp
auto update(const output_vector_t& z)
    -> ctrlpp::expected<void, luenberger_update_error>;
```

Corrects the state estimate with measurement: x = x + L(z - Cx).

The step is rejected **before the state is assigned**, so a rejected step leaves the state bitwise unchanged and the caller may retry with the next sample. The observer carries no covariance, so there are two causes rather than the three the covariance filters have, and each is an exact domain condition rather than a tuning preference.

| Condition | Error | Why it is a separate cause |
|---|---|---|
| the carried state estimate is already non-finite | `luenberger_update_error::non_finite_state` | The fault is upstream of the measurement, so it is reported ahead of it: a caller told the measurement is bad would replace a working sensor while the real fault sits in the prediction that poisoned the state |
| the supplied measurement has a non-finite component | `luenberger_update_error::non_finite_measurement` | The correction is the single fused expression x + L(z - Cx), so every state component whose gain row is nonzero becomes non-finite and the observer carries that state forward with no path back |

`predict` is deliberately not fallible. Its input is a command the caller already owns and the plant already took, so refusing it would leave the filter with no propagation for a step that happened. A prediction that poisons the carried estimate is reported by [`health`](#health) instead, which the next `update` latches.

### state

```cpp
auto state() const -> const state_vector_t&;
```

Returns the current state estimate.

### health

```cpp
auto health() const -> luenberger_health;
```

Returns the persistent state-health status, one of `luenberger_health::ok` or `luenberger_health::non_finite_estimate`. It answers a question a per-call result cannot, because the question outlives the call: whether the carried estimate is still degraded from a step several samples ago. The status starts at `ok` and latches to `non_finite_estimate` the first time a step finds the carried state already non-finite, which is how a poisoned `predict` becomes visible; [`reset`](#reset) clears it, because it replaces the very state the status describes. A **rejected measurement does not set it**: the rejection mutates nothing, so it leaves the observer healthy. The query carries no discard warning; asking it is optional.

### set_gain

```cpp
void set_gain(const gain_matrix_t& L);
```

Updates the observer gain.

### set_model

```cpp
void set_model(system_t sys);
```

Replaces the state-space model.

### reset

```cpp
void reset(const state_vector_t& x0);
```

Resets the state estimate to a new initial value. This is the one operation that clears a latched [`health`](#health) status, because it replaces the state that status describes.

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

    // Place observer poles at 0.3 +/- 0.1j (fast convergence)
    std::array<std::complex<Scalar>, NX> poles = {
        std::complex<Scalar>{0.3, 0.1},
        std::complex<Scalar>{0.3, -0.1}
    };
    auto L_opt = ctrlpp::place_observer<Scalar, NX, NY>(sys.A, sys.C, poles);

    Eigen::Vector2d x0_est = Eigen::Vector2d::Zero();
    ctrlpp::luenberger_observer<Scalar, NX, NU, NY> obs(sys, *L_opt, x0_est);

    Eigen::Vector2d x_true;
    x_true << 1.0, 0.5;

    for (int k = 0; k < 50; ++k) {
        Eigen::Matrix<Scalar, 1, 1> u;
        u << 0.0;

        Eigen::Matrix<Scalar, 1, 1> z = sys.C * x_true;
        x_true = sys.A * x_true + sys.B * u;

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

## See Also

- [kalman](kalman.md)<br/> optimal stochastic observer
- [observer-policy](observer-policy.md)<br/> concept satisfied by this type
- [place](../control/place.md)<br/> pole placement for observer gain design
- [background/kalman](../../background/kalman.md)<br/> observer theory and comparison with Kalman filter
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> composing observers with controllers
