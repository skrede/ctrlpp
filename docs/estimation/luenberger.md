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
void update(const output_vector_t& z);
```

Corrects the state estimate with measurement: x = x + L(z - Cx).

### state

```cpp
auto state() const -> const state_vector_t&;
```

Returns the current state estimate.

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

Resets the state estimate to a new initial value.

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
        obs.update(z);

        auto est = obs.state();
        std::cout << k * 0.1 << "," << x_true(0) << "," << est(0) << "\n";
    }
}
```

## See Also

- [kalman](kalman.md) -- optimal stochastic observer
- [observer-policy](observer-policy.md) -- concept satisfied by this type
- [place](../control/place.md) -- pole placement for observer gain design
