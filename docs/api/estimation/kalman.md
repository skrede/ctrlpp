# kalman_filter

Linear discrete-time Kalman filter with Joseph-form covariance update. Operates on a `discrete_state_space` model (A, B, C, D matrices) and provides optimal state estimation for linear systems with Gaussian noise. The Joseph form ensures numerical stability of the covariance update even with finite-precision arithmetic.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::kalman_filter<Scalar, NX, NU, NY>` | `#include <ctrlpp/estimation/kalman.h>` |
| `ctrlpp::kalman_filter<Scalar, NX, NU, NY>` | `#include <ctrlpp/kalman.h>` (convenience) |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
| `NX` | `std::size_t` | State dimension |
| `NU` | `std::size_t` | Input dimension |
| `NY` | `std::size_t` | Measurement dimension |

## Type Aliases

```cpp
using state_vector_t    = Eigen::Matrix<Scalar, nx, 1>;
using input_vector_t    = Eigen::Matrix<Scalar, nu, 1>;
using output_vector_t   = Eigen::Matrix<Scalar, ny, 1>;
using cov_matrix_t      = Eigen::Matrix<Scalar, nx, nx>;
using meas_cov_matrix_t = Eigen::Matrix<Scalar, ny, ny>;
using system_t          = discrete_state_space<Scalar, NX, NU, NY>;
```

## Config (`kalman_config`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | Process noise covariance |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance |
| `x0` | `Vector<Scalar, NX>` | Zero | Initial state estimate |
| `P0` | `Matrix<Scalar, NX, NX>` | Identity | Initial error covariance |

## Constructor

```cpp
kalman_filter(system_t sys, kalman_config<Scalar, NX, NU, NY> config);
```

Constructs the filter from a discrete state-space model and configuration. Uses C++20 designated initializers for config.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates state and covariance one step forward: x = Ax + Bu, P = APA' + Q.

### update

```cpp
auto update(const output_vector_t& z)
    -> ctrlpp::expected<void, kalman_update_error>;
```

Incorporates a measurement via the Kalman gain. Uses Joseph-form covariance update: P = (I-KC)P(I-KC)' + KRK'. Computes and stores the innovation and NIS.

The step is rejected **before any member is assigned**, so a rejected step leaves the state, the covariance, the innovation and the NIS bitwise unchanged and the caller may retry with the next sample. Each cause is an exact domain condition, not a tuning preference: a non-finite operand makes every downstream product non-finite, so no estimate can be formed at all.

| Condition | Error | Why it is a separate cause |
|---|---|---|
| the carried state estimate is already non-finite | `kalman_update_error::non_finite_state` | The fault is upstream of the measurement, so it is reported ahead of it: a caller told the measurement is bad would replace a working sensor while the real fault sits in the prediction that poisoned the state |
| the carried covariance is already non-finite | `kalman_update_error::non_finite_covariance` | The covariance recursion is driven by the model and by `Q` and `R`, never by the measurement, so the caller fixes a different input |
| the supplied measurement has a non-finite component | `kalman_update_error::non_finite_measurement` | The gain carries it into the state, which is the filter's carried memory, so one such sample destroys the estimate permanently |

`predict` is deliberately not fallible. Its input is a command the caller already owns and the plant already took, so refusing it would leave the filter with no propagation for a step that happened. A prediction that poisons the carried estimate is reported by [`health`](#health) instead, which the next `update` latches.

### state

```cpp
const state_vector_t& state() const;
```

Returns the current state estimate.

### covariance

```cpp
const cov_matrix_t& covariance() const;
```

Returns the current error covariance.

### innovation

```cpp
const output_vector_t& innovation() const;
```

Returns the most recent measurement innovation (z - Cx).

### nis

```cpp
Scalar nis() const;
```

Returns the Normalized Innovation Squared from the last update (innovation^T S^-1 innovation), chi-square distributed with dof = NY under a consistent filter.

### health

```cpp
kalman_health health() const;
```

Returns the persistent state-health status, one of `kalman_health::ok` or `kalman_health::non_finite_estimate`. It answers a question a per-call result cannot, because the question outlives the call: whether the carried estimate is still degraded from a step several samples ago. The status starts at `ok` and latches to `non_finite_estimate` the first time a step finds the carried state or covariance already non-finite, which is how a poisoned `predict` becomes visible. A **rejected measurement does not set it**: the rejection mutates nothing, so it leaves the filter healthy. The query carries no discard warning; asking it is optional.

### is_steady_state

```cpp
bool is_steady_state(Scalar tol = Scalar{1e-10}) const;
```

Returns true if the covariance has converged (relative change below tolerance).

### reset_covariance

```cpp
void reset_covariance(const cov_matrix_t& P0);
```

Resets the covariance to a new initial value.

### set_model

```cpp
void set_model(system_t sys);
```

Replaces the state-space model (for adaptive or switched systems).

### set_noise

```cpp
void set_noise(cov_matrix_t Q, meas_cov_matrix_t R);
```

Updates the process and measurement noise covariances.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'true', '' using 1:3 with lines title 'estimate'"

#include <ctrlpp/estimation/kalman.h>
#include <ctrlpp/model/state_space.h>

#include <Eigen/Dense>

#include <iostream>
#include <random>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2, NU = 1, NY = 1;

    // Constant-velocity model: x = [position, velocity]
    ctrlpp::discrete_state_space<Scalar, NX, NU, NY> sys{};
    double dt = 0.1;
    sys.A << 1.0, dt, 0.0, 1.0;
    sys.B << 0.5 * dt * dt, dt;
    sys.C << 1.0, 0.0;
    sys.D << 0.0;

    ctrlpp::kalman_config<Scalar, NX, NU, NY> cfg{
        .Q = (Eigen::Matrix2d() << 0.01, 0.0, 0.0, 0.1).finished(),
        .R = (Eigen::Matrix<Scalar, 1, 1>() << 1.0).finished(),
        .x0 = Eigen::Vector2d::Zero(),
        .P0 = Eigen::Matrix2d::Identity() * 10.0
    };

    ctrlpp::kalman_filter<Scalar, NX, NU, NY> kf(sys, cfg);

    std::mt19937 rng(42);
    std::normal_distribution<> noise(0.0, 1.0);

    Eigen::Vector2d x_true;
    x_true << 0.0, 1.0;  // start at origin, moving at 1 m/s

    for (int k = 0; k < 50; ++k) {
        Eigen::Matrix<Scalar, 1, 1> u = Eigen::Matrix<Scalar, 1, 1>::Zero();
        x_true = sys.A * x_true + sys.B * u;
        Eigen::Matrix<Scalar, 1, 1> z;
        z << x_true(0) + noise(rng);

        kf.predict(u);
        if(!kf.update(z))
        {
            std::cerr << "Kalman filter rejected the measurement at step " << k << "\n";
            return 1;
        }

        auto est = kf.state();
        std::cout << k * dt << "," << x_true(0) << "," << est(0) << "\n";
    }
}
```

## See Also

- [luenberger](luenberger.md)<br/> fixed-gain observer alternative
- [ekf](ekf.md)<br/> nonlinear extension via linearization
- [observer-policy](observer-policy.md)<br/> concept satisfied by this type
- [guides/intro/your-first-estimator](../../guides/intro/your-first-estimator.md)<br/> introductory Kalman filter tutorial
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> composing observers with controllers
- [background/kalman](../../background/kalman.md)<br/> mathematical derivation
