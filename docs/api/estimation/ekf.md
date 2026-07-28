# ekf

Extended Kalman Filter with compile-time dispatch between analytical and numerical Jacobians. The EKF linearizes nonlinear dynamics and measurement models around the current state estimate at each step. If the dynamics or measurement model satisfies the `differentiable_dynamics` or `differentiable_measurement` concept (provides a `jacobian_x` / `jacobian` method), the analytical Jacobian is used at zero overhead. Otherwise, numerical central differences are computed automatically.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::ekf<Scalar, NX, NU, NY, Dynamics, Measurement>` | `#include <ctrlpp/estimation/ekf.h>` |

No convenience header exists for this type. Use the categorical path.

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
| `NX` | `std::size_t` | State dimension |
| `NU` | `std::size_t` | Input dimension |
| `NY` | `std::size_t` | Measurement dimension |
| `Dynamics` | satisfies `dynamics_model<Scalar, NX, NU>` | Callable: `(Vector<NX>, Vector<NU>) -> Vector<NX>` |
| `Measurement` | satisfies `measurement_model<Scalar, NX, NY>` | Callable: `(Vector<NX>) -> Vector<NY>` |

## Type Aliases

```cpp
using state_vector_t    = Vector<Scalar, NX>;
using input_vector_t    = Vector<Scalar, NU>;
using output_vector_t   = Vector<Scalar, NY>;
using cov_matrix_t      = Matrix<Scalar, NX, NX>;
using meas_cov_matrix_t = Matrix<Scalar, NY, NY>;
```

## Config (`ekf_config`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | Process noise covariance |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance |
| `x0` | `Vector<Scalar, NX>` | Zero | Initial state estimate |
| `P0` | `Matrix<Scalar, NX, NX>` | Identity | Initial error covariance |
| `numerical_eps` | `Scalar` | `sqrt(eps)` | Perturbation step for numerical Jacobians (only used when analytical Jacobians are not provided) |

## Constructor

```cpp
ekf(Dynamics dynamics, Measurement measurement, ekf_config<Scalar, NX, NU, NY> config);
```

CTAD deduction guide available: template parameters are deduced from the config type.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates state through the dynamics model and linearizes to propagate covariance. If `Dynamics` satisfies `differentiable_dynamics`, calls `dynamics.jacobian_x(x, u)` for the state Jacobian F. Otherwise, computes F via numerical central differences.

### update

```cpp
auto update(const output_vector_t& z)
    -> ctrlpp::expected<void, ekf_update_error>;
```

Incorporates a measurement. Computes the measurement Jacobian H (analytically if `differentiable_measurement` is satisfied, numerically otherwise), then performs the standard Kalman gain computation with Joseph-form covariance update.

The step is rejected **before any member is assigned**, so a rejected step leaves the state, the covariance, the innovation and the NIS bitwise unchanged and the caller may retry with the next sample. Each cause is an exact domain condition, not a tuning preference.

| Condition | Error | Why it is a separate cause |
|---|---|---|
| the carried state estimate is already non-finite | `ekf_update_error::non_finite_state` | The measurement Jacobian is evaluated **at** the carried state, so a poisoned state makes the linearization meaningless before the measurement is used at all |
| the carried covariance is already non-finite | `ekf_update_error::non_finite_covariance` | The covariance recursion is driven by the linearized dynamics and by `Q` and `R`, never by the measurement |
| the supplied measurement has a non-finite component | `ekf_update_error::non_finite_measurement` | The gain carries it into the state, which is the filter's carried memory, so one such sample destroys the estimate permanently |

`predict` is deliberately not fallible. Its input is a command the caller already owns and the plant already took, so refusing it would leave the filter with no propagation for a step that happened. A prediction that poisons the carried estimate is reported by [`health`](#health) instead, which the next `update` latches.

### state

```cpp
const state_vector_t& state() const;
```

### covariance

```cpp
const cov_matrix_t& covariance() const;
```

### innovation

```cpp
const output_vector_t& innovation() const;
```

### nis

```cpp
Scalar nis() const;
```

Normalized Innovation Squared from the last update (innovation^T S^-1 innovation), chi-square distributed with dof = NY under a consistent filter.

### health

```cpp
ekf_health health() const;
```

Returns the persistent state-health status, one of `ekf_health::ok` or `ekf_health::non_finite_estimate`. It answers a question a per-call result cannot, because the question outlives the call: whether the carried estimate is still degraded from a step several samples ago. The status starts at `ok` and latches to `non_finite_estimate` the first time a step finds the carried state or covariance already non-finite, which is how a poisoned `predict` or a dynamics model that returned a non-finite state becomes visible. A **rejected measurement does not set it**: the rejection mutates nothing, so it leaves the filter healthy. The query carries no discard warning; asking it is optional.

## Jacobian Dispatch

The EKF selects Jacobian computation at compile time via `if constexpr`:

```cpp
// Analytical (zero overhead)<br/> if model provides jacobian_x:
if constexpr (differentiable_dynamics<Dynamics, Scalar, NX, NU>)
    F = dynamics.jacobian_x(x, u);
// Numerical (automatic fallback):
else
    F = numerical_jacobian_x(dynamics, x, u, eps);
```

To provide analytical Jacobians, add a `jacobian_x(x, u)` method to your dynamics callable and/or a `jacobian(x)` method to your measurement callable. See `model/differentiable_dynamics.h` and `model/differentiable_measurement.h` for the concept definitions.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'true', '' using 1:3 with lines title 'estimate'"

#include <ctrlpp/estimation/ekf.h>

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <random>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2, NU = 1, NY = 1;
    constexpr Scalar dt = 0.01;

    // Simple pendulum: x = [theta, omega]
    auto dynamics = [dt](const ctrlpp::Vector<Scalar, NX>& x,
                         const ctrlpp::Vector<Scalar, NU>& u) -> ctrlpp::Vector<Scalar, NX> {
        ctrlpp::Vector<Scalar, NX> xn;
        xn(0) = x(0) + dt * x(1);
        xn(1) = x(1) + dt * (-9.81 * std::sin(x(0)) + u(0));
        return xn;
    };

    auto measurement = [](const ctrlpp::Vector<Scalar, NX>& x) -> ctrlpp::Vector<Scalar, NY> {
        ctrlpp::Vector<Scalar, NY> z;
        z(0) = x(0);  // observe angle
        return z;
    };

    ctrlpp::ekf_config<Scalar, NX, NU, NY> cfg{
        .Q = (Eigen::Matrix2d() << 0.001, 0.0, 0.0, 0.01).finished(),
        .R = (Eigen::Matrix<Scalar, 1, 1>() << 0.1).finished(),
        .x0 = Eigen::Vector2d::Zero(),
        .P0 = Eigen::Matrix2d::Identity()
    };

    ctrlpp::ekf filter(dynamics, measurement, cfg);

    std::mt19937 rng(42);
    std::normal_distribution<> noise(0.0, std::sqrt(0.1));

    ctrlpp::Vector<Scalar, NX> x_true;
    x_true << 0.5, 0.0;

    for (int k = 0; k < 500; ++k) {
        ctrlpp::Vector<Scalar, NU> u = ctrlpp::Vector<Scalar, NU>::Zero();
        x_true = dynamics(x_true, u);

        ctrlpp::Vector<Scalar, NY> z;
        z(0) = x_true(0) + noise(rng);

        filter.predict(u);
        if(!filter.update(z))
        {
            std::cerr << "EKF rejected the measurement at step " << k << "\n";
            return 1;
        }

        auto est = filter.state();
        std::cout << k * dt << "," << x_true(0) << "," << est(0) << "\n";
    }
}
```

## See Also

- [ukf](ukf.md)<br/> sigma-point alternative avoiding explicit Jacobians
- [kalman](kalman.md)<br/> linear Kalman filter for LTI systems
- [observer-policy](observer-policy.md)<br/> concept satisfied by this type
- [background/ekf-ukf](../../background/ekf-ukf.md)<br/> EKF derivation and comparison
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> composing observers with controllers
