# ekf

Extended Kalman Filter with compile-time dispatch between analytical and numerical Jacobians. The EKF linearises nonlinear dynamics and measurement models around the current state estimate at each step. If the dynamics or measurement model satisfies the `differentiable_dynamics` or `differentiable_measurement` concept (provides a `jacobian_x` / `jacobian` method), the analytical Jacobian is used at zero overhead. Otherwise, numerical central differences are computed automatically.

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

## Config

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
struct ekf_config
{
    Matrix<Scalar, NX, NX> Q;             // process noise covariance (default: identity)
    Matrix<Scalar, NY, NY> R;             // measurement noise covariance (default: identity)
    Vector<Scalar, NX>     x0;            // initial state estimate (default: zero)
    Matrix<Scalar, NX, NX> P0;            // initial covariance (default: identity)
    Scalar numerical_eps;                  // perturbation for numerical Jacobians (default: sqrt(eps))
};
```

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

Propagates state through the dynamics model and linearises to propagate covariance. If `Dynamics` satisfies `differentiable_dynamics`, calls `dynamics.jacobian_x(x, u)` for the state Jacobian F. Otherwise, computes F via numerical central differences.

### update

```cpp
void update(const output_vector_t& z);
```

Incorporates a measurement. Computes the measurement Jacobian H (analytically if `differentiable_measurement` is satisfied, numerically otherwise), then performs the standard Kalman gain computation with Joseph-form covariance update.

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

### nees

```cpp
Scalar nees() const;
```

Normalised Estimation Error Squared from the last update.

## Jacobian Dispatch

The EKF selects Jacobian computation at compile time via `if constexpr`:

```cpp
// Analytical (zero overhead) -- if model provides jacobian_x:
if constexpr (differentiable_dynamics<Dynamics, Scalar, NX, NU>)
    F = dynamics.jacobian_x(x, u);
// Numerical (automatic fallback):
else
    F = numerical_jacobian_x(dynamics, x, u, eps);
```

To provide analytical Jacobians, add a `jacobian_x(x, u)` method to your dynamics callable and/or a `jacobian(x)` method to your measurement callable. See `model/differentiable_dynamics.h` and `model/differentiable_measurement.h` for the concept definitions.

## Usage Example

```cpp
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
        filter.update(z);

        auto est = filter.state();
        std::cout << k * dt << "," << x_true(0) << "," << est(0) << "\n";
    }
}
```

## See Also

- [ukf](ukf.md) -- sigma-point alternative avoiding explicit Jacobians
- [kalman](kalman.md) -- linear Kalman filter for LTI systems
- [observer-policy](observer-policy.md) -- concept satisfied by this type
- [reference/ekf-ukf-theory](../reference/ekf-ukf-theory.md) -- EKF derivation and comparison
