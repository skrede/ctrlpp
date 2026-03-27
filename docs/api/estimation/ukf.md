# ukf

Unscented Kalman Filter with swappable sigma point strategies and configurable gain decomposition. The UKF avoids explicit Jacobian computation by propagating a deterministic set of sigma points through the nonlinear dynamics and measurement models, then recovering mean and covariance from the transformed points. All sigma point storage uses `std::array` to avoid heap allocation.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::ukf<Scalar, NX, NU, NY, Dynamics, Measurement, Strategy>` | `#include <ctrlpp/estimation/ukf.h>` |

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
| `Strategy` | satisfies `sigma_point_strategy<Scalar, NX>` | Sigma point generator (default: `merwe_sigma_points<Scalar, NX>`) |

## Type Aliases

```cpp
using state_vector_t    = Vector<Scalar, NX>;
using input_vector_t    = Vector<Scalar, NU>;
using output_vector_t   = Vector<Scalar, NY>;
using cov_matrix_t      = Matrix<Scalar, NX, NX>;
using meas_cov_matrix_t = Matrix<Scalar, NY, NY>;
```

## Config (`ukf_config`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | Process noise covariance |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance |
| `x0` | `Vector<Scalar, NX>` | Zero | Initial state estimate |
| `P0` | `Matrix<Scalar, NX, NX>` | Identity | Initial error covariance |
| `decomposition` | `gain_decomposition` | `ldlt` | Kalman gain decomposition method (`ldlt` or `qr`) |

## Constructors

```cpp
ukf(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config);

ukf(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config,
    typename Strategy::options_t strategy_options);

ukf(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config,
    Strategy strategy);
```

CTAD deduction guide available: deduces to `merwe_sigma_points` as default strategy.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Generates sigma points from current (x, P), propagates them through dynamics, and recovers predicted mean and covariance.

### update

```cpp
void update(const output_vector_t& z);
```

Generates sigma points, transforms through measurement model, computes innovation covariance S and cross-covariance Pxz, then applies the Kalman gain correction.

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

## Supporting Types

### merwe_sigma_points

The default sigma point strategy. Generates 2N+1 sigma points using the scaled unscented transform.

```cpp
template <typename Scalar>
struct merwe_options
{
    Scalar alpha{Scalar{1e-3}};  // spread around mean (small -> tight)
    Scalar beta{Scalar{2}};     // prior knowledge (2 optimal for Gaussian)
    Scalar kappa{Scalar{0}};    // secondary scaling
};
```

Header: `#include <ctrlpp/estimation/sigma_points/merwe_sigma_points.h>`

### julier_sigma_points

Alternative strategy with a single kappa parameter. Generates 2N+1 sigma points.

Header: `#include <ctrlpp/estimation/sigma_points/julier_sigma_points.h>`

### gain_decomposition

```cpp
enum class gain_decomposition { ldlt, qr };
```

Selects the decomposition for K = Pxz * S^{-1}. LDLT is faster for well-conditioned S; QR is more robust for ill-conditioned problems.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'true', '' using 1:3 with lines title 'estimate'"

#include <ctrlpp/estimation/ukf.h>

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <random>

int main()
{
    using Scalar = double;
    constexpr std::size_t NX = 2, NU = 1, NY = 1;
    constexpr Scalar dt = 0.01;

    auto dynamics = [dt](const ctrlpp::Vector<Scalar, NX>& x,
                         const ctrlpp::Vector<Scalar, NU>& u) -> ctrlpp::Vector<Scalar, NX> {
        ctrlpp::Vector<Scalar, NX> xn;
        xn(0) = x(0) + dt * x(1);
        xn(1) = x(1) + dt * (-9.81 * std::sin(x(0)) + u(0));
        return xn;
    };

    auto measurement = [](const ctrlpp::Vector<Scalar, NX>& x) -> ctrlpp::Vector<Scalar, NY> {
        ctrlpp::Vector<Scalar, NY> z;
        z(0) = x(0);
        return z;
    };

    ctrlpp::ukf_config<Scalar, NX, NU, NY> cfg{
        .Q = (Eigen::Matrix2d() << 0.001, 0.0, 0.0, 0.01).finished(),
        .R = (Eigen::Matrix<Scalar, 1, 1>() << 0.1).finished(),
        .x0 = Eigen::Vector2d::Zero(),
        .P0 = Eigen::Matrix2d::Identity()
    };

    ctrlpp::ukf filter(dynamics, measurement, cfg);

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

- [ekf](ekf.md) -- linearisation-based alternative
- [manifold-ukf](manifold-ukf.md) -- UKF on SO(3) manifold for attitude estimation
- [observer-policy](observer-policy.md) -- concept satisfied by this type
- [background/ekf-ukf](../../background/ekf-ukf.md) -- unscented transform derivation
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md) -- composing observers with controllers
