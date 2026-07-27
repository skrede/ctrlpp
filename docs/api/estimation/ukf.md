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

## Construction

```cpp
ukf(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config);

ukf(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config,
    Strategy strategy);

[[nodiscard]] static auto try_create(Dynamics dynamics, Measurement measurement,
                                     ukf_config<Scalar, NX, NU, NY> config,
                                     typename Strategy::options_t strategy_options)
    -> ctrlpp::expected<ukf, filter_error>;

ukf(Dynamics dynamics, Measurement measurement, ukf_config<Scalar, NX, NU, NY> config,
    typename Strategy::options_t strategy_options);
```

The first two forms cannot fail: the default-strategy form uses the strategy's own in-domain defaults, and the pre-built-strategy form receives a strategy that was validated where it was constructed.

The options-aggregate form is the one that builds the strategy inside the filter, so it is the only path on which a sigma-point parameter set can be out of domain. `try_create` is that path's primary API: it forwards the strategy's rejection verbatim as a `filter_error` (from `<ctrlpp/estimation/estimation_types.h>`). The matching constructor is the exception-gated convenience wrapper: with `merwe_sigma_points` it throws in an exceptions-enabled build and does not compile on an exception-free build, where `try_create` is the construction path to use. As a static member of a class template, `try_create` requires explicit template arguments, so name the filter type first:

```cpp
using filter_t = ctrlpp::ukf<double, NX, NU, NY, Dynamics, Measurement>;

auto filter = filter_t::try_create(Dynamics{}, Measurement{}, cfg,
                                   ctrlpp::merwe_options<double>{.alpha = 1e-3, .beta = 2.0, .kappa = 0.0});
if(!filter)
    return filter.error();
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

Generates sigma points, transforms through measurement model, computes innovation covariance S and cross-covariance Pxz, then applies the Kalman gain correction. The covariance update uses the algebraically complete minimum mean-square-error reduction `P = P - K*S*K^T`, with `S = Pzz + R` and `K = Pxz*S^{-1}`. Since `K*S*K^T = K*Pxz^T`, this term is exactly the uncertainty the measurement removes, and no extra `K*R*K^T` term is added. The sigma points feeding this update are built from a permutation-correct covariance square root, so the reduction stays symmetric positive semidefinite.

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

Returns the Normalized Innovation Squared from the last update, `innovation^T S^-1 innovation`, chi-square distributed with dof = NY under a consistent filter. Useful for online consistency monitoring.

### health

```cpp
ukf_health health() const;
```

Returns the filter-health status, one of `ukf_health::ok` or `ukf_health::covariance_repaired`. The status starts at `ok` and latches to `covariance_repaired` the first time a non-positive-definite covariance had to be repaired to the nearest symmetric positive definite matrix during sigma-point generation, signaling that the estimate has entered a numerically degraded regime.

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

#### Parameter domain

The strategy is constructed through a fallible factory, because two of its three parameters have an exact admissible domain:

```cpp
[[nodiscard]] static auto try_create(options_t opts = options_t{})
    -> ctrlpp::expected<merwe_sigma_points, filter_error>;
```

The scaling term is `lambda = alpha^2 (n + kappa) - n`, so the weight denominator `n + lambda` and the squared sigma-point offset scale `gamma^2 = n + lambda` both collapse to `alpha^2 (n + kappa)`. The divisor and the radicand are the same expression, which fixes the domain exactly. There is no tolerance and no fitted constant.

| Rejection | Enumerator | Why |
|-----------|------------|-----|
| `alpha` is NaN, infinite, zero, or negative | `filter_error::non_positive_sigma_spread` | It divides the weight denominator. A zero spread makes the weights infinite; a negative one makes them finite but wrong, since only `alpha^2` is used |
| `kappa` is NaN or infinite, or `n + kappa <= 0` | `filter_error::non_positive_scaling_radicand` | The sum sits under the square root that scales the offsets and inside the same denominator, so a non-positive sum yields non-finite offsets or non-finite weights |

`beta` is not validated: it enters only the additive prior-kurtosis term of the first covariance weight and carries no domain restriction of this kind.

Default construction (`merwe_sigma_points<Scalar, NX>{}`) cannot fail: the default `alpha` is finite and positive, and with the default zero `kappa` the sum `n + kappa` reduces to `NX`, which is required to be positive. The options constructor is the exception-gated wrapper over `try_create` and is available only in an exceptions-enabled build. `so3_merwe_sigma_points` forwards the same check unchanged, since it lifts this strategy's tangent-space points onto SO(3) and inherits its weights.

Choosing good default values across dimension, scale, and scalar tier is a separate question from admissibility, and these checks do not address it.

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

- [ekf](ekf.md)<br/> linearization-based alternative
- [manifold-ukf](manifold-ukf.md)<br/> UKF on SO(3) manifold for attitude estimation
- [observer-policy](observer-policy.md)<br/> concept satisfied by this type
- [background/ekf-ukf](../../background/ekf-ukf.md)<br/> unscented transform derivation
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> composing observers with controllers
