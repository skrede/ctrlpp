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

## Construction

```cpp
static auto create(Dynamics dynamics, Measurement measurement,
                   ekf_config<Scalar, NX, NU, NY> config)
    -> ctrlpp::expected<ekf, filter_error>;
```

`create` is the only construction path; there is no public constructor. A
factory that sat beside one would validate nothing, since any caller could
bypass it. The CTAD deduction guide was removed with the constructor it named,
so the template arguments are written out:

```cpp
auto filter = ctrlpp::ekf<double, NX, NU, NY, Dynamics, Measurement>::create(
    Dynamics{}, Measurement{}, cfg);
if(!filter)
    return filter.error();
```

Rejections, checked in the order the config aggregate declares its fields (the
fields are independent, so the order is declaration order rather than a claim
about causality). Each matrix is checked for finiteness, then for symmetry, then for definiteness:

| Enumerator | Condition | What it prevents |
|------------|-----------|------------------|
| `filter_error::non_finite_process_noise` | `Q` has a non-finite entry | `Q` is added to the propagated covariance, so an infinite entry gives `P = Inf`, a gain formed from `Inf/Inf` and therefore `NaN`, and a non-finite estimate at the first step |
| `filter_error::asymmetric_process_noise` | `Q`'s asymmetry exceeds the rounding bound | a covariance is symmetric by definition; triangles that differ by more than rounding describe two different matrices |
| `filter_error::indefinite_process_noise` | `Q` has an eigenvalue below `-floor` | a matrix that is not positive semidefinite is not a covariance, so the predicted covariance it is added to stops being one |
| `filter_error::non_finite_measurement_noise` | `R` has a non-finite entry | `R` is added to the innovation covariance the gain solve is posed against, so that solve is meaningless for every measurement |
| `filter_error::asymmetric_measurement_noise` | `R`'s asymmetry exceeds the rounding bound | the innovation covariance the gain is posed against would carry it |
| `filter_error::indefinite_measurement_noise` | `R` has an eigenvalue below `-floor` | the innovation covariance the gain is posed against stops being a covariance |
| `filter_error::non_finite_initial_state` | `x0` has a non-finite component | `x0` seeds the carried estimate, which is the filter's memory, so every later estimate is formed from it |
| `filter_error::non_finite_initial_covariance` | `P0` has a non-finite entry | `P0` seeds the covariance recursion, which has no mechanism that returns a non-finite covariance to a finite one |
| `filter_error::asymmetric_initial_covariance` | `P0`'s asymmetry exceeds the rounding bound | there is no Gaussian prior with an asymmetric covariance |
| `filter_error::indefinite_initial_covariance` | `P0` has an eigenvalue below `-floor` | there is no Gaussian prior with an indefinite covariance |

The bound is `max(N eps max|M_ij|, N denorm_min)` for both checks. As a
definiteness floor it is the order of the symmetric eigensolver's backward
error, so an eigenvalue within it of zero is indistinguishable from zero; the
absolute term keeps it from underflowing to zero for a subnormal-scale matrix.
An asymmetry `max|M_ij - M_ji|` within it is rounding drift: the matrix is
accepted, the filter stores its symmetrized form (each differing pair replaced
by its midpoint), and the definiteness verdict is taken on that symmetrized
matrix. This is the same rule every filter in the estimation module applies: a
matrix that is not symmetric positive semidefinite is not a covariance,
whatever the filter does with it. The convention and its count are in
[Rounding drift is accepted, structure is refused](../../guides/patterns/numerical-behavior.md#rounding-drift-is-accepted-structure-is-refused).

A **singular** positive semidefinite matrix is **accepted**: a point-mass
`P0 = 0`, a rank-one `Q` and a noiseless `R = 0` are legitimate configurations,
because the recursion only adds and multiplies them and inverts the innovation
covariance, never `R` alone. An ill-conditioned configuration, with entries many
orders of magnitude apart, is accepted as well: conditioning is a
numerical-behavior question, not a domain violation.

`numerical_eps` is **not** validated. It is the finite-difference step used only
when the dynamics or measurement model is not analytically differentiable, and
its own domain condition (finite and strictly positive, since the
central-difference stencil divides by it) is separately owned work.

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
