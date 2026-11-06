# ukf

Unscented Kalman Filter with swappable sigma-point strategies and a selectable gain decomposition. The UKF needs no Jacobian: it passes a deterministic set of sigma points through the nonlinear dynamics and measurement models and recovers the mean and covariance from the transformed points. All sigma-point storage is `std::array`, so no step allocates.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::ukf<Scalar, NX, NU, NY, Dynamics, Measurement, Strategy>` | `#include <ctrlpp/estimation/ukf.h>` |
| `ctrlpp::ukf<Scalar, NX, NU, NY, Dynamics, Measurement, Strategy>` | `#include <ctrlpp/ukf.h>` (convenience) |
| `ctrlpp::ukf_config`, `ctrlpp::gain_decomposition`, `ctrlpp::ukf_update_error`, `ctrlpp::ukf_health` | `#include <ctrlpp/estimation/ukf_types.h>` (included by `ukf.h`) |
| `ctrlpp::sigma_point_strategy`, `ctrlpp::sigma_result` | `#include <ctrlpp/estimation/sigma_points/sigma_point_strategy.h>` (included by `ukf.h`) |
| `ctrlpp::merwe_sigma_points`, `ctrlpp::merwe_options` | `#include <ctrlpp/estimation/sigma_points/merwe_sigma_points.h>` (included by `ukf.h`) |
| `ctrlpp::julier_sigma_points`, `ctrlpp::julier_options` | `#include <ctrlpp/estimation/sigma_points/julier_sigma_points.h>` (not included by `ukf.h`) |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type (e.g. `double`, `float`) |
| `NX` | `std::size_t` | State dimension |
| `NU` | `std::size_t` | Input dimension |
| `NY` | `std::size_t` | Measurement dimension |
| `Dynamics` | satisfies `dynamics_model<Dynamics, Scalar, NX, NU>` | Callable: `(Vector<NX>, Vector<NU>) -> Vector<NX>` |
| `Measurement` | satisfies `measurement_model<Measurement, Scalar, NX, NY>` | Callable: `(Vector<NX>) -> Vector<NY>` |
| `Strategy` | satisfies `sigma_point_strategy<Strategy, Scalar, NX>` | Sigma-point generator (default: `merwe_sigma_points<Scalar, NX>`) |

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
| `decomposition` | `gain_decomposition` | `ldlt` | Factorization the gain `K = Pxz S^-1` is solved with (see [gain_decomposition](#gain_decomposition)) |

## Construction

```cpp
static ctrlpp::expected<ukf, filter_error> create(Dynamics dynamics, Measurement measurement,
                                                  ukf_config<Scalar, NX, NU, NY> config,
                                                  Strategy strategy = Strategy{});

static ctrlpp::expected<ukf, filter_error> create(Dynamics dynamics, Measurement measurement,
                                                  ukf_config<Scalar, NX, NU, NY> config,
                                                  typename Strategy::options_t strategy_options);
```

`create` is the only construction path; there is no public constructor. The
header declares a deduction guide from `(dynamics, measurement, config)`, but no
public constructor matches it, so class template argument deduction does not
produce a filter: the template arguments are written out (with a lambda model,
as `decltype(model)`). The sigma-point strategies are built through
`try_create`, beside a default constructor that cannot fail.

The first overload takes a strategy already built, so the strategy was validated
where it was built; omit it and the strategy's own in-domain defaults are used.
Both shipped strategies, `merwe_sigma_points` and `julier_sigma_points`, are
built through their own `try_create`, so either works with this overload and
with the next one.

The second overload builds the strategy inside the filter, so it is the only
path on which a sigma-point parameter set can be out of domain. It forwards the
strategy's rejection verbatim, **before** validating the configuration: a
strategy that cannot be built leaves nothing for the configuration to be a
configuration of. It requires the strategy to expose a `try_create`.

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

When the unpivoted Cholesky factorization the sigma points are built from fails
on a singular `P0`, the first sigma-point generation repairs it and `health()`
reports `covariance_repaired`; that is a reported numerical event, not a
refusal.

As a static member of a class template it requires explicit template arguments,
so name the filter type first:

```cpp
using filter_t = ctrlpp::ukf<double, NX, NU, NY, Dynamics, Measurement>;

auto filter = filter_t::create(Dynamics{}, Measurement{}, cfg,
                               ctrlpp::merwe_options<double>{.alpha = 1e-3, .beta = 2.0, .kappa = 0.0});
if(!filter)
    return filter.error();
```

`merwe_sigma_points` is the default strategy template argument. A filter on
`julier_sigma_points` names it as the seventh template argument and includes
its header.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Generates sigma points from the current `(x, P)`, propagates each through the dynamics with the same `u`, and recovers the predicted mean `sum Wm_i f(X_i, u)` and covariance `sum Wc_i (f(X_i, u) - x)(f(X_i, u) - x)' + Q`, symmetrized. A sigma-point generation that had to repair `P` latches `covariance_repaired` from `ok`.

`predict` then judges what it produced and latches [`health`](#health): `non_finite_estimate` when the predicted state or covariance has a non-finite entry, otherwise `covariance_repaired` (from `ok`) when the predicted covariance is not positive definite. Definiteness is decided by the unpivoted Cholesky factorization of the predicted covariance, which is the factorization the Merwe and Julier generations attempt first on that same matrix, so `predict` reports exactly the repair the next sigma-point generation will make. A committed [`update`](#update) judges the covariance it leaves the same way, so `covariance()` never reports an indefinite matrix beside `health() == ok`, after either call. A strategy that factors otherwise still gets this verdict on definiteness. The cost is one `NX x NX` Cholesky factorization per step and its factor in `predict`'s frame; the arithmetic of the prediction itself is unchanged.

### update

```cpp
ctrlpp::expected<void, ukf_update_error> update(const output_vector_t& z);
```

A refused step leaves the state, the covariance, the innovation and the NIS bitwise unchanged, so the caller may retry with the next sample. The first three causes are checked before anything is computed; the fourth is checked on the computed result, which is then discarded together with any health change the step made. Each cause is an exact domain condition, not a tuning preference.

| Condition | Error | Why it is a separate cause |
|---|---|---|
| the carried state estimate is already non-finite | `ukf_update_error::non_finite_state` | The sigma points are generated **around** the carried state, so every one of them is non-finite before the measurement is used at all |
| the carried covariance is already non-finite | `ukf_update_error::non_finite_covariance` | The covariance is what the sigma-point spread factors, and it is driven by the model and by `Q` and `R`, never by the measurement |
| the supplied measurement has a non-finite component | `ukf_update_error::non_finite_measurement` | The gain carries it into the state, which is the filter's carried memory, so one such sample destroys the estimate permanently |
| the updated state, covariance, innovation or NIS is non-finite although the carried estimate and the measurement were finite | `ukf_update_error::non_finite_result` | The step's own arithmetic failed: the measurement model returned a non-finite value at a sigma point, a weighted sum overflowed, or the gain solve divided by zero |

`predict` is deliberately not fallible. Its input is a command the caller already owns and the plant already took, so refusing it would leave the filter with no propagation for a step that happened. A prediction that poisons the carried estimate, or predicts a covariance that is not positive definite, is reported by [`health`](#health) instead, which `predict` latches itself; the next `update` then refuses to start from a non-finite estimate. An update refused with `non_finite_result` restores the health it began with, so the rung `predict` latched survives the refusal.

On a step that runs, `update` generates a fresh sigma set from the predicted `(x, P)` (it does not reuse the points `predict` propagated), passes each point through the measurement model, and forms `z_pred = sum Wm_i h(X_i)`, `Pzz = sum Wc_i (h(X_i) - z_pred)(h(X_i) - z_pred)'`, `S = Pzz + R` and `Pxz = sum Wc_i (X_i - x)(h(X_i) - z_pred)'`. The gain is `K = Pxz S^-1`, the state becomes `x + K(z - z_pred)` and the covariance `P - K*S*K^T`, symmetrized. With `K = Pxz S^-1`, `K*S*K^T = K*Pxz^T`, so the reduction is the uncertainty the measurement removes and no `K*R*K^T` term is added; the code relies on nothing more than that identity. The reduction is taken from the covariance the sigma set represents, which is `S S^T` for the square root `S` the points were built from. When the carried covariance was not positive definite and had to be repaired, that is the repaired matrix, so a repaired update carries the repaired covariance less the gain term, which with a nonnegative center weight is positive semidefinite up to rounding, and `health()` reports `covariance_repaired`.

A committed update then judges the covariance it leaves by the same unpivoted Cholesky factorization `predict` uses, and latches `covariance_repaired` from `ok` when it fails. The reduction can leave an indefinite posterior from a positive definite prior when the center covariance weight is negative, as with Julier at a negative `kappa` (Julier's recommended `3 - n` is negative above three states) and with Merwe at its default options (see [Limitations](#limitations)). With a nonnegative center weight the joint sigma-point covariance of state and measurement is positive semidefinite, and so is the posterior; a negative one removes that guarantee, so with a unit prior and a measurement whose linear part reads one state, the posterior variance along that state is `1 - 1/S`, negative once the sigma-point `Pzz` pulls `S` below 1. The judgment costs one `NX x NX` Cholesky factorization per committed update and changes no value: the state, covariance, innovation and NIS are the ones the update computed. An update refused with `non_finite_result` still restores the health it began with.

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

Returns the normalized innovation squared of the last committed update, `innovation' S^-1 innovation`, with `S` the sigma-point innovation covariance above solved by a column-pivoting QR factorization whatever `decomposition` says. For a consistent filter it is approximately chi-square distributed with `NY` degrees of freedom. It is zero before the first update.

### health

```cpp
ukf_health health() const;
```

Returns the persistent state-health status, one of `ukf_health::ok`, `ukf_health::covariance_repaired` or `ukf_health::non_finite_estimate`. It answers a question a per-call result cannot, because the question outlives the call: whether the carried estimate is still degraded from a step several samples ago. The enumerators are ordered by severity and the status latches, never returning to a lower rung on its own. It moves to `covariance_repaired` the first time `predict` or a committed `update` produces a covariance that is not positive definite or a sigma-point generation has to repair one to the nearest symmetric positive definite matrix, signaling a numerically degraded but still usable estimate; it moves to `non_finite_estimate` the first time the carried state or covariance is non-finite, which is not recoverable by any further step. `predict` latches both rungs on what it produced and a committed `update` latches `covariance_repaired` on what it leaves, so a poisoned or indefinite estimate is visible as soon as the call that produced it returns rather than at the next sigma-point generation. An update refused with `non_finite_state` or `non_finite_covariance` latches `non_finite_estimate`; one refused with `non_finite_measurement` leaves the status unchanged, and one refused with `non_finite_result` restores the status it began with. Nothing clears it.

## Supporting Types

### merwe_sigma_points

The default sigma-point strategy. Generates `2 NX + 1` sigma points using the scaled unscented transform.

```cpp
template <typename Scalar>
struct merwe_options
{
    Scalar alpha{static_cast<Scalar>(1e-3)};  // spread around mean (small -> tight)
    Scalar beta{Scalar{2}};                   // prior knowledge (2 optimal for Gaussian)
    Scalar kappa{Scalar{0}};                  // secondary scaling
};
```

Header: `#include <ctrlpp/estimation/sigma_points/merwe_sigma_points.h>`

#### Parameter domain

The strategy is constructed through a fallible factory, because two of its three parameters have an exact admissible domain:

```cpp
static ctrlpp::expected<merwe_sigma_points, filter_error> try_create(options_t opts = options_t{});
```

The scaling term is `lambda = alpha^2 (n + kappa) - n`, so the weight denominator `n + lambda` and the squared sigma-point offset scale `gamma^2 = n + lambda` both collapse to `alpha^2 (n + kappa)`. The divisor and the radicand are the same expression, which fixes the domain exactly. There is no tolerance and no fitted constant.

| Rejection | Enumerator | Why |
|-----------|------------|-----|
| `alpha` is NaN, infinite, zero, or negative | `filter_error::non_positive_sigma_spread` | It divides the weight denominator. A zero spread makes the weights infinite; a negative one makes them finite but wrong, since only `alpha^2` is used |
| `kappa` is NaN or infinite, or `n + kappa <= 0` | `filter_error::non_positive_scaling_radicand` | The sum sits under the square root that scales the offsets and inside the same denominator, so a non-positive sum yields non-finite offsets or non-finite weights |

`beta` is not validated: it enters only the additive prior-kurtosis term of the first covariance weight and carries no domain restriction of this kind.

Default construction (`merwe_sigma_points<Scalar, NX>{}`) cannot fail: the default `alpha` is finite and positive, and with the default zero `kappa` the sum `n + kappa` reduces to `NX`, which is required to be positive. `so3_merwe_sigma_points` forwards the same check unchanged, since it lifts this strategy's tangent-space points onto SO(3) and inherits its weights.

The points are `x` and `x +/- gamma S_i` for the columns `S_i` of the lower Cholesky factor of `P`, `gamma = sqrt(n + lambda)`, with weights `Wm_0 = lambda / (n + lambda)`, `Wc_0 = Wm_0 + 1 - alpha^2 + beta` and `Wm_i = Wc_i = 1 / (2 (n + lambda))`. A `P` the Cholesky factorization refuses is repaired by clamping its eigenvalues up to `N eps max|P_ij|` and the result sets `spd_repaired`. Admissibility is all `try_create` checks; whether the defaults suit a given problem is a separate question, stated under [Limitations](#limitations).

### julier_sigma_points

Alternative strategy with a single kappa parameter. Generates `2 NX + 1` sigma points, `x` and `x +/- sqrt(n + kappa) S_i`, with `Wm_0 = Wc_0 = kappa / (n + kappa)` and `Wm_i = Wc_i = 1 / (2 (n + kappa))`, repairing `P` as the Merwe strategy does.

```cpp
template <typename Scalar>
struct julier_options
{
    Scalar kappa{Scalar{0}};
};
```

Header: `#include <ctrlpp/estimation/sigma_points/julier_sigma_points.h>`

#### Parameter domain

The strategy is constructed through a fallible factory, because kappa has an exact admissible domain:

```cpp
static ctrlpp::expected<julier_sigma_points, filter_error> try_create(options_t opts = options_t{});
```

The sum `n + kappa` is both the radicand of the offset scale `gamma = sqrt(n + kappa)` and the weight denominator, which fixes the domain exactly. There is no tolerance and no fitted constant.

| Rejection | Enumerator | Why |
|-----------|------------|-----|
| `kappa` is NaN or infinite, or `n + kappa <= 0` | `filter_error::non_positive_scaling_radicand` | The sum sits under the square root that scales the offsets and divides every weight, so a non-positive sum yields non-finite offsets or non-finite weights |

A negative `kappa` with `n + kappa > 0` is inside the domain. Julier's own choice `kappa = 3 - n` is such a value for `n > 3`; it gives a negative center weight, which the method permits.

Default construction (`julier_sigma_points<Scalar, NX>{}`) cannot fail: with the default zero `kappa` the sum `n + kappa` reduces to `NX`, which is required to be positive. There is no public constructor from options; a non-default `kappa` goes through `try_create`, or through the options overload of the filter's `create`, which forwards the same rejection.

### Custom strategies

A strategy satisfies `sigma_point_strategy<S, Scalar, NX>` (`ctrlpp/estimation/sigma_points/sigma_point_strategy.h`) by exposing a static `num_points`, a member type `options_t` and a const `generate(x, P)` returning `sigma_result<Scalar, NX, num_points>`. The options overload of `create` additionally calls `S::try_create(options_t)`, so a strategy without one works only with the first overload:

```cpp
template <typename Scalar, std::size_t NX, std::size_t NumPoints>
struct sigma_result
{
    sigma_result();   // spd_repaired = false; the other fields are left for the strategy to fill

    std::array<Vector<Scalar, NX>, NumPoints> points;
    std::array<Scalar, NumPoints> Wm;
    std::array<Scalar, NumPoints> Wc;
    bool spd_repaired;
    Matrix<Scalar, NX, NX> factor;
};
```

`factor` is the matrix `S` with `S S^T` equal to the covariance the points were built from. The filter reads it only when `spd_repaired` is set, as the base of the covariance reduction, so a custom strategy that repairs the covariance and sets the flag must fill `factor` with the repaired square root.

### gain_decomposition

```cpp
enum class gain_decomposition { ldlt, qr };
```

Selects how `K = Pxz S^-1` is solved. `ldlt` uses Eigen's pivoted LDLT factorization of `S`; `qr` uses a column-pivoting Householder QR factorization of `S'`, or a division when `NY = 1`. Neither tests `S` for rank, and the NIS is solved with the QR factorization in both cases.

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
        .R = (Eigen::Matrix<Scalar, 1, 1>() << 0.1).finished()
    };

    using filter_t = ctrlpp::ukf<Scalar, NX, NU, NY, decltype(dynamics), decltype(measurement)>;
    auto created = filter_t::create(dynamics, measurement, cfg);
    if(!created)
    {
        std::cerr << "invalid UKF configuration\n";
        return 1;
    }
    auto& filter = *created;

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
            std::cerr << "UKF rejected the measurement\n";
            return 1;
        }

        auto est = filter.state();
        std::cout << k * dt << "," << x_true(0) << "," << est(0) << "\n";
    }
}
```

## Limitations

### The Merwe defaults are not tuned for dimension, scale or precision

`merwe_options` defaults to `alpha = 1e-3`, `beta = 2` and `kappa = 0`, the same three constants for every state dimension, every state scale and both `float` and `double`; `try_create` checks only that they are admissible. With `kappa = 0` these defaults give `Wm_0 = 1 - 1/alpha^2`, about `-1e6`, `Wc_0 = 4 - alpha^2 - 1/alpha^2`, also about `-1e6`, and offset weights `1 / (2 alpha^2 NX)`, while the offsets themselves are `1e-3 sqrt(NX)` standard deviations. The mean and covariance are therefore sums of terms near `1e6` in magnitude that cancel, and a state component whose magnitude exceeds its standard deviation by a factor near `1e-3 sqrt(NX) / eps` keeps few significant bits of the offset added to it.

Measured on two two-state problems: the pendulum of the usage example in `float` ends with a position variance of `0.026` against `0.012` in `double`, and a constant-velocity filter in `float` whose position sits at `1e4` with a unit-scale spread reports a position variance of `2e5` against `0.16` and latches `covariance_repaired`. With `alpha = 0.5` both `float` runs come within half a percent of `double`. The negative center weight also removes the guarantee, stated under [update](#update), that a nonlinear measurement model leaves a positive semidefinite posterior. A caller with a `float` filter, a state far from the origin relative to its spread, or a high-dimensional state sets `alpha`, `beta` and `kappa` through `merwe_options` for that problem.

## See Also

- [ekf](ekf.md)<br/> linearization-based alternative
- [manifold-ukf](manifold-ukf.md)<br/> UKF on SO(3) manifold for attitude estimation
- [observer-policy](observer-policy.md)<br/> concept satisfied by this type
- [background/ukf](../../background/ukf.md)<br/> unscented transform derivation
- [guides/estimation/observer-controller](../../guides/estimation/observer-controller.md)<br/> composing observers with controllers
