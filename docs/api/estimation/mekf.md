# mekf

Multiplicative Extended Kalman Filter for attitude estimation on SO(3). Maintains a two-track state representation: a nominal state of `4 + NB` components (unit quaternion and bias) with an error-state covariance of dimension `3 + NB` (a 3D rotation error in the tangent space and the bias error). The multiplicative formulation avoids the quaternion norm constraint by working with 3D rotation error vectors that are composed back onto the nominal quaternion on the right, `q = q_nominal * exp(delta)`, via the SO(3) exponential map. After each update the covariance is carried into the corrected frame by the reset Jacobian `G`.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::mekf<Scalar, NB, NY, Measurement>` | `#include <ctrlpp/estimation/mekf.h>` |
| `ctrlpp::mekf<Scalar, NB, NY, Measurement>` | `#include <ctrlpp/mekf.h>` (convenience) |
| `ctrlpp::mekf_config`, `ctrlpp::mekf_update_error`, `ctrlpp::mekf_health`, `ctrlpp::mekf_measurement_model`, `ctrlpp::differentiable_mekf_measurement` | `#include <ctrlpp/estimation/mekf_types.h>` (included by `mekf.h`) |

`mekf` satisfies `ObserverPolicy` and `CovarianceObserver` ([observer-policy](observer-policy.md)).

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type (e.g. `double`, `float`) |
| `NB` | `std::size_t`, `NB >= 3` | Bias dimension (3 for gyro bias alone) |
| `NY` | `std::size_t`, `NY > 0` | Measurement dimension |
| `Measurement` | satisfies `mekf_measurement_model<Measurement, Scalar, NB, NY>` | Callable: `(const Quaternion<Scalar>&, const Vector<Scalar, NB>&) -> Vector<Scalar, NY>` |

A measurement model that also provides `jacobian(q, b)` returning `Matrix<Scalar, NY, 3 + NB>` satisfies `differentiable_mekf_measurement`, and the filter uses that Jacobian instead of central differences. The Jacobian is taken with respect to the right-multiplicative rotation error and then the bias, the same error state the covariance describes.

The error-state dimension is NE = 3 + NB (3 for rotation + NB for bias).

### Bias dimension lower bound

`NB >= 3` is a hard requirement, enforced by a `static_assert` on both `mekf_config` and `mekf`. The propagation corrects the measured angular rate with the **leading three elements** of the bias vector, a fixed-width slice:

```cpp
Vector<Scalar, 3> omega_corr = omega - m_b.template head<3>();
```

Instantiating with `NB < 3` would read past the end of an `NB`-element vector on every `predict()`, which is why it does not compile. Both assertions carry the reason, so the diagnostic names the three-element gyro-bias slice rather than only the bound. A bias vector longer than three elements is accepted: elements beyond the leading three are carried in the state and the covariance but are not consumed by the rate correction.

## Type Aliases

```cpp
using state_vector_t = Vector<Scalar, 4 + NB>;  // quaternion (w,x,y,z) + bias
using input_vector_t = Vector<Scalar, 3>;         // angular velocity (gyroscope)
using output_vector_t = Vector<Scalar, NY>;
using cov_matrix_t   = Matrix<Scalar, NE, NE>;    // 6x6 for NB=3
using meas_cov_t     = Matrix<Scalar, NY, NY>;
```

## Config (`mekf_config`)

Where NE = 3 + NB (3 rotation dimensions + NB bias dimensions).

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NE, NE>` | Identity | Discrete error-state process noise covariance, added once per `predict` call (see [predict](#predict)) |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance |
| `q0` | `Eigen::Quaternion<Scalar>` | Identity | Initial quaternion estimate, unit within rounding |
| `b0` | `Vector<Scalar, NB>` | Zero | Initial gyroscope bias estimate |
| `P0` | `Matrix<Scalar, NE, NE>` | Identity | Initial error-state covariance |
| `dt` | `Scalar` | `0.01` | Default time step for the single-argument `predict()` overload |
| `numerical_eps` | `Scalar` | `cbrt(eps)` | Absolute central-difference step for the numerical measurement Jacobian: each rotation column perturbs the quaternion by `exp(+/- numerical_eps * e_i)` on the right, each bias column perturbs the bias component by `numerical_eps` rounded to a representable difference (read only when the model provides no `jacobian`) |

## Construction

### create

```cpp
static ctrlpp::expected<mekf, filter_error> create(Measurement measurement, mekf_config<Scalar, NB, NY> config);
```

`create` is the only construction path. There is no non-fallible constructor, so a degenerate configuration is a value the caller has to inspect and never a filter that quietly stands in for one; there is no CTAD deduction guide either, since nothing is left to deduce from. Validates the initial quaternion first, before the normalization that seeds the filter state: a `q0` with zero or non-finite norm is rejected with `filter_error::degenerate_quaternion` (from `<ctrlpp/estimation/estimation_types.h>`), since normalizing such a quaternion produces NaN and silently poisons the whole filter state. A `q0` whose norm differs from one by more than the rounding bound `max(4 eps, 4 denorm_min)` is rejected with `filter_error::non_unit_quaternion`: a rotation quaternion is unit by definition, so `(2, 0, 0, 0)` names no attitude. A norm defect within the bound is rounding drift, the residue of the arithmetic that formed `q0`, and is normalized away; every quaternion produced by a working-precision normalization, an angle-axis conversion, `so3::exp` or a composition of such rotations is within it. The convention and its count are in [Rounding drift is accepted, structure is refused](../../guides/patterns/numerical-behavior.md#rounding-drift-is-accepted-structure-is-refused). As a static member of a class template, `create` requires explicit template arguments, e.g. `mekf<double, 3, 3, Measurement>::create(m, cfg)`.

After the quaternion, `create` checks the remaining fields in the order the config declares them, each matrix for finiteness, then for symmetry, then for definiteness, as every Kalman-family factory does:

| Enumerator | Condition |
|------------|-----------|
| `filter_error::non_finite_process_noise` / `asymmetric_process_noise` / `indefinite_process_noise` | `Q` has a non-finite entry / an asymmetry beyond the rounding bound / an eigenvalue below `-floor` |
| `filter_error::non_finite_measurement_noise` / `asymmetric_measurement_noise` / `indefinite_measurement_noise` | `R` has a non-finite entry / an asymmetry beyond the rounding bound / an eigenvalue below `-floor` |
| `filter_error::non_finite_initial_bias` | `b0` has a non-finite component |
| `filter_error::non_finite_initial_covariance` / `asymmetric_initial_covariance` / `indefinite_initial_covariance` | `P0` has a non-finite entry / an asymmetry beyond the rounding bound / an eigenvalue below `-floor` |

The bound is `max(N eps max|M_ij|, N denorm_min)` for both checks; as a definiteness floor it is the order of the symmetric eigensolver's backward error, kept from underflowing to zero. An asymmetry `max|M_ij - M_ji|` within it is rounding drift: the matrix is accepted, stored as its symmetrized form (each differing pair replaced by its midpoint), and judged for definiteness as that matrix. A matrix that is not symmetric positive semidefinite is not a covariance; the convention and its count are in [Rounding drift is accepted, structure is refused](../../guides/patterns/numerical-behavior.md#rounding-drift-is-accepted-structure-is-refused). A singular positive semidefinite matrix, such as `Q = 0` for a constant bias or a point-mass `P0`, is accepted. `dt` and `numerical_eps` are not validated here.

## Methods

### predict

```cpp
void predict(const input_vector_t& omega);
void predict(const input_vector_t& omega, Scalar dt);
```

The single-argument overload uses `dt` from the configuration. With the bias-corrected rate `w = omega - b.head<3>()`, the nominal attitude becomes `q * exp(w dt)`, normalized; the bias is held. The error-state covariance becomes `F P F' + Q`, symmetrized, where `F` is the identity except in its first three rows: the attitude block is `C'` with `C` the rotation matrix of `exp(w dt)`, and the attitude-bias block is `-dt` times the leading `3 x NB` identity, so only the first three bias components drive the attitude error.

**`Q` is a discrete covariance added unchanged once per call, for any `dt`.** It is not scaled by the step, so the process noise per unit time grows as `dt` shrinks. A caller whose sample period varies, or who calls the two-argument overload with a `dt` other than the one `Q` was discretized for, must account for that in `Q` itself; the filter does not convert a continuous-time noise density. `dt` is not validated: a non-finite `dt` poisons the attitude, which the next `update` reports and [`health`](#health) latches.

### update

```cpp
ctrlpp::expected<void, mekf_update_error> update(const output_vector_t& z);
```

A refused step leaves the attitude, the bias, the covariance and the innovation bitwise unchanged, preserves the attitude quaternion's unit norm exactly, and the caller may retry with the next sample. The first three causes are checked before anything is computed; the fourth is checked on the computed result, which is then discarded. Each cause is an exact domain condition, not a tuning preference.

| Condition | Error | Why it is a separate cause |
|---|---|---|
| the carried attitude quaternion or bias is already non-finite | `mekf_update_error::non_finite_state` | The measurement and its Jacobian are evaluated **at** the nominal state, so a poisoned nominal state makes both meaningless before the measurement is used |
| the carried error-state covariance is already non-finite | `mekf_update_error::non_finite_covariance` | The covariance recursion is driven by the propagation Jacobian and by `Q` and `R`, never by the measurement |
| the supplied measurement has a non-finite component | `mekf_update_error::non_finite_measurement` | The gain carries it into the multiplicative correction, so one such sample makes the attitude quaternion non-finite and no later normalization recovers it |
| the updated attitude, bias, covariance or innovation is non-finite although the carried estimate and the measurement were finite | `mekf_update_error::non_finite_result` | The step's own arithmetic failed: the measurement model or its Jacobian returned a non-finite value at the nominal state, a product overflowed, or the innovation covariance was zero with one output |

`predict` is deliberately not fallible. Its input is a command the caller already owns and the plant already took, so refusing it would leave the filter with no propagation for a step that happened. A prediction that poisons the carried estimate is reported by [`health`](#health) instead, which the next `update` latches.

On a step that runs: the innovation is `z - h(q, b)` at the nominal state; the measurement Jacobian `H` is `measurement.jacobian(q, b)` when `differentiable_mekf_measurement` is satisfied and central differences otherwise. Then `S = HPH' + R`, and the gain `K = PH'S^-1` is solved with a column-pivoting QR factorization of `S` (a division when `NY = 1`). The error-state correction `dx = K(z - h(q, b))` is applied as `q = q * exp(dx_att)`, normalized, and `b = b + dx_bias`. The covariance takes the Joseph form `(I-KH)P(I-KH)' + KRK'` and is then carried into the corrected frame by `G P G'`, symmetrized, where `G` is the identity with its attitude block `I - 0.5 skew(dx_att)`.

### state

```cpp
const state_vector_t& state() const;
```

Returns the full state vector: [w, x, y, z, bias_0, ..., bias_{NB-1}].

### covariance

```cpp
const cov_matrix_t& covariance() const;
```

### innovation

```cpp
const output_vector_t& innovation() const;
```

### attitude

```cpp
Eigen::Quaternion<Scalar> attitude() const;
```

Returns the current quaternion estimate.

### bias

```cpp
const Vector<Scalar, NB>& bias() const;
```

Returns the current gyroscope bias estimate.

### health

```cpp
mekf_health health() const;
```

Returns the persistent state-health status, one of `mekf_health::ok` or `mekf_health::non_finite_estimate`. It answers a question a per-call result cannot, because the question outlives the call: whether the carried estimate is still degraded from a step several samples ago. The status starts at `ok` and latches to `non_finite_estimate` when an `update` is refused with `non_finite_state` or `non_finite_covariance`, which is how a poisoned gyro rate or a non-finite timestep reaching `predict` becomes visible. Nothing clears it. An update refused with `non_finite_measurement` or `non_finite_result` does not set it, since neither leaves a non-finite estimate behind.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'attitude error (rad)'"

#include <ctrlpp/estimation/mekf.h>
#include <ctrlpp/lie/so3.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <random>

int main()
{
    using Scalar = double;
    constexpr std::size_t NB = 3, NY = 3;

    // Accelerometer measurement model: rotated gravity
    auto accel_model = [](const Eigen::Quaterniond& q,
                          const ctrlpp::Vector<Scalar, NB>& /*bias*/) -> ctrlpp::Vector<Scalar, NY> {
        Eigen::Vector3d g_world{0.0, 0.0, 9.81};
        return q.toRotationMatrix().transpose() * g_world;
    };

    ctrlpp::mekf_config<Scalar, NB, NY> cfg{
        .Q = (ctrlpp::Vector<Scalar, 6>() << 1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8).finished().asDiagonal(),
        .R = ctrlpp::Matrix<Scalar, 3, 3>::Identity() * 0.1,
        .P0 = ctrlpp::Matrix<Scalar, 6, 6>::Identity() * 0.01
    };

    using filter_t = ctrlpp::mekf<Scalar, NB, NY, decltype(accel_model)>;
    auto created = filter_t::create(accel_model, cfg);
    if(!created)
    {
        std::cerr << "invalid MEKF configuration\n";
        return 1;
    }
    auto& filter = *created;

    std::mt19937 rng(42);
    std::normal_distribution<> gyro_noise(0.0, 0.01);
    std::normal_distribution<> accel_noise(0.0, 0.1);

    Eigen::Quaterniond q_true = Eigen::Quaterniond::Identity();
    Eigen::Vector3d omega_true{0.1, 0.0, 0.05};
    Eigen::Vector3d step = omega_true * 0.01;

    for (int k = 0; k < 500; ++k) {
        // Simulate true rotation
        q_true = (q_true * ctrlpp::so3::exp(step)).normalized();

        // Noisy gyro
        Eigen::Vector3d gyro = omega_true;
        for (int i = 0; i < 3; ++i) gyro(i) += gyro_noise(rng);

        // Noisy accelerometer
        Eigen::Vector3d g_body = q_true.toRotationMatrix().transpose() * Eigen::Vector3d{0, 0, 9.81};
        for (int i = 0; i < 3; ++i) g_body(i) += accel_noise(rng);

        filter.predict(gyro);
        if(!filter.update(g_body))
        {
            std::cerr << "MEKF rejected the measurement\n";
            return 1;
        }

        auto q_est = filter.attitude();
        Eigen::Vector3d err = ctrlpp::so3::log(q_true.conjugate() * q_est);
        std::cout << k * 0.01 << "," << err.norm() << "\n";
    }
}
```

## See Also

- [manifold-ukf](manifold-ukf.md)<br/> sigma-point alternative on SO(3)
- [complementary-filter](complementary-filter.md)<br/> lightweight IMU fusion
- [so3](../lie/so3.md)<br/> SO(3) quaternion utilities used by MEKF
- [background/attitude-estimation](../../background/attitude-estimation.md)<br/> MEKF derivation
