# mekf

Multiplicative Extended Kalman Filter for attitude estimation on SO(3). Maintains a two-track state representation: a 7D nominal state (quaternion + gyroscope bias) with a 6D error-state covariance in the tangent space. The multiplicative formulation avoids the quaternion norm constraint by working with 3D rotation error vectors that are composed back onto the nominal quaternion via the SO(3) exponential map. The mandatory post-update covariance reset via the frame-change Jacobian G is the key correctness concern.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::mekf<Scalar, NB, NY, Measurement>` | `#include <ctrlpp/estimation/mekf.h>` |

No convenience header exists for this type. Use the categorical path.

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
| `NB` | `std::size_t`, `NB >= 3` | Bias dimension (3 for gyro bias alone) |
| `NY` | `std::size_t` | Measurement dimension |
| `Measurement` | satisfies `mekf_measurement_model<Scalar, NB, NY>` | Callable: `(Quaternion<Scalar>, Vector<NB>) -> Vector<NY>` |

The error-state dimension is NE = 3 + NB (3 for rotation + NB for bias).

### Bias dimension lower bound

`NB >= 3` is a hard requirement, enforced by a `static_assert` on both `mekf_config` and `mekf`. The propagation corrects the measured angular rate with the **leading three elements** of the bias vector, a fixed-width slice:

```cpp
Vector<Scalar, 3> omega_corr = omega - b_.template head<3>();
```

Instantiating with `NB < 3` therefore reads past the end of an `NB`-element vector on every `predict()`. Both assertions carry the reason, so the diagnostic names the three-element gyro-bias slice rather than only the bound. A bias vector longer than three elements is accepted: elements beyond the leading three are carried in the state and the covariance but are not consumed by the rate correction.

## Type Aliases

```cpp
using state_vector_t = Vector<Scalar, 4 + NB>;  // quaternion (w,x,y,z) + bias
using input_vector_t = Vector<Scalar, 3>;         // angular velocity (gyroscope)
using output_vector_t = Vector<Scalar, NY>;
using cov_matrix_t   = Matrix<Scalar, NE, NE>;    // 6x6 for NB=3
```

## Config (`mekf_config`)

Where NE = 3 + NB (3 rotation dimensions + NB bias dimensions).

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NE, NE>` | Identity | Error-state process noise covariance |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance |
| `q0` | `Eigen::Quaternion<Scalar>` | Identity | Initial quaternion estimate |
| `b0` | `Vector<Scalar, NB>` | Zero | Initial gyroscope bias estimate |
| `P0` | `Matrix<Scalar, NE, NE>` | Identity | Initial error-state covariance |
| `dt` | `Scalar` | `0.01` | Default time step for the single-argument `predict()` overload |
| `numerical_eps` | `Scalar` | `sqrt(eps)` | Perturbation for numerical measurement Jacobians |

## Construction

### create

```cpp
static auto create(Measurement measurement, mekf_config<Scalar, NB, NY> config)
    -> ctrlpp::expected<mekf, filter_error>;
```

`create` is the only construction path. There is no non-fallible constructor, so a degenerate configuration is a value the caller has to inspect and never a filter that quietly stands in for one; there is no CTAD deduction guide either, since nothing is left to deduce from. Validates the initial quaternion before the normalization that seeds the filter state: a `q0` with zero or non-finite norm is rejected with `filter_error::degenerate_quaternion` (from `<ctrlpp/estimation/estimation_types.h>`), since normalizing such a quaternion produces NaN and silently poisons the whole filter state. Any finite nonzero `q0` is accepted and normalized. As a static member of a class template, `create` requires explicit template arguments, e.g. `mekf<double, 3, 3, Measurement>::create(m, cfg)`.

## Methods

### predict

```cpp
void predict(const input_vector_t& omega);
void predict(const input_vector_t& omega, Scalar dt);
```

Propagates the nominal quaternion by integrating bias-corrected angular velocity via the SO(3) exponential map, and propagates the error-state covariance through the linearized dynamics.

### update

```cpp
auto update(const output_vector_t& z)
    -> ctrlpp::expected<void, mekf_update_error>;
```

The step is rejected **before any member is assigned**, so a rejected step leaves the attitude, the bias, the covariance and the innovation bitwise unchanged, preserves the attitude quaternion's unit norm exactly, and the caller may retry with the next sample. Each cause is an exact domain condition, not a tuning preference.

| Condition | Error | Why it is a separate cause |
|---|---|---|
| the carried attitude quaternion or bias is already non-finite | `mekf_update_error::non_finite_state` | The measurement and its Jacobian are evaluated **at** the nominal state, so a poisoned nominal state makes both meaningless before the measurement is used |
| the carried error-state covariance is already non-finite | `mekf_update_error::non_finite_covariance` | The covariance recursion is driven by the propagation Jacobian and by `Q` and `R`, never by the measurement |
| the supplied measurement has a non-finite component | `mekf_update_error::non_finite_measurement` | The gain carries it into the multiplicative correction, so one such sample makes the attitude quaternion non-finite and no later normalization recovers it |

`predict` is deliberately not fallible. Its input is a command the caller already owns and the plant already took, so refusing it would leave the filter with no propagation for a step that happened. A prediction that poisons the carried estimate is reported by [`health`](#health) instead, which the next `update` latches.

On a step that runs: computes the measurement Jacobian H (analytically if `differentiable_mekf_measurement` is satisfied, numerically otherwise), applies the Kalman gain to obtain a 6D error-state correction, and composes the rotation correction multiplicatively onto the nominal quaternion. The mandatory frame-change Jacobian G = I - 0.5 * skew(delta_att) is applied to the covariance after the Joseph-form update.

### state

```cpp
auto state() const -> const state_vector_t&;
```

Returns the full state vector: [w, x, y, z, bias_0, ..., bias_{NB-1}].

### covariance

```cpp
auto covariance() const -> const cov_matrix_t&;
```

### innovation

```cpp
auto innovation() const -> const output_vector_t&;
```

### attitude

```cpp
auto attitude() const -> Eigen::Quaternion<Scalar>;
```

Returns the current quaternion estimate.

### bias

```cpp
auto bias() const -> const Vector<Scalar, NB>&;
```

Returns the current gyroscope bias estimate.

### health

```cpp
auto health() const -> mekf_health;
```

Returns the persistent state-health status, one of `mekf_health::ok` or `mekf_health::non_finite_estimate`. It answers a question a per-call result cannot, because the question outlives the call: whether the carried estimate is still degraded from a step several samples ago. The status starts at `ok` and latches to `non_finite_estimate` the first time a step finds the carried nominal state or covariance already non-finite, which is how a poisoned gyro rate or a non-finite timestep reaching `predict` becomes visible. A **rejected measurement does not set it**: the rejection mutates nothing, so it leaves the filter healthy. The query carries no discard warning; asking it is optional.

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
        .Q = ctrlpp::Matrix<Scalar, 6, 6>::Identity() * 0.001,
        .R = ctrlpp::Matrix<Scalar, 3, 3>::Identity() * 0.1,
        .q0 = Eigen::Quaterniond::Identity(),
        .b0 = ctrlpp::Vector<Scalar, 3>::Zero(),
        .P0 = ctrlpp::Matrix<Scalar, 6, 6>::Identity() * 0.01,
        .dt = 0.01
    };

    ctrlpp::mekf filter(accel_model, cfg);

    std::mt19937 rng(42);
    std::normal_distribution<> gyro_noise(0.0, 0.01);
    std::normal_distribution<> accel_noise(0.0, 0.1);

    Eigen::Quaterniond q_true = Eigen::Quaterniond::Identity();
    Eigen::Vector3d omega_true{0.1, 0.0, 0.05};

    for (int k = 0; k < 500; ++k) {
        // Simulate true rotation
        q_true = (q_true * ctrlpp::so3::exp(omega_true * 0.01)).normalized();

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
