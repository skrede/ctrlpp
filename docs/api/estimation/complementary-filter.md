# complementary_filter

Mahony nonlinear complementary filter for attitude estimation on SO(3). Provides computationally lightweight quaternion-based orientation estimation from IMU (gyro + accelerometer) or MARG (gyro + accelerometer + magnetometer) sensor data. The filter fuses high-frequency gyroscope integration with low-frequency accelerometer/magnetometer corrections using a PI controller on the rotation error, making it well-suited for resource-constrained embedded systems where a full Kalman filter is too expensive. It carries no covariance.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::complementary_filter<Scalar>`, `ctrlpp::cf_config`, `ctrlpp::cf_update_error`, `ctrlpp::cf_health` | `#include <ctrlpp/estimation/complementary_filter.h>` |
| `ctrlpp::complementary_filter<Scalar>` | `#include <ctrlpp/complementary_filter.h>` (convenience) |

`complementary_filter` satisfies `ObserverPolicy` ([observer-policy](observer-policy.md)), not `CovarianceObserver`.

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` | Floating-point type (e.g. `double`, `float`) |

## Type Aliases

```cpp
using state_vector_t  = Vector<Scalar, 7>;   // quaternion (w,x,y,z) + bias (3)
using input_vector_t  = Vector<Scalar, 3>;    // angular velocity (gyroscope)
using output_vector_t = Vector<Scalar, 3>;    // accelerometer reading
```

## Config

```cpp
template <ctrlpp_floating_scalar Scalar>
struct cf_config
{
    Scalar k_p{Scalar{2}};
    Scalar k_i{static_cast<Scalar>(0.005)};
    Scalar dt{static_cast<Scalar>(0.01)};
    Eigen::Quaternion<Scalar> q0{Eigen::Quaternion<Scalar>::Identity()};
};
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `k_p` | `Scalar` | `2.0` | Proportional correction gain. Higher values trust the accelerometer and magnetometer more. |
| `k_i` | `Scalar` | `0.005` | Integral correction gain for online gyroscope bias estimation. |
| `dt` | `Scalar` | `0.01` | Time step used by the single-argument `update` of the ObserverPolicy interface; the sensor overloads take their own `dt`. |
| `q0` | `Eigen::Quaternion<Scalar>` | identity | Initial orientation estimate, unit within rounding. |

## Construction

### create

```cpp
static ctrlpp::expected<complementary_filter, filter_error> create(cf_config<Scalar> config);
```

`create` is the only construction path. There is no non-fallible constructor, so a degenerate configuration is a value the caller has to inspect and never a filter that quietly stands in for one; there is no CTAD deduction guide either, since nothing is left to deduce from. Validates the initial quaternion before it seeds the filter state: a `q0` with zero or non-finite norm is rejected with `filter_error::degenerate_quaternion` (from `<ctrlpp/estimation/estimation_types.h>`), since normalizing such a quaternion produces NaN and silently poisons the filter state. A `q0` whose norm differs from one by more than the rounding bound `max(4 eps, 4 denorm_min)` is rejected with `filter_error::non_unit_quaternion`: a rotation quaternion is unit by definition, so `(2, 0, 0, 0)` names no attitude. A norm defect within the bound is rounding drift and is normalized away at construction, so the gravity and magnetic correction terms, which treat the stored quaternion as a unit rotation, start from one; every quaternion produced by a working-precision normalization, an angle-axis conversion, `so3::exp` or a composition of such rotations is within it. The convention and its count are in [Rounding drift is accepted, structure is refused](../../guides/patterns/numerical-behavior.md#rounding-drift-is-accepted-structure-is-refused). As a static member of a class template, `create` requires explicit template arguments, e.g. `complementary_filter<double>::create(cfg)`. The quaternion is the only field `create` checks: `k_p`, `k_i` and `dt` are taken as given, and a non-finite configured `dt` is refused with `non_finite_timestep` by each single-argument `update`. The bias starts at zero.

## Methods

### update (IMU: 6&ndash;DOF)

```cpp
ctrlpp::expected<void, cf_update_error> update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, Scalar dt);
```

IMU update fusing gyroscope angular velocity with the accelerometer's gravity reference. With `a` the accelerometer direction and `v = R(q)' e_z` the world `+z` axis seen in the body frame, the error is `e = a x v`; a level filter at rest therefore expects the reading along body `+z`, the specific force an accelerometer reports. The step then updates the bias and integrates the corrected rate over the full step:

```cpp
bias -= k_i * e * dt;
omega_c = gyro - bias + k_p * e;
q = (q * so3::exp(dt * omega_c)).normalized();
```

### update (MARG: 9&ndash;DOF)

```cpp
ctrlpp::expected<void, cf_update_error> update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel,
                                               const Vector<Scalar, 3>& mag, Scalar dt);
```

MARG update adding a magnetometer correction to the IMU update. The magnetometer direction `m` is rotated into the world frame, `h = R(q) m`; its reference is `b = [sqrt(h_x^2 + h_y^2), 0, h_z]`, the same field with its horizontal part on the world x axis; the error term `m x R(q)' b` is added to the gravity term. The magnetic term is a full cross product, not a projection onto heading, so it also acts on roll and pitch. A reading with no direction drops its own term, as below.

### Rejection, and what is deliberately not a rejection

A refused step leaves the attitude and the bias bitwise unchanged and the caller may retry with the next sample. The first three causes are checked before anything is computed; the fourth is checked on the computed result, which is then discarded. The filter carries no covariance, so the carried estimate is the attitude quaternion together with the gyro bias. Each cause is an exact domain condition, not a tuning preference.

| Condition | Error | Why it is a separate cause |
|---|---|---|
| the carried attitude quaternion or bias is already non-finite | `cf_update_error::non_finite_state` | Both correction terms are computed from the carried attitude's rotation matrix, so the fault is upstream of the sensors and is reported ahead of them |
| a supplied sensor vector (rate, acceleration, or magnetic field where present) has a non-finite component | `cf_update_error::non_finite_measurement` | Every one of them reaches the quaternion integration, whose output is the filter's carried memory, so one such sample destroys the attitude permanently |
| the supplied integration step is non-finite | `cf_update_error::non_finite_timestep` | It names a broken clock rather than a broken sensor, so the caller fixes a different input. It poisons the integration just as surely: the tangent vector is the step times the corrected rate |
| the updated attitude or bias is non-finite although every input was finite | `cf_update_error::non_finite_result` | The step's own arithmetic overflowed, for example a gain or step large enough that the bias or the tangent vector leaves the representable range |

**A sensor reading with no direction removes its own correction term and nothing else.** The two correction terms are independent sums in Mahony's law, so each is applied when its own vector has a direction and omitted when it does not. **The rate is integrated either way.** The body rotated whether or not the accelerometer could say which way is down, so a step that dropped the integration would discard a rotation that happened and still report success. Reporting the missing correction through the failure channel would be the opposite error: the step did exactly what the algorithm prescribes, and a caller who learns that channel carries non-failures will eventually ignore a real one.

**There is no magnitude threshold on a sensor reading.** A reading is expressed in the caller's units, so an absolute floor would give the same physical acceleration different treatment depending on whether it is reported in g or in millimeters per second squared. Only the reading's *direction* enters either correction, and a direction is scale-free. The coefficients are divided by their largest magnitude before the norm is formed, so the squared norm always lies in `[1, 3]` and neither end of the exponent range is reachable: every finite nonzero reading has a direction, **including one whose squared norm would underflow to zero or overflow to infinity if the norm were formed directly**. The only reading with no direction is the exactly zero vector, which carries none at all rather than a small one. This mirrors `so3::normalize`, which repairs the identical defect on the attitude quaternion.

### ObserverPolicy Interface

```cpp
void predict(const input_vector_t& u);   // stores gyro for next update
ctrlpp::expected<void, cf_update_error> update(const output_vector_t& z);   // calls update(gyro_buf, z, dt)
```

These wrappers satisfy the `ObserverPolicy` concept for composition with controllers. The single-measurement form **forwards** the IMU overload's result rather than swallowing it, so the caller sees the same cause it would have seen through the three-argument form.

`predict` is deliberately not fallible: it only buffers the rate the caller supplied, and a non-finite rate is caught by the very next `update` with `cf_update_error::non_finite_measurement`.

### state

```cpp
const state_vector_t& state() const;
```

Returns the 7D state: [w, x, y, z, bias_x, bias_y, bias_z].

### attitude

```cpp
Eigen::Quaternion<Scalar> attitude() const;
```

Returns the quaternion estimate.

### bias

```cpp
const Vector<Scalar, 3>& bias() const;
```

Returns the estimated gyroscope bias.

### health

```cpp
cf_health health() const;
```

Returns the persistent state-health status, one of `cf_health::ok` or `cf_health::non_finite_estimate`. It answers a question a per-call result cannot, because the question outlives the call: whether the carried attitude is still degraded from a step several samples ago. The status starts at `ok` and latches to `non_finite_estimate` when an update is refused with `non_finite_state`. Nothing clears it. The other three refusals do not set it, since none leaves a non-finite estimate behind.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'attitude error (rad)'"

#include <ctrlpp/estimation/complementary_filter.h>
#include <ctrlpp/lie/so3.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <random>

int main()
{
    using Scalar = double;
    constexpr Scalar dt = 0.01;

    auto created = ctrlpp::complementary_filter<Scalar>::create(ctrlpp::cf_config<Scalar>{});
    if(!created)
    {
        std::cerr << "invalid complementary filter configuration\n";
        return 1;
    }
    auto& filter = *created;

    std::mt19937 rng(42);
    std::normal_distribution<> gyro_noise(0.0, 0.01);
    std::normal_distribution<> accel_noise(0.0, 0.05);

    Eigen::Quaterniond q_true = Eigen::Quaterniond::Identity();
    Eigen::Vector3d omega{0.1, 0.0, 0.05};
    Eigen::Vector3d step = omega * dt;

    for (int k = 0; k < 1000; ++k) {
        q_true = (q_true * ctrlpp::so3::exp(step)).normalized();

        // Noisy gyro
        Eigen::Vector3d gyro = omega;
        for (int i = 0; i < 3; ++i) gyro(i) += gyro_noise(rng);

        // Noisy accelerometer (gravity in body frame)
        Eigen::Vector3d accel = q_true.toRotationMatrix().transpose() * Eigen::Vector3d{0, 0, 9.81};
        for (int i = 0; i < 3; ++i) accel(i) += accel_noise(rng);

        if(!filter.update(gyro, accel, dt))
        {
            std::cerr << "complementary filter rejected the sample\n";
            return 1;
        }

        auto q_est = filter.attitude();
        Eigen::Vector3d err = ctrlpp::so3::log(q_true.conjugate() * q_est);
        std::cout << k * dt << "," << err.norm() << "\n";
    }
}
```

## See Also

- [mekf](mekf.md)<br/> optimal attitude estimation with bias tracking
- [manifold-ukf](manifold-ukf.md)<br/> sigma-point filter on SO(3)
- [so3](../lie/so3.md)<br/> SO(3) quaternion utilities
- [background/attitude-estimation](../../background/attitude-estimation.md)<br/> complementary filter derivation
