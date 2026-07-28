# complementary_filter

Mahony nonlinear complementary filter for attitude estimation on SO(3). Provides computationally lightweight quaternion-based orientation estimation from IMU (gyro + accelerometer) or MARG (gyro + accelerometer + magnetometer) sensor data. The filter fuses high-frequency gyroscope integration with low-frequency accelerometer/magnetometer corrections using a PI controller on the rotation error, making it well-suited for resource-constrained embedded systems where a full Kalman filter is too expensive.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::complementary_filter<Scalar>` | `#include <ctrlpp/estimation/complementary_filter.h>` |

No convenience header exists for this type. Use the categorical path.

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |

## Type Aliases

```cpp
using state_vector_t  = Vector<Scalar, 7>;   // quaternion (w,x,y,z) + bias (3)
using input_vector_t  = Vector<Scalar, 3>;    // angular velocity (gyroscope)
using output_vector_t = Vector<Scalar, 3>;    // accelerometer reading
```

## Config

```cpp
template <typename Scalar>
struct cf_config
{
    Scalar k_p;                           // proportional gain (default: 2.0)
    Scalar k_i;                           // integral gain for bias estimation (default: 0.005)
    Scalar dt;                            // default time step (default: 0.01)
    Eigen::Quaternion<Scalar> q0;         // initial quaternion (default: identity)
};
```

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `k_p` | `Scalar` | `2.0` | Proportional correction gain. Higher values trust the accelerometer more. |
| `k_i` | `Scalar` | `0.005` | Integral correction gain for online gyroscope bias estimation. |
| `dt` | `Scalar` | `0.01` | Default time step used by the ObserverPolicy interface. |
| `q0` | `Eigen::Quaternion<Scalar>` | identity | Initial orientation estimate. |

## Construction

### create

```cpp
static auto create(cf_config<Scalar> config)
    -> ctrlpp::expected<complementary_filter, filter_error>;
```

`create` is the only construction path. There is no non-fallible constructor, so a degenerate configuration is a value the caller has to inspect and never a filter that quietly stands in for one; there is no CTAD deduction guide either, since nothing is left to deduce from. Validates the initial quaternion before it seeds the filter state: a `q0` with zero or non-finite norm is rejected with `filter_error::degenerate_quaternion` (from `<ctrlpp/estimation/estimation_types.h>`), since normalizing such a quaternion produces NaN and silently poisons the filter state. Any finite nonzero `q0` is accepted and normalized onto the unit sphere at construction; the gravity and magnetic correction terms treat the stored quaternion as a unit rotation. As a static member of a class template, `create` requires explicit template arguments, e.g. `complementary_filter<double>::create(cfg)`.

## Methods

### update (IMU: 6&ndash;DOF)

```cpp
auto update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, Scalar dt)
    -> ctrlpp::expected<void, cf_update_error>;
```

Natural IMU update fusing gyroscope angular velocity with accelerometer gravity reference. Computes the rotation error between the expected and measured gravity direction in body frame via cross product, applies PI correction to the gyroscope, and integrates the corrected angular velocity via the SO(3) exponential map.

### update (MARG: 9&ndash;DOF)

```cpp
auto update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel,
            const Vector<Scalar, 3>& mag, Scalar dt)
    -> ctrlpp::expected<void, cf_update_error>;
```

MARG update adding magnetometer heading correction to the IMU update. The magnetic field reference direction is computed in the world frame, and the magnetometer error is combined with the accelerometer error for full 3D orientation correction. Falls back to IMU-only update if the magnetometer reading is degenerate.

### Rejection, and what is deliberately not a rejection

Both overloads reject **before any member is assigned**, so a rejected step leaves the attitude and the bias bitwise unchanged and the caller may retry with the next sample. The filter carries no covariance, so the carried estimate is the attitude quaternion together with the gyro bias. Each cause is an exact domain condition, not a tuning preference.

| Condition | Error | Why it is a separate cause |
|---|---|---|
| the carried attitude quaternion or bias is already non-finite | `cf_update_error::non_finite_state` | Both correction terms are computed from the carried attitude's rotation matrix, so the fault is upstream of the sensors and is reported ahead of them |
| a supplied sensor vector (rate, acceleration, or magnetic field where present) has a non-finite component | `cf_update_error::non_finite_measurement` | Every one of them reaches the quaternion integration, whose output is the filter's carried memory, so one such sample destroys the attitude permanently |
| the supplied integration step is non-finite | `cf_update_error::non_finite_timestep` | It names a broken clock rather than a broken sensor, so the caller fixes a different input. It poisons the integration just as surely: the tangent vector is the step times the corrected rate |

**A sensor reading with no direction removes its own correction term and nothing else.** The two correction terms are independent sums in Mahony's law, so each is applied when its own vector has a direction and omitted when it does not. **The rate is integrated either way.** The body rotated whether or not the accelerometer could say which way is down, so a step that dropped the integration would discard a rotation that happened and still report success. Reporting the missing correction through the failure channel would be the opposite error: the step did exactly what the algorithm prescribes, and a caller who learns that channel carries non-failures will eventually ignore a real one.

**There is no magnitude threshold on a sensor reading.** A reading is expressed in the caller's units, so an absolute floor would give the same physical acceleration different treatment depending on whether it is reported in g or in millimetres per second squared. Only the reading's *direction* enters either correction, and a direction is scale-free. The coefficients are divided by their largest magnitude before the norm is formed, so the squared norm always lies in `[1, 3]` and neither end of the exponent range is reachable: every finite nonzero reading has a direction, **including one whose squared norm would underflow to zero or overflow to infinity if the norm were formed directly**. The only reading with no direction is the exactly zero vector, which carries none at all rather than a small one. This mirrors `so3::normalize`, which repairs the identical defect on the attitude quaternion.

### ObserverPolicy Interface

```cpp
void predict(const input_vector_t& u);   // stores gyro for next update
auto update(const output_vector_t& z)    // calls update(gyro_buf, z, dt)
    -> ctrlpp::expected<void, cf_update_error>;
```

These wrappers satisfy the `ObserverPolicy` concept for composition with controllers. The single-measurement form **forwards** the IMU overload's result rather than swallowing it, so the caller sees the same cause it would have seen through the three-argument form.

`predict` is deliberately not fallible: it only buffers the rate the caller supplied, and a non-finite rate is caught by the very next `update` with `cf_update_error::non_finite_measurement`.

### state

```cpp
auto state() const -> const state_vector_t&;
```

Returns the 7D state: [w, x, y, z, bias_x, bias_y, bias_z].

### attitude

```cpp
auto attitude() const -> Eigen::Quaternion<Scalar>;
```

Returns the quaternion estimate.

### bias

```cpp
auto bias() const -> const Vector<Scalar, 3>&;
```

Returns the estimated gyroscope bias.

### health

```cpp
auto health() const -> cf_health;
```

Returns the persistent state-health status, one of `cf_health::ok` or `cf_health::non_finite_estimate`. It answers a question a per-call result cannot, because the question outlives the call: whether the carried attitude is still degraded from a step several samples ago. The status starts at `ok` and latches to `non_finite_estimate` the first time a step finds the carried attitude or bias already non-finite. A **rejected update does not set it**: the rejection mutates nothing, so it leaves the filter healthy. The query carries no discard warning; asking it is optional.

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

    ctrlpp::cf_config<Scalar> cfg{
        .k_p = 2.0,
        .k_i = 0.005,
        .dt = dt,
        .q0 = Eigen::Quaterniond::Identity()
    };

    ctrlpp::complementary_filter filter(cfg);

    std::mt19937 rng(42);
    std::normal_distribution<> gyro_noise(0.0, 0.01);
    std::normal_distribution<> accel_noise(0.0, 0.05);

    Eigen::Quaterniond q_true = Eigen::Quaterniond::Identity();
    Eigen::Vector3d omega{0.1, 0.0, 0.05};

    for (int k = 0; k < 1000; ++k) {
        q_true = (q_true * ctrlpp::so3::exp(omega * dt)).normalized();

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
