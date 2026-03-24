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

## Constructor

```cpp
explicit complementary_filter(cf_config<Scalar> config);
```

CTAD deduction guide available.

## Methods

### update (IMU -- 6-DOF)

```cpp
void update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel, Scalar dt);
```

Natural IMU update fusing gyroscope angular velocity with accelerometer gravity reference. Computes the rotation error between the expected and measured gravity direction in body frame via cross product, applies PI correction to the gyroscope, and integrates the corrected angular velocity via the SO(3) exponential map.

### update (MARG -- 9-DOF)

```cpp
void update(const Vector<Scalar, 3>& gyro, const Vector<Scalar, 3>& accel,
            const Vector<Scalar, 3>& mag, Scalar dt);
```

MARG update adding magnetometer heading correction to the IMU update. The magnetic field reference direction is computed in the world frame, and the magnetometer error is combined with the accelerometer error for full 3D orientation correction. Falls back to IMU-only update if the magnetometer reading is degenerate.

### ObserverPolicy Interface

```cpp
void predict(const input_vector_t& u);   // stores gyro for next update
void update(const output_vector_t& z);   // calls update(gyro_buf, z, dt)
```

These wrappers satisfy the `ObserverPolicy` concept for composition with controllers.

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

## Usage Example

```cpp
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

        filter.update(gyro, accel, dt);

        auto q_est = filter.attitude();
        Eigen::Vector3d err = ctrlpp::so3::log(q_true.conjugate() * q_est);
        std::cout << k * dt << "," << err.norm() << "\n";
    }
}
```

## See Also

- [mekf](mekf.md) -- optimal attitude estimation with bias tracking
- [manifold-ukf](manifold-ukf.md) -- sigma-point filter on SO(3)
- [reference/attitude-estimation-theory](../reference/attitude-estimation-theory.md) -- complementary filter derivation
