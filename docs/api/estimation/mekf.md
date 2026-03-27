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
| `NB` | `std::size_t` | Bias dimension (typically 3 for gyro bias) |
| `NY` | `std::size_t` | Measurement dimension |
| `Measurement` | satisfies `mekf_measurement_model<Scalar, NB, NY>` | Callable: `(Quaternion<Scalar>, Vector<NB>) -> Vector<NY>` |

The error-state dimension is NE = 3 + NB (3 for rotation + NB for bias).

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

## Constructor

```cpp
mekf(Measurement measurement, mekf_config<Scalar, NB, NY> config);
```

CTAD deduction guide available.

## Methods

### predict

```cpp
void predict(const input_vector_t& omega);
void predict(const input_vector_t& omega, Scalar dt);
```

Propagates the nominal quaternion by integrating bias-corrected angular velocity via the SO(3) exponential map, and propagates the error-state covariance through the linearised dynamics.

### update

```cpp
void update(const output_vector_t& z);
```

Computes the measurement Jacobian H (analytically if `differentiable_mekf_measurement` is satisfied, numerically otherwise), applies the Kalman gain to obtain a 6D error-state correction, and composes the rotation correction multiplicatively onto the nominal quaternion. The mandatory frame-change Jacobian G = I - 0.5 * skew(delta_att) is applied to the covariance after the Joseph-form update.

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
        filter.update(g_body);

        auto q_est = filter.attitude();
        Eigen::Vector3d err = ctrlpp::so3::log(q_true.conjugate() * q_est);
        std::cout << k * 0.01 << "," << err.norm() << "\n";
    }
}
```

## See Also

- [manifold-ukf](manifold-ukf.md) -- sigma-point alternative on SO(3)
- [complementary-filter](complementary-filter.md) -- lightweight IMU fusion
- [so3](../lie/so3.md) -- SO(3) quaternion utilities used by MEKF
- [background/attitude-estimation](../../background/attitude-estimation.md) -- MEKF derivation
