# manifold_ukf

Unscented Kalman Filter on the SO(3) manifold with geodesic mean computation. Sigma points are generated in the tangent space (Lie algebra so(3)), mapped onto the manifold via the exponential map, propagated through the dynamics model, and their statistics are recovered using the logarithmic map and iterative geodesic (Frechet) mean computation. This avoids the quaternion norm constraint and singularities of Euler angle representations.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::manifold_ukf<Scalar, NY, Dynamics, Measurement, Strategy>` | `#include <ctrlpp/estimation/manifold_ukf.h>` |

No convenience header exists for this type. Use the categorical path.

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
| `NY` | `std::size_t` | Measurement dimension |
| `Dynamics` | satisfies `manifold_ukf_dynamics_model<Scalar>` | Callable: `(Quaternion<Scalar>, Vector<3>) -> Quaternion<Scalar>` |
| `Measurement` | satisfies `manifold_ukf_measurement_model<Scalar, NY>` | Callable: `(Quaternion<Scalar>) -> Vector<NY>` |
| `Strategy` | satisfies `manifold_sigma_point_strategy<Scalar>` | Manifold sigma point generator (default: `so3_merwe_sigma_points<Scalar>`) |

The state dimension is always 3 (tangent space of SO(3)). The nominal state is a unit quaternion.

## Type Aliases

```cpp
using state_vector_t  = Vector<Scalar, 4>;     // quaternion as (w,x,y,z)
using input_vector_t  = Vector<Scalar, 3>;      // angular velocity
using output_vector_t = Vector<Scalar, NY>;
using cov_matrix_t    = Matrix<Scalar, 3, 3>;   // tangent-space covariance
```

## Config (`manifold_ukf_config`)

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, 3, 3>` | Identity | Process noise covariance in tangent space |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance |
| `q0` | `Eigen::Quaternion<Scalar>` | Identity | Initial quaternion estimate |
| `P0` | `Matrix<Scalar, 3, 3>` | Identity | Initial tangent-space covariance |
| `geodesic_mean_max_iter` | `std::size_t` | `30` | Maximum iterations for geodesic mean computation |
| `geodesic_mean_tol` | `Scalar` | `1e-9` | Convergence tolerance for geodesic mean |

## Construction

### create

```cpp
static auto create(Dynamics dynamics, Measurement measurement,
                                     manifold_ukf_config<Scalar, NY> config,
                                     Strategy strategy = Strategy{})
    -> ctrlpp::expected<manifold_ukf, filter_error>;
```

`create` is the only construction path. There is no non-fallible constructor, so a degenerate configuration is a value the caller has to inspect and never a filter that quietly stands in for one; there is no CTAD deduction guide either, since nothing is left to deduce from. Validates the initial quaternion before the normalization that seeds the filter state: a `q0` with zero or non-finite norm is rejected with `filter_error::degenerate_quaternion` (from `<ctrlpp/estimation/estimation_types.h>`), since normalizing such a quaternion produces NaN and silently poisons the whole filter state. Any finite nonzero `q0` is accepted and normalized. As a static member of a class template, `create` requires explicit template arguments, e.g. `manifold_ukf<double, 3, Dynamics, Measurement>::create(d, m, cfg)`.

## Methods

### predict

```cpp
void predict(const input_vector_t& omega);
```

Generates manifold sigma points around the current quaternion, propagates each through the dynamics model, and computes the geodesic (Frechet) mean and tangent-space covariance of the propagated set. The time step is embedded in the dynamics model, so no `dt` is passed here. If the geodesic-mean iteration exhausts its budget without converging, the filter latches a degraded status reported through `health()`.

### update

```cpp
void update(const output_vector_t& z);
```

Generates manifold sigma points, transforms through measurement model, computes innovation covariance S and manifold cross-covariance Pxz (using log-map deviations), then applies the Kalman gain correction via the exponential map. After the correction, the tangent-space covariance is transported into the corrected frame with the reset Jacobian G = I - 0.5 * skew(delta_phi), mirroring the MEKF reset, so the reported covariance stays anchored to the updated attitude rather than the pre-correction frame.

### state

```cpp
const state_vector_t& state() const;
```

Returns the quaternion as a 4-vector (w, x, y, z).

### covariance

```cpp
const cov_matrix_t& covariance() const;
```

Returns the 3x3 tangent-space covariance.

### innovation

```cpp
const output_vector_t& innovation() const;
```

### attitude

```cpp
const Eigen::Quaternion<Scalar>& attitude() const;
```

Returns the quaternion estimate directly.

### health

```cpp
manifold_ukf_health health() const;
```

Returns the filter-health status, one of `manifold_ukf_health::ok` or `manifold_ukf_health::mean_not_converged`. The status starts at `ok` and latches to `mean_not_converged` the first time a geodesic (Karcher) mean iteration exhausts its iteration budget without meeting `geodesic_mean_tol`. A latched status signals that the predicted attitude may be a non-converged iterate rather than the true geodesic mean.

## Supporting Types

### so3_merwe_sigma_points

Default manifold sigma point strategy. Generates 2*3+1 = 7 sigma points in the tangent space, maps them onto SO(3) via the exponential map composed with the nominal quaternion.

Header: `#include <ctrlpp/estimation/sigma_points/so3_sigma_points.h>`

It shares `merwe_options` with the tangent-space strategy it lifts, and forwards that strategy's parameter-domain check unchanged through its own fallible factory:

```cpp
static auto try_create(options_t opts = options_t{})
    -> ctrlpp::expected<so3_merwe_sigma_points, filter_error>;
```

See [ukf](ukf.md) for the two rejections and the reason each parameter has an exact domain. This strategy keeps the `try_` prefix because it has a real non-fallible counterpart: default construction takes in-domain defaults and cannot fail, so the prefix distinguishes two genuine alternatives rather than restating the return type.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'attitude error (rad)'"

#include <ctrlpp/estimation/manifold_ukf.h>
#include <ctrlpp/lie/so3.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <cmath>
#include <iostream>
#include <random>

int main()
{
    using Scalar = double;
    constexpr std::size_t NY = 3;
    constexpr Scalar dt = 0.01;

    // Gyroscope integration dynamics
    auto dynamics = [dt](const Eigen::Quaterniond& q,
                         const Eigen::Vector3d& omega) -> Eigen::Quaterniond {
        return (q * ctrlpp::so3::exp(omega * dt)).normalized();
    };

    // Accelerometer measurement (gravity in body frame)
    auto measurement = [](const Eigen::Quaterniond& q) -> ctrlpp::Vector<Scalar, NY> {
        return q.toRotationMatrix().transpose() * Eigen::Vector3d{0, 0, 9.81};
    };

    ctrlpp::manifold_ukf_config<Scalar, NY> cfg{
        .Q = Eigen::Matrix3d::Identity() * 0.001,
        .R = Eigen::Matrix3d::Identity() * 0.1,
        .q0 = Eigen::Quaterniond::Identity(),
        .P0 = Eigen::Matrix3d::Identity() * 0.01
    };

    ctrlpp::manifold_ukf filter(dynamics, measurement, cfg);

    std::mt19937 rng(42);
    std::normal_distribution<> noise(0.0, 0.1);

    Eigen::Quaterniond q_true = Eigen::Quaterniond::Identity();
    Eigen::Vector3d omega{0.1, 0.05, 0.02};

    for (int k = 0; k < 500; ++k) {
        q_true = (q_true * ctrlpp::so3::exp(omega * dt)).normalized();

        Eigen::Vector3d g_body = q_true.toRotationMatrix().transpose() * Eigen::Vector3d{0, 0, 9.81};
        for (int i = 0; i < 3; ++i) g_body(i) += noise(rng);

        filter.predict(omega);
        filter.update(g_body);

        auto q_est = filter.attitude();
        Eigen::Vector3d err = ctrlpp::so3::log(q_true.conjugate() * q_est);
        std::cout << k * dt << "," << err.norm() << "\n";
    }
}
```

## See Also

- [mekf](mekf.md)<br/> multiplicative EKF alternative for attitude estimation
- [ukf](ukf.md)<br/> Euclidean UKF for non-manifold problems
- [complementary-filter](complementary-filter.md)<br/> lightweight sensor fusion
- [so3](../lie/so3.md)<br/> SO(3) quaternion utilities used by manifold UKF
- [background/attitude-estimation](../../background/attitude-estimation.md)<br/> manifold filtering theory
