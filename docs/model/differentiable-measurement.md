# differentiable_measurement

C++23 concept refining `measurement_model` by additionally requiring an analytic Jacobian of the measurement function with respect to state. When satisfied, estimators use the analytic Jacobian directly instead of numerical differentiation.

## Header

| Form | Header |
|------|--------|
| `differentiable_measurement` (concept) | `#include <ctrlpp/model/differentiable_measurement.h>` |

## Concept Definition

```cpp
template <typename M, typename Scalar, std::size_t NX, std::size_t NY>
concept differentiable_measurement =
    measurement_model<M, Scalar, NX, NY>
    && requires(const M& m, const Vector<Scalar, NX>& x) {
        { m.jacobian(x) } -> std::convertible_to<Matrix<Scalar, NY, NX>>;
    };
```

## What Satisfies It

A type satisfies `differentiable_measurement` when it satisfies `measurement_model` and additionally provides:

1. `jacobian(x)` returning `Matrix<Scalar, NY, NX>` -- the Jacobian of `h(x)` with respect to `x`.

When a measurement model does **not** satisfy this concept, the EKF and MHE fall back to central finite-difference Jacobians.

## Example Model Implementation

```cpp
#include "ctrlpp/model/differentiable_measurement.h"

#include <Eigen/Dense>

#include <iostream>

struct linear_measurement
{
    Eigen::RowVector2d C{1.0, 0.0};

    Eigen::Matrix<double, 1, 1> operator()(const Eigen::Vector2d& x) const
    {
        return Eigen::Matrix<double, 1, 1>{C * x};
    }

    Eigen::Matrix<double, 1, 2> jacobian(const Eigen::Vector2d&) const
    {
        return C;
    }
};

static_assert(ctrlpp::differentiable_measurement<linear_measurement, double, 2, 1>);

int main()
{
    linear_measurement sensor;
    Eigen::Vector2d x(3.0, 1.5);

    auto z = sensor(x);
    auto H = sensor.jacobian(x);

    std::cout << "Measurement: " << z[0] << "\n"
              << "Jacobian H: " << H << "\n";
}
```

## See Also

- [measurement-model](measurement-model.md) -- base measurement concept
- [differentiable-dynamics](differentiable-dynamics.md) -- analogous concept for dynamics
- [estimation/ekf](../estimation/ekf.md) -- EKF preferring analytic Jacobians
