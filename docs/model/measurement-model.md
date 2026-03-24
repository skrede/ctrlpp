# measurement_model

C++20 concept constraining measurement callables. Used throughout the library to type-check measurement functions passed to estimators (EKF, UKF, particle filter, MHE, NMHE).

## Header

| Form | Header |
|------|--------|
| `measurement_model` (concept) | `#include <ctrlpp/model/measurement_model.h>` |

## Concept Definition

```cpp
template <typename M, typename Scalar, std::size_t NX, std::size_t NY>
concept measurement_model = requires(const M& m,
                                     const Vector<Scalar, NX>& x) {
    { m(x) } -> std::convertible_to<Vector<Scalar, NY>>;
};
```

## What Satisfies It

A type `M` satisfies `measurement_model<M, Scalar, NX, NY>` when it provides a const-callable `operator()` that:

1. Accepts a state vector `Vector<Scalar, NX>`.
2. Returns something convertible to `Vector<Scalar, NY>` (the measurement prediction).

## Example Model Implementation

```cpp
#include "ctrlpp/model/measurement_model.h"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

// Range-bearing measurement from a 2D position state
struct range_bearing
{
    Eigen::Vector2d sensor_pos{0.0, 0.0};

    Eigen::Vector2d operator()(const Eigen::Vector2d& x) const
    {
        Eigen::Vector2d delta = x - sensor_pos;
        double range = delta.norm();
        double bearing = std::atan2(delta[1], delta[0]);
        return {range, bearing};
    }
};

static_assert(ctrlpp::measurement_model<range_bearing, double, 2, 2>);

// Simple linear measurement (position only)
struct position_sensor
{
    Eigen::Matrix<double, 1, 1> operator()(const Eigen::Vector2d& x) const
    {
        return Eigen::Matrix<double, 1, 1>{x[0]};
    }
};

static_assert(ctrlpp::measurement_model<position_sensor, double, 2, 1>);

int main()
{
    range_bearing sensor{.sensor_pos = {0.0, 0.0}};
    Eigen::Vector2d target(3.0, 4.0);

    auto z = sensor(target);
    std::cout << "Range: " << z[0] << "  Bearing: " << z[1] << " rad\n";
}
```

## See Also

- [differentiable-measurement](differentiable-measurement.md) -- extended concept with Jacobian
- [dynamics-model](dynamics-model.md) -- dynamics model concept
- [estimation/ekf](../estimation/ekf.md) -- EKF using measurement_model
