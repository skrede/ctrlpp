# dynamics_model

C++23 concept constraining discrete-time dynamics callables. Used throughout the library to type-check dynamics functions passed to estimators (EKF, UKF, particle filter) and controllers (NMPC, NMHE).

## Header

| Form | Header |
|------|--------|
| `dynamics_model` (concept) | `#include <ctrlpp/model/dynamics_model.h>` |

## Concept Definition

```cpp
template <typename D, typename Scalar, std::size_t NX, std::size_t NU>
concept dynamics_model = requires(const D& d,
                                  const Vector<Scalar, NX>& x,
                                  const Vector<Scalar, NU>& u) {
    { d(x, u) } -> std::convertible_to<Vector<Scalar, NX>>;
};
```

## What Satisfies It

A type `D` satisfies `dynamics_model<D, Scalar, NX, NU>` when it provides a const-callable `operator()` that:

1. Accepts a state vector `Vector<Scalar, NX>` and an input vector `Vector<Scalar, NU>`.
2. Returns something convertible to `Vector<Scalar, NX>` (the next state).

Any callable works: a struct with `operator()`, a lambda, or a function pointer (when wrapped).

## Example Model Implementation

```cpp
#include "ctrlpp/model/dynamics_model.h"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

// Struct satisfying the concept
struct pendulum
{
    double dt{0.01};
    double g{9.81};
    double l{1.0};

    Eigen::Vector2d operator()(const Eigen::Vector2d& x,
                               const Eigen::Matrix<double, 1, 1>& u) const
    {
        double theta = x[0];
        double omega = x[1];
        double alpha = -g / l * std::sin(theta) + u[0];
        return {theta + omega * dt, omega + alpha * dt};
    }
};

static_assert(ctrlpp::dynamics_model<pendulum, double, 2, 1>);

// Lambda also satisfies the concept
int main()
{
    auto linear_dynamics = [](const Eigen::Vector2d& x,
                              const Eigen::Matrix<double, 1, 1>& u)
        -> Eigen::Vector2d
    {
        Eigen::Matrix2d A;
        A << 1.0, 0.01, 0.0, 1.0;
        Eigen::Vector2d B(0.0, 0.01);
        return A * x + B * u[0];
    };

    static_assert(ctrlpp::dynamics_model<decltype(linear_dynamics), double, 2, 1>);

    pendulum pend{.dt = 0.01};
    Eigen::Vector2d x(0.5, 0.0);
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < 100; ++k)
    {
        x = pend(x, u);
        if(k % 25 == 0)
            std::cout << "k=" << k << "  x=[" << x.transpose() << "]\n";
    }
}
```

## See Also

- [differentiable-dynamics](differentiable-dynamics.md) -- extended concept with Jacobians
- [measurement-model](measurement-model.md) -- measurement model concept
- [estimation/ekf](../estimation/ekf.md) -- EKF using dynamics_model
- [mpc/nmpc](../mpc/nmpc.md) -- NMPC using dynamics_model
