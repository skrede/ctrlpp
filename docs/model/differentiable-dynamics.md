# differentiable_dynamics

C++23 concept refining `dynamics_model` by additionally requiring analytic Jacobians with respect to state and input. When a dynamics model satisfies this concept, estimators (EKF, MHE) use the analytic Jacobians directly instead of falling back to numerical differentiation.

## Header

| Form | Header |
|------|--------|
| `differentiable_dynamics` (concept) | `#include <ctrlpp/model/differentiable_dynamics.h>` |

## Concept Definition

```cpp
template <typename D, typename Scalar, std::size_t NX, std::size_t NU>
concept differentiable_dynamics =
    dynamics_model<D, Scalar, NX, NU>
    && requires(const D& d,
                const Vector<Scalar, NX>& x,
                const Vector<Scalar, NU>& u) {
        { d.jacobian_x(x, u) } -> std::convertible_to<Matrix<Scalar, NX, NX>>;
        { d.jacobian_u(x, u) } -> std::convertible_to<Matrix<Scalar, NX, NU>>;
    };
```

## What Satisfies It

A type satisfies `differentiable_dynamics` when it satisfies `dynamics_model` and additionally provides:

1. `jacobian_x(x, u)` returning `Matrix<Scalar, NX, NX>` -- the Jacobian of `f(x, u)` with respect to `x`.
2. `jacobian_u(x, u)` returning `Matrix<Scalar, NX, NU>` -- the Jacobian of `f(x, u)` with respect to `u`.

When a dynamics model does **not** satisfy this concept, the EKF and MHE automatically compute Jacobians via central finite differences using the `numerical_eps` configuration parameter.

## Example Model Implementation

```cpp
#include <ctrlpp/model/differentiable_dynamics.h>

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

struct linear_dynamics
{
    Eigen::Matrix2d A;
    Eigen::Vector2d B;

    Eigen::Vector2d operator()(const Eigen::Vector2d& x,
                               const Eigen::Matrix<double, 1, 1>& u) const
    {
        return A * x + B * u[0];
    }

    Eigen::Matrix2d jacobian_x(const Eigen::Vector2d&,
                               const Eigen::Matrix<double, 1, 1>&) const
    {
        return A;
    }

    Eigen::Matrix<double, 2, 1> jacobian_u(const Eigen::Vector2d&,
                                            const Eigen::Matrix<double, 1, 1>&) const
    {
        return B;
    }
};

static_assert(ctrlpp::differentiable_dynamics<linear_dynamics, double, 2, 1>);

int main()
{
    linear_dynamics dyn{
        .A = (Eigen::Matrix2d() << 1.0, 0.01, 0.0, 1.0).finished(),
        .B = Eigen::Vector2d(0.0, 0.01)};

    Eigen::Vector2d x(1.0, 0.0);
    Eigen::Matrix<double, 1, 1> u;
    u << 0.5;

    auto F = dyn.jacobian_x(x, u);
    auto G = dyn.jacobian_u(x, u);

    std::cout << "F (df/dx) =\n" << F << "\n"
              << "G (df/du) = " << G.transpose() << "\n";
}
```

## See Also

- [dynamics-model](dynamics-model.md) -- base dynamics concept
- [estimation/ekf](../estimation/ekf.md) -- EKF preferring analytic Jacobians
