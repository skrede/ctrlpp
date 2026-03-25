# constraint_model

C++23 concepts for nonlinear inequality constraints used in NMPC and NMHE. Defines both path constraints (applied at each time step) and terminal constraints (applied at the final time step).

## Header

| Form | Header |
|------|--------|
| `constraint_model` (concept) | `#include <ctrlpp/model/constraint_model.h>` |
| `terminal_constraint_model` (concept) | `#include <ctrlpp/model/constraint_model.h>` |

## Concept Definitions

### constraint_model (path constraints)

```cpp
template <typename G, typename Scalar,
          std::size_t NX, std::size_t NU, std::size_t NC>
concept constraint_model = requires(const G& g,
                                    const Vector<Scalar, NX>& x,
                                    const Vector<Scalar, NU>& u) {
    { g(x, u) } -> std::convertible_to<Vector<Scalar, NC>>;
};
```

Represents NC nonlinear path inequality constraints `g(x, u) <= 0` applied at each stage of the prediction horizon.

### terminal_constraint_model

```cpp
template <typename H, typename Scalar,
          std::size_t NX, std::size_t NTC>
concept terminal_constraint_model = requires(const H& h,
                                             const Vector<Scalar, NX>& x) {
    { h(x) } -> std::convertible_to<Vector<Scalar, NTC>>;
};
```

Represents NTC nonlinear terminal inequality constraints `h(x_N) <= 0` applied at the final prediction step.

## What Satisfies Them

**constraint_model:** Any const-callable taking `(Vector<NX>, Vector<NU>)` and returning `Vector<NC>`. The constraint is satisfied when all elements of the returned vector are non-positive.

**terminal_constraint_model:** Any const-callable taking `(Vector<NX>)` and returning `Vector<NTC>`. Same non-positivity convention.

## Example Model Implementation

```cpp
#include <ctrlpp/model/constraint_model.h>

#include <Eigen/Dense>

#include <iostream>

// Path constraint: keep position within a circle of radius R
struct circle_constraint
{
    double R{5.0};

    Eigen::Matrix<double, 1, 1> operator()(
        const Eigen::Vector2d& x,
        const Eigen::Matrix<double, 1, 1>&) const
    {
        // g(x,u) = x'x - R^2 <= 0
        return Eigen::Matrix<double, 1, 1>{x.squaredNorm() - R * R};
    }
};

static_assert(ctrlpp::constraint_model<circle_constraint, double, 2, 1, 1>);

// Terminal constraint: must end near the origin
struct terminal_origin
{
    double tol{0.5};

    Eigen::Matrix<double, 1, 1> operator()(const Eigen::Vector2d& x) const
    {
        // h(x_N) = ||x|| - tol <= 0
        return Eigen::Matrix<double, 1, 1>{x.norm() - tol};
    }
};

static_assert(ctrlpp::terminal_constraint_model<terminal_origin, double, 2, 1>);

int main()
{
    circle_constraint g{.R = 5.0};
    terminal_origin h{.tol = 0.5};

    Eigen::Vector2d x(3.0, 4.0);
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    auto gval = g(x, u);
    auto hval = h(x);

    std::cout << "Path constraint g(x,u) = " << gval[0]
              << " (feasible: " << (gval[0] <= 0.0) << ")\n";
    std::cout << "Terminal constraint h(x) = " << hval[0]
              << " (feasible: " << (hval[0] <= 0.0) << ")\n";
}
```

## See Also

- [mpc/nmpc](../mpc/nmpc.md) -- NMPC using constraint models
- [mpc/nmhe](../mpc/nmhe.md) -- NMHE using constraint models
- [dynamics-model](dynamics-model.md) -- dynamics model concept
