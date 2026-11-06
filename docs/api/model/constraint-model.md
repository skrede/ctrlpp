# constraint_model

C++20 concepts describing nonlinear inequality constraints: path constraints `g(x, u) <= 0` and terminal constraints `h(x_N) <= 0`. No library template requires these concepts. They state the shape of the callables `nmpc_config` accepts in `path_constraint_cb` and `terminal_constraint_cb`, which are stored as `std::function`. NMHE's `path_constraint_cb` takes the state only, `g(x) <= 0`, so it does not have the `constraint_model` shape.

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

Represents NC nonlinear path inequality constraints `g(x, u) <= 0`, the shape of NMPC's per-stage path constraint.

### terminal_constraint_model

```cpp
template <typename H, typename Scalar,
          std::size_t NX, std::size_t NTC>
concept terminal_constraint_model = requires(const H& h,
                                             const Vector<Scalar, NX>& x) {
    { h(x) } -> std::convertible_to<Vector<Scalar, NTC>>;
};
```

Represents NTC nonlinear terminal inequality constraints `h(x_N) <= 0`, the shape of NMPC's terminal constraint.

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

- [mpc/nmpc](../mpc/nmpc.md)<br/> NMPC using constraint models
- [mpc/nmhe](../mpc/nmhe.md)<br/> NMHE, whose path constraint takes the state only
- [dynamics-model](dynamics-model.md)<br/> dynamics model concept
