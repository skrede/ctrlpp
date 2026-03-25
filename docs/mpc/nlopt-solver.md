# nlopt_solver

> **Requires NLopt.** Enable with `-DCTRLPP_BUILD_NLOPT=ON` when configuring CMake.

NLopt-based NLP solver backend for nonlinear MPC and NMHE. Satisfies the `nlp_solver` concept and can be injected as the `Solver` template parameter of `nmpc` and `nmhe`. Supports multiple algorithms (SLSQP, MMA, COBYLA, ISRES), automatic partitioning of equality and inequality constraints, and finite-difference gradient computation for constraint Jacobians.

## Header

| Form | Header |
|------|--------|
| `nlopt_solver<Scalar>` | `#include <ctrlpp/mpc/nlopt_solver.h>` |

```cpp
template <typename Scalar>
class nlopt_solver;
```

`Scalar` must be `double` (NLopt operates in double precision only).

## nlp_solver Concept

`nlopt_solver` satisfies the `nlp_solver` concept defined in `<ctrlpp/mpc/nlp_solver.h>`:

```cpp
template <typename S>
concept nlp_solver = requires { typename S::scalar_type; }
    && requires(S solver,
                const nlp_problem<typename S::scalar_type>& prob,
                const nlp_update<typename S::scalar_type>& upd) {
        { solver.setup(prob) } -> std::same_as<void>;
        { solver.solve(upd) } -> std::same_as<nlp_result<typename S::scalar_type>>;
    };
```

## Supporting Types

### nlopt_algorithm

```cpp
enum class nlopt_algorithm : std::uint8_t {
    slsqp,   // Sequential Least Squares Programming (gradient-based, supports equality)
    mma,     // Method of Moving Asymptotes (gradient-based, NO equality constraints)
    cobyla,  // Constrained Optimization BY Linear Approximations (derivative-free)
    isres    // Improved Stochastic Ranking Evolution Strategy (global, derivative-free)
};
```

### nlopt_settings

```cpp
template <typename Scalar>
struct nlopt_settings {
    nlopt_algorithm algorithm{nlopt_algorithm::slsqp};
    Scalar ftol_rel{1e-6};
    Scalar xtol_rel{1e-6};
    int max_eval{500};
    Scalar max_time{0};          // 0 = no time limit
    Scalar constraint_tol{1e-8};
};
```

## Constructors

```cpp
explicit nlopt_solver(nlopt_settings<Scalar> settings = {});
```

Constructs the solver with the given settings. Defaults to SLSQP with standard tolerances.

## Methods

### setup

```cpp
void setup(const nlp_problem<Scalar>& problem);
```

Configures the NLopt optimizer from an NLP problem definition. Partitions constraints into equality and inequality groups automatically. Throws `std::invalid_argument` if MMA is selected with equality constraints.

### solve

```cpp
auto solve(const nlp_update<Scalar>& update) -> nlp_result<Scalar>;
```

Solves the NLP from the initial guess in `update.x0`. Returns the solution, objective value, solver status, iteration count, solve time, and maximum constraint violation.

## Usage Example

```cpp
#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

struct spring_mass
{
    double dt{0.05};
    double k{1.0};
    double b{0.1};

    Eigen::Vector2d operator()(const Eigen::Vector2d& x,
                               const Eigen::Matrix<double, 1, 1>& u) const
    {
        double pos = x[0];
        double vel = x[1];
        double acc = -k * pos - b * vel + u[0];
        return {pos + vel * dt, vel + acc * dt};
    }
};

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;

    spring_mass dynamics;

    ctrlpp::nlopt_settings<double> solver_settings{
        .algorithm = ctrlpp::nlopt_algorithm::slsqp,
        .ftol_rel = 1e-8,
        .max_eval = 1000};

    ctrlpp::nmpc_config<double, NX, NU> cfg{
        .horizon = 20,
        .Q = Eigen::Vector2d(5.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity() * 0.1,
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-2.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(2.0)};

    ctrlpp::nmpc<double, NX, NU, ctrlpp::nlopt_solver<double>, spring_mass>
        controller(dynamics, cfg);

    Eigen::Vector2d x(2.0, 0.0);
    Eigen::Vector2d x_ref(0.0, 0.0);

    for(int k = 0; k < 80; ++k)
    {
        auto u_opt = controller.solve(x, x_ref);
        if(!u_opt)
            break;

        std::cout << "k=" << k << "  x=[" << x.transpose()
                  << "]  u=" << (*u_opt)[0] << "\n";
        x = dynamics(x, *u_opt);
    }
}
```

## See Also

- [nmpc](nmpc.md) -- nonlinear MPC using NLopt
- [nmhe](nmhe.md) -- nonlinear MHE using NLopt
- [osqp-solver](osqp-solver.md) -- OSQP QP solver backend
- [guides/mpc/solver-injection](../guides/mpc/solver-injection.md) -- solver injection guide
- [reference/mpc-theory](../reference/mpc-theory.md) -- NLP formulation and theory
