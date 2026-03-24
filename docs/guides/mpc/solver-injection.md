# Concept-Based Solver Injection

ctrlpp's MPC and NMPC controllers accept their optimisation backend as a
template parameter. The backend must satisfy one of two C++23 concepts:
`qp_solver` for linear MPC, or `nlp_solver` for nonlinear MPC. This design
gives you compile-time solver dispatch with zero virtual overhead.

## The Concepts

### qp_solver

```cpp
template <typename S>
concept qp_solver = requires { typename S::scalar_type; }
    && requires(S solver,
                const qp_problem<typename S::scalar_type>& prob,
                const qp_update<typename S::scalar_type>& upd) {
    { solver.setup(prob) } -> std::same_as<void>;
    { solver.solve(upd) } -> std::same_as<qp_result<typename S::scalar_type>>;
};
```

A QP solver stores problem structure in `setup()` and solves updated instances
in `solve()`. The separation allows warm-starting between MPC calls.

### nlp_solver

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

The NLP interface mirrors the QP interface. The `nlp_problem` carries cost
and constraint functions rather than matrices.

## Using the Default Solvers

ctrlpp provides two solver backends: OSQP for QP and NLopt for NLP.

### Linear MPC with OSQP

> **Requires OSQP.** Enable with `-DCTRLPP_BUILD_OSQP=ON`.

```cpp
#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/propagate.h"

#include <Eigen/Dense>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr double dt = 0.1;

    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        .A = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        .B = (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(),
        .C = Eigen::Matrix2d::Identity(),
        .D = Eigen::Matrix<double, 2, 1>::Zero()};

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 20,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
        .Qf = std::nullopt,
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-1.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(1.0),
        .x_min = Eigen::Vector2d(
            -std::numeric_limits<double>::infinity(), -2.0),
        .x_max = Eigen::Vector2d(
            std::numeric_limits<double>::infinity(), 2.0)};

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> controller(sys, cfg);

    Eigen::Vector2d x(5.0, 0.0);

    for (double t = 0.0; t < 10.0; t += dt)
    {
        auto u_opt = controller.solve(x);
        if (!u_opt)
        {
            std::cerr << "Solve failed\n";
            return EXIT_FAILURE;
        }
        x = ctrlpp::propagate(sys, x, *u_opt);
        std::cout << t << "," << x[0] << "," << x[1] << "\n";
    }
}
```

### Nonlinear MPC with NLopt

> **Requires NLopt.** Enable with `-DCTRLPP_BUILD_NLOPT=ON`.

Swap the solver template parameter and provide nonlinear dynamics:

```cpp
#include "ctrlpp/nmpc.h"
#include "ctrlpp/mpc/nlopt_solver.h"

// ... define dynamics_model and nmpc_config ...

ctrlpp::nmpc<double, NX, NU, ctrlpp::nlopt_solver> controller(dynamics, cfg);
```

The rest of the API is identical: `controller.solve(x)` returns
`std::optional<Vector>`.

## Writing a Custom Solver

Any type satisfying `qp_solver` or `nlp_solver` can be injected. Here is a
skeleton for a custom QP solver:

```cpp
struct my_qp_solver
{
    using scalar_type = double;

    void setup(const ctrlpp::qp_problem<double>& problem)
    {
        // Store problem structure (H, A matrices, bounds)
        // Allocate workspace
    }

    auto solve(const ctrlpp::qp_update<double>& update)
        -> ctrlpp::qp_result<double>
    {
        // Solve QP with updated linear terms / bounds
        // Return solution vector, status, cost, timing
    }
};

static_assert(ctrlpp::qp_solver<my_qp_solver>);
```

Then use it:

```cpp
ctrlpp::mpc<double, NX, NU, my_qp_solver> controller(sys, cfg);
```

The `static_assert` catches concept mismatches at compile time with clear
diagnostics.

## Why Concepts Instead of Virtual Dispatch

- **Zero overhead** -- the solver type is known at compile time. The compiler
  inlines the solve call.
- **Compile-time checking** -- a concept violation produces a clear error
  message pointing to the missing method.
- **No heap allocation** -- the solver lives inside the MPC object with no
  pointer indirection.

## Next Steps

- [OSQP Solver API](../../mpc/osqp-solver.md) -- OSQP wrapper details
- [NLopt Solver API](../../mpc/nlopt-solver.md) -- NLopt wrapper details
- [MPC API Reference](../../mpc/mpc.md) -- linear MPC interface
- [NMPC API Reference](../../mpc/nmpc.md) -- nonlinear MPC interface
- [MPC Theory](../../reference/mpc-theory.md) -- QP/NLP formulation details
