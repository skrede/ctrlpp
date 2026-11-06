# Concept-Based Solver Injection

ctrlpp's MPC and MHE classes accept their optimization backend as a template
parameter. `mpc` and `mhe` require the C++20 concept `qp_solver`, and
`nmpc_dynamic` and `nmhe` require `nlp_solver`. The compile-time-horizon `nmpc`
places no concept on its solver, but sets it up against the fixed-dimension
`nlp_problem_static` (see [nmpc](../../api/mpc/nmpc.md)). This design gives you
compile-time solver dispatch with no virtual call.

## The Concepts

### qp_solver

```cpp
template <typename S>
concept qp_solver = requires { typename S::scalar_type; }
    && requires(S solver, const qp_update<typename S::scalar_type>& upd) {
        { solver.solve(upd) } -> std::same_as<qp_result<typename S::scalar_type>>;
    }
    && requires(S solver, const qp_problem<typename S::scalar_type>& prob) {
        { solver.setup(prob).has_value() } -> std::convertible_to<bool>;
    };
```

A QP solver stores problem structure in `setup()` and solves updated instances
in `solve()`. The separation allows warm-starting between MPC calls. `setup()`
is fallible: it returns a `ctrlpp::expected<void, E>` over the backend's own
setup-error enum, and a backend with no setup failure mode returns a trivially
succeeding one.

### nlp_solver

```cpp
template <typename S>
concept nlp_solver = requires { typename S::scalar_type; }
    && requires(S solver, const nlp_update<typename S::scalar_type>& upd) {
        { solver.solve(upd) } -> std::same_as<nlp_result<typename S::scalar_type>>;
    }
    && requires(S solver, const nlp_problem<typename S::scalar_type>& prob) {
        { solver.setup(prob).has_value() } -> std::convertible_to<bool>;
    };
```

The NLP interface mirrors the QP interface. The `nlp_problem` carries cost
and constraint functions rather than matrices.

## Using the Shipped Solvers

ctrlpp ships four backends, each behind its own CMake option:

| Backend | Concept | Serves |
|---------|---------|--------|
| [osqp_solver](../../api/mpc/osqp-solver.md) | `qp_solver` | `mpc`, `mhe` |
| [argmin_qp_solver](../../api/mpc/argmin-qp-solver.md) | `qp_solver` | `mpc`, `mhe` |
| [nlopt_solver](../../api/mpc/nlopt-solver.md) | `nlp_solver` | `nmpc_dynamic`, `nmhe` |
| [argmin_solver](../../api/mpc/argmin-solver.md) | `nlp_solver` with a runtime dimension | `nmpc_dynamic`, `nmhe`; with a fixed dimension, `nmpc` |

### Linear MPC with OSQP

> **Requires OSQP.** Enable with `-DCTRLPP_BUILD_OSQP=ON`.

```cpp
// Usage: ./solver_injection | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'position', '' using 1:3 with lines title 'velocity'"
#include <ctrlpp/mpc.h>
#include <ctrlpp/mpc/osqp_solver.h>
#include <ctrlpp/model/propagate.h>

#include <Eigen/Dense>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr double dt = 0.1;

    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(),
        Eigen::Matrix2d::Identity(),
        Eigen::Matrix<double, 2, 1>::Zero()};

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 20,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-1.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(1.0),
        .x_min = Eigen::Vector2d(
            -std::numeric_limits<double>::infinity(), -2.0),
        .x_max = Eigen::Vector2d(
            std::numeric_limits<double>::infinity(), 2.0)};

    auto created = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(sys, cfg);
    if (!created)
    {
        std::cerr << "MPC configuration rejected\n";
        return EXIT_FAILURE;
    }
    auto controller = *std::move(created);

    Eigen::Vector2d x(5.0, 0.0);

    for (double t = 0.0; t < 10.0; t += dt)
    {
        auto u_opt = controller.solve(x);
        if (!u_opt)
        {
            std::cerr << "Solve failed\n";
            return EXIT_FAILURE;
        }
        x = ctrlpp::propagate(sys, x, u_opt->input);
        std::cout << t << "," << x[0] << "," << x[1] << "\n";
    }
}
```

`argmin_qp_solver` (`-DCTRLPP_BUILD_ARGMIN=ON`) is injected the same way:
`ctrlpp::mpc<double, NX, NU, ctrlpp::argmin_qp_solver>`. Its setup-error type
differs from OSQP's, so the `construction_error` alternative a refused setup
lands in changes with it.

### Nonlinear MPC with NLopt

> **Requires NLopt.** Enable with `-DCTRLPP_BUILD_NLOPT=ON`.

Swap the solver template parameter and provide nonlinear dynamics. With NLopt
the controller is the runtime-horizon `nmpc_dynamic`, which takes its horizon
from `cfg.horizon`:

<!-- snippet: exceptions -->
```cpp
#include <ctrlpp/nmpc.h>
#include <ctrlpp/mpc/nlopt_solver.h>

#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <utility>
#include <iostream>

struct pendulum_dynamics
{
    double dt{0.05};

    Eigen::Vector2d operator()(const Eigen::Vector2d& x,
                               const Eigen::Matrix<double, 1, 1>& u) const
    {
        return {x[0] + x[1] * dt, x[1] + (9.81 * std::sin(x[0]) + u[0]) * dt};
    }
};

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;

    pendulum_dynamics dynamics;
    ctrlpp::nmpc_config<double, NX, NU> cfg{
        .horizon = 20,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Constant(0.01),
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-5.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(5.0)};

    auto created = ctrlpp::nmpc_dynamic<double, NX, NU, ctrlpp::nlopt_solver<double>,
                                        pendulum_dynamics>::create(dynamics, cfg);
    if (!created)
        return EXIT_FAILURE;
    auto controller = *std::move(created);

    Eigen::Vector2d x(1.0, 0.0);
    for (int k = 0; k < 60; ++k)
    {
        auto u_opt = controller.solve(x);
        if (!u_opt)
            return EXIT_FAILURE;
        x = dynamics(x, u_opt->input);
        std::cout << k << "," << x[0] << "," << x[1] << "\n";
    }
}
```

### Nonlinear MPC with argmin

> **Requires argmin.** Enable with `-DCTRLPP_BUILD_ARGMIN=ON`.

`argmin_solver` backs the same `nmpc_dynamic`, and the default
compile-time-horizon `nmpc` once its decision dimension is fixed to
`NV = (NH + 1) * NX + NH * NU`:

```cpp
using solver = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, NV>;
auto created = ctrlpp::nmpc<double, NX, NU, NH, solver, pendulum_dynamics>::create(dynamics, cfg);
```

`cfg.horizon` must equal `NH`. The [argmin_solver](../../api/mpc/argmin-solver.md)
page carries a complete program.

The rest of the API is identical: `controller.solve(x)` returns
`ctrlpp::expected<solve_output<Scalar, NU>, solver_error>`, where the success
branch carries `->input` plus a soft `status` and the error branch carries a
`solver_error`.

## Writing a Custom Solver

Any type satisfying `qp_solver` or `nlp_solver` can be injected. Here is a
skeleton for a custom QP solver:

```cpp
#include <ctrlpp/expected.h>
#include <ctrlpp/mpc/qp_solver.h>

#include <cstdint>

enum class my_setup_error : std::uint8_t
{
    rejected
};

struct my_qp_solver
{
    using scalar_type = double;

    ctrlpp::expected<void, my_setup_error> setup(const ctrlpp::qp_problem<double>& /*problem*/)
    {
        // Store the problem (P, q, A, l, u) and allocate the workspace;
        // report a failure as ctrlpp::unexpected(my_setup_error::rejected)
        return {};
    }

    ctrlpp::qp_result<double> solve(const ctrlpp::qp_update<double>& /*update*/)
    {
        ctrlpp::qp_result<double> result;
        // Solve QP with updated linear terms / bounds
        // Fill in solution vector, status, cost, timing
        return result;
    }
};

static_assert(ctrlpp::qp_solver<my_qp_solver>);
```

Then use it. The two-argument `create` default-constructs the solver; the
three-argument form takes a preconfigured instance:

```cpp
auto created = ctrlpp::mpc<double, NX, NU, my_qp_solver>::create(sys, cfg);
```

A refused `setup()` refuses `create`: `created.error()` then holds the
`my_setup_error` alternative of the controller's `construction_error`, so the
backend's own cause reaches the caller and no controller exists whose every
solve would fail.

The `static_assert` catches concept mismatches at compile time with clear
diagnostics.

## Why Concepts Instead of Virtual Dispatch

- **No virtual call**<br/>the solver type is known at compile time, so the
  compiler can inline the solve call.
- **Compile-time checking**<br/>a concept violation produces a clear error
  message pointing to the missing method.
- **No indirection**<br/>the solver is a member of the controller, not a
  pointer to an interface. Whether a solve allocates depends on the backend and
  the controller; the [RT-safety matrix](../../rt-safety-matrix.md) records it.

## Next Steps

- [OSQP Solver API](../../api/mpc/osqp-solver.md)<br/> OSQP wrapper details
- [argmin QP Solver API](../../api/mpc/argmin-qp-solver.md)<br/> argmin QP backend details
- [NLopt Solver API](../../api/mpc/nlopt-solver.md)<br/> NLopt wrapper details
- [argmin Solver API](../../api/mpc/argmin-solver.md)<br/> argmin NLP backend details
- [MPC API Reference](../../api/mpc/mpc.md)<br/> linear MPC interface
- [NMPC API Reference](../../api/mpc/nmpc.md)<br/> nonlinear MPC interface
- [MPC Theory](../../background/mpc.md)<br/> QP formulation details
- [NMPC Theory](../../background/nmpc.md)<br/> NLP formulation details
