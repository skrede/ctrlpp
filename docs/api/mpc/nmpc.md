# nmpc

> **Requires NLopt.** Enable with `-DCTRLPP_BUILD_NLOPT=ON` when configuring CMake.

Nonlinear Model Predictive Controller using NLP optimization with multiple shooting. Supports arbitrary nonlinear dynamics via the `dynamics_model` concept, custom stage and terminal cost overrides, nonlinear path and terminal constraints with optional soft relaxation, and automatic warm-starting via trajectory shifting.

## Header and Alias

The public `nmpc` name is the **compile-time-horizon** static controller (the
horizon `NH` is a template parameter), which is allocation-free at the solver's
fixed-N floor and is the default. `nmpc_dynamic` is the **runtime-horizon**
controller (horizon taken from `nmpc_config::horizon`) — the opt-in soft-RT path
this NLopt example uses.

| Form | Header |
|------|--------|
| `nmpc<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>` (default, compile-time horizon) | `#include <ctrlpp/nmpc.h>` |
| `nmpc_dynamic<Scalar, NX, NU, Solver, Dynamics, NC, NTC>` (opt-in, runtime horizon) | `#include <ctrlpp/nmpc.h>` |

```cpp
// Default: compile-time horizon NH — allocation-free at the solver's fixed-N floor.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NH,
          typename Solver,
          dynamics_model<Scalar, NX, NU> Dynamics,
          std::size_t NC = 0, std::size_t NTC = 0>
using nmpc = nmpc_static<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>;

// Opt-in: runtime horizon from nmpc_config::horizon — soft-RT.
template <typename Scalar, std::size_t NX, std::size_t NU,
          nlp_solver Solver,
          dynamics_model<Scalar, NX, NU> Dynamics,
          std::size_t NC = 0, std::size_t NTC = 0>
class nmpc_dynamic;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type (`double`, `float`) |
| `NX` | `>= 1` | State dimension |
| `NU` | `>= 1` | Input dimension |
| `Solver` | satisfies `nlp_solver` | NLP solver backend. See [nlopt-solver](nlopt-solver.md). |
| `Dynamics` | satisfies `dynamics_model<Scalar, NX, NU>` | Discrete-time dynamics callable. See [dynamics-model](../model/dynamics-model.md). |
| `NC` | `>= 0` | Number of path constraints (default 0) |
| `NTC` | `>= 0` | Number of terminal constraints (default 0) |

## nmpc_config

Configuration struct `nmpc_config<Scalar, NX, NU, NC, NTC>` passed at construction.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `horizon` | `int` | `1` | Prediction horizon length N |
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | State cost weight |
| `R` | `Matrix<Scalar, NU, NU>` | Identity | Input cost weight |
| `Qf` | `optional<Matrix<Scalar, NX, NX>>` | none | Terminal cost weight |
| `u_min` | `optional<Vector<Scalar, NU>>` | none | Element-wise lower input bound |
| `u_max` | `optional<Vector<Scalar, NU>>` | none | Element-wise upper input bound |
| `x_min` | `optional<Vector<Scalar, NX>>` | none | Element-wise lower state bound |
| `x_max` | `optional<Vector<Scalar, NX>>` | none | Element-wise upper state bound |
| `du_max` | `optional<Vector<Scalar, NU>>` | none | Maximum input rate of change per step |
| `stage_cost` | `optional<function<Scalar(x, u)>>` | none | Custom stage cost override (replaces Q/R quadratic) |
| `terminal_cost` | `optional<function<Scalar(x)>>` | none | Custom terminal cost override (replaces Qf quadratic) |
| `path_constraint` | `optional<function<Vector<NC>(x, u)>>` | none | Path constraint g(x,u) <= 0 element-wise |
| `terminal_constraint` | `optional<function<Vector<NTC>(x)>>` | none | Terminal constraint h(x_N) <= 0 element-wise |
| `soft_constraints` | `bool` | `true` | Soften constraints with slack variables |
| `path_penalty` | `Vector<Scalar, NC>` | `1e4` each | Per-constraint L1 penalty for path constraints |
| `terminal_penalty` | `Vector<Scalar, NTC>` | `1e4` each | Per-constraint L1 penalty for terminal constraints |

## Construction

### nmpc (compile-time horizon)

```cpp
nmpc(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config);

nmpc(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config,
     Solver solver);
```

Constructs the controller from a dynamics model and configuration. Builds the NLP formulation and initializes the solver. The horizon is the template parameter `NH`, so it carries a compile-time domain: `NH == 0` fails to compile, because a zero horizon leaves no input to apply and makes the warm-start shift and the input offset read outside the decision vector.

`config.horizon` must equal `NH`, and the configuration must be slack-free (`soft_constraints` off, or no path/terminal constraint set), so the runtime decision dimension equals the compile-time `NV = (NH+1)*NX + NH*NU`. Both preconditions are checked unconditionally, in a release build as well as a debug build, before any of the dependent problem data is built. A violation is reported through the controller's existing error channel: setup is marked incomplete, so every subsequent `solve` returns `solver_error::setup_incomplete`.

Calling the formulation factory directly surfaces the reason instead:

| Condition | `nlp_formulation_error` |
|-----------|--------------------------|
| `config.horizon != NH` | `horizon_mismatch` |
| a configuration that would add slack decision variables | `slack_not_supported` |

### nmpc_dynamic (runtime horizon)

```cpp
static auto try_create(Dynamics dynamics,
                       const nmpc_config<Scalar, NX, NU, NC, NTC>& config)
    -> expected<nmpc_dynamic, controller_construction_error>;

static auto try_create(Dynamics dynamics,
                       const nmpc_config<Scalar, NX, NU, NC, NTC>& config,
                       Solver solver)
    -> expected<nmpc_dynamic, controller_construction_error>;
```

`try_create` is the construction path for the runtime-horizon controller. The two-argument form default-constructs the solver and chains into the three-argument form; both move the dynamics and the solver in before the NLP is posed.

```cpp
auto created = ctrlpp::nmpc_dynamic<double, NX, NU, Solver, Dynamics>::try_create(
    dynamics, cfg);
if (!created)
    return handle(created.error());
auto controller = *std::move(created);
```

The prediction horizon is the only runtime quantity that scales the posed problem, so it is validated once here, before any dimension product is formed and before any storage is reserved. Rejections, checked in order:

| Condition | `controller_construction_error` |
|-----------|----------------------------------|
| `config.horizon <= 0` | `non_positive_horizon` |
| `config.horizon` above the representable bound of the derived dimensions | `horizon_overflow` |

`horizon` stays a signed `int` deliberately: a mistaken negative value remains representable as negative and is therefore rejectable, whereas an unsigned field would turn the same mistake into an enormous allocation. The overflow bound is a representability condition on `int`, not a chosen ceiling. At the worst-case configuration the horizon `N` scales the decision vector as `N*(NX + NU + NC) + NX + NTC` and the constraint rows as `N*(NX + 2*NU + NC) + NX + NTC`, so both stay representable exactly when `horizon <= (INT_MAX - (NX + NTC)) / (NX + 2*NU + NC)`.

Two alternatives are deliberately not implemented. Validating at the first solve would surface a configuration error at the first control step, the worst possible moment. Clamping the horizon to one would turn a caller mistake into a silently different controller.

The matching throwing constructors remain available when the compiler has exception support (`CTRLPP_HAS_EXCEPTIONS`); they delegate to `try_create` and throw the `bad_expected_access` of the active `ctrlpp::expected` target on a rejected configuration. On an exception-free build they are compiled out and `try_create` is the only construction path.

## Failure contract

`nmpc` shares the exact soft-constraint and failure contract documented for
[mpc](mpc.md). Every `solve` overload returns
`ctrlpp::expected<solve_output<Scalar, NU>, solver_error>`:

- The **success branch** holds a `solve_output` whose `input` field is the first
  control input and whose `status` field is a soft `solve_result_status`
  (`converged`, `solved_inaccurate`, or `budget_exhausted`). A budget-limited
  solve still returns its best iterate tagged `budget_exhausted`.
- The **error branch** holds a `solver_error` (`infeasible`, `invalid_problem`, or
  `setup_incomplete`) and no input.

The input is reached explicitly through `->input`; there is no implicit conversion
to `Vector<Scalar, NU>`. On the error branch the internal previous-input record is
not updated and no hidden fallback is applied. Use `set_applied_input` to record
the input the caller actually commanded.

## Methods

### solve (regulation)

```cpp
ctrlpp::expected<solve_output<Scalar, NU>, solver_error>
solve(const Vector<Scalar, NX>& x0);
```

Solves the NLP for regulating state to the origin. Returns the success branch with
the first input and a soft status, or the error branch on failure. If solver setup
failed at construction (reported through the solver's `try_setup`), every `solve`
overload returns the error branch with `solver_error::setup_incomplete`.

### solve (constant reference)

```cpp
ctrlpp::expected<solve_output<Scalar, NU>, solver_error>
solve(const Vector<Scalar, NX>& x0, const Vector<Scalar, NX>& x_ref);
```

Solves the NLP for tracking a constant reference.

### solve (trajectory reference)

```cpp
ctrlpp::expected<solve_output<Scalar, NU>, solver_error>
solve(const Vector<Scalar, NX>& x0, std::span<const Vector<Scalar, NX>> x_ref);
```

Solves the NLP for tracking a time-varying reference trajectory. References shorter
than the horizon are back-filled from the last supplied element; an empty span is
rejected via the error branch (`solver_error::invalid_problem`).

### set_applied_input

```cpp
void set_applied_input(const Vector<Scalar, NU>& u);
```

Records the control input the caller actually commanded. The internal previous
input is updated only on a successful solve; this accessor keeps that reference
consistent when the caller commands a different input.

### trajectory

```cpp
ctrlpp::expected<std::pair<std::vector<Vector<Scalar, NX>>,
                           std::vector<Vector<Scalar, NU>>>,
                 solver_error> trajectory() const;
```

Returns the full predicted state and input trajectories from the last solve. It is
guarded: before the first valid solve it returns the error branch
(`solver_error::setup_incomplete`) rather than stale data.

### diagnostics

```cpp
mpc_diagnostics<Scalar> diagnostics() const;
```

Returns solver diagnostics including constraint violation metrics and total slack.

## Usage Example

```cpp
// gnuplot: plot "< ./nmpc_pendulum" using 1:2 with lines title "theta"
#include <ctrlpp/nmpc.h>
#include <ctrlpp/mpc/nlopt_solver.h>

#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iostream>

// Simple pendulum dynamics: x = [theta, omega]
struct pendulum_dynamics
{
    double dt{0.05};
    double g{9.81};
    double l{1.0};

    Eigen::Vector2d operator()(const Eigen::Vector2d& x,
                               const Eigen::Matrix<double, 1, 1>& u) const
    {
        double theta = x[0];
        double omega = x[1];
        // theta = 0 is the unstable upright equilibrium under this convention.
        double alpha = g / l * std::sin(theta) + u[0];
        return {theta + omega * dt, omega + alpha * dt};
    }
};

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;

    pendulum_dynamics dynamics{.dt = 0.05};

    ctrlpp::nmpc_config<double, NX, NU> cfg{
        .horizon = 30,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity() * 0.01,
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-5.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(5.0)};

    // Runtime horizon (cfg.horizon) with the NLopt solver -> the opt-in nmpc_dynamic.
    auto created = ctrlpp::nmpc_dynamic<double, NX, NU, ctrlpp::nlopt_solver<double>,
                                        pendulum_dynamics>::try_create(dynamics, cfg);
    if(!created)
    {
        // created.error() carries the controller_construction_error.
        return 1;
    }
    auto controller = *std::move(created);

    Eigen::Vector2d x(1.0, 0.0);  // Start at 1 radian
    Eigen::Vector2d x_ref(0.0, 0.0);  // Drive toward the upright vertical (theta = 0)

    for(int k = 0; k < 100; ++k)
    {
        auto u_opt = controller.solve(x, x_ref);
        if(!u_opt)
        {
            // u_opt.error() carries the solver_error describing the failure.
            std::cerr << "NMPC solve failed at step " << k << "\n";
            return EXIT_FAILURE;
        }

        std::cout << "k=" << k
                  << "  theta=" << x[0]
                  << "  omega=" << x[1]
                  << "  u=" << u_opt->input[0] << "\n";

        x = dynamics(x, u_opt->input);
    }
}
```

## See Also

- [mpc](mpc.md)<br/> linear MPC
- [nlopt-solver](nlopt-solver.md)<br/> NLopt NLP solver backend
- [dynamics-model](../model/dynamics-model.md)<br/> dynamics model concept
- [constraint-model](../model/constraint-model.md)<br/> constraint model concepts
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
- [background/mpc](../../background/mpc.md)<br/> MPC theory and background
