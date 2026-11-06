# argmin_solver

> **Requires argmin.** Enable with `-DCTRLPP_BUILD_ARGMIN=ON` when configuring CMake. argmin is header-only.

Header-only NLP solver backend for the nonlinear controllers and estimators, built on argmin's SQP, L-BFGS-B and
augmented-Lagrangian policies. With a runtime decision dimension it satisfies `nlp_solver` and backs `nmpc_dynamic` and
`nmhe`; with the decision dimension fixed at compile time it backs the default, compile-time-horizon `nmpc`, and with a
fixed constraint cap as well its steady-state solve allocates nothing (see the
[RT-safety matrix](../../rt-safety-matrix.md) for the measured coverage).

| Form | Header |
|------|--------|
| `argmin_solver<Scalar, Policy, Constrained, NV, MaxM, Budget>` | `#include <ctrlpp/mpc/argmin_solver.h>` |
| the policy tags, `argmin_settings`, `argmin_timed_settings`, `warm_start_mode`, the budget tags | `#include <ctrlpp/mpc/argmin_policies.h>` |
| the bridges `argmin_problem` and `argmin_constrained_problem` | `#include <ctrlpp/mpc/argmin_problem.h>` |
| `argmin_setup_error` | `#include <ctrlpp/mpc/nlp_types.h>` |

`argmin_solver.h` includes the other three.

## Template Parameters

```cpp
template <typename Scalar, typename Policy, bool Constrained = true, int NV = Eigen::Dynamic,
          int MaxM = Eigen::Dynamic, argmin_budget Budget = argmin_step_budget>
class argmin_solver;
```

| Parameter | Description |
|-----------|-------------|
| `Scalar` | Floating-point type of the problem; the policies ctrlpp names are argmin's `double` policies |
| `Policy` | A policy tag from `argmin_policies.h`, below |
| `Constrained` | `true` binds `argmin_constrained_problem`, which passes the constraint rows to argmin; `false` binds `argmin_problem`, which passes only the cost and the variable bounds |
| `NV` | Decision dimension. `Eigen::Dynamic` binds `nlp_problem<Scalar>`; a fixed `NV` binds `nlp_problem_static<Scalar, NV>` and rebinds the policy to fixed-size storage |
| `MaxM` | Cap on the partitioned constraint rows. A fixed cap makes argmin's multiplier storage inline |
| `Budget` | `argmin_step_budget` (default) or `argmin_step_and_time_budget` |

`problem_type` names the problem the solver sets up against: `nlp_problem<Scalar>` for `NV = Eigen::Dynamic`,
`nlp_problem_static<Scalar, NV>` otherwise. The compile-time-horizon `nmpc` needs the second, with
`NV = (NH + 1) * NX + NH * NU`.

### Policies

| Tag | argmin policy | Constraints it takes |
|-----|---------------|----------------------|
| `argmin_slsqp` | `kraft_slsqp_policy<>`, Kraft's SLSQP | equalities, inequalities, bounds |
| `argmin_nw_sqp` | `nw_sqp_policy<>`, the line-search SQP of Nocedal and Wright | equalities, inequalities, bounds |
| `argmin_nw_sqp_static<NV>` | `nw_sqp_policy<double, NV>` | as `argmin_nw_sqp` |
| `argmin_filter_slsqp` | `filter_slsqp_policy<>`, SLSQP with filter acceptance in place of the merit function | equalities, inequalities, bounds |
| `argmin_filter_nw_sqp` | `filter_nw_sqp_policy<>`, the Nocedal and Wright SQP with filter acceptance | equalities, inequalities, bounds |
| `argmin_lbfgsb` | `lbfgsb_policy<>`, L-BFGS-B | bounds only |
| `argmin_auglag<Inner = argmin_lbfgsb>` | `augmented_lagrangian_policy<Inner::algorithm>` | equalities and inequalities in the augmented objective, bounds to the inner solver |

A fixed `NV` rebinds every policy to `NV` through its `rebind` member, so `argmin_nw_sqp` with a fixed `NV` already
yields the type `argmin_nw_sqp_static<NV>` names. `rebind_argmin_algorithm_t<Algorithm, N>` is that mapping.
`argmin_auglag` compiles with `argmin_lbfgsb` as its inner solver; at the pinned argmin revision it does not compile
with one of the SQP tags as the inner solver.

The solver runs argmin with `argmin_ctrlpp_convergence`, a convergence policy of a gradient criterion, relative
objective and step criteria, and a stall criterion.

## Config

`settings_type` is `argmin_settings<Scalar>` under `argmin_step_budget` and `argmin_timed_settings<Scalar>` under
`argmin_step_and_time_budget`.

```cpp
template <typename Scalar>
struct argmin_settings {
    Scalar ftol_rel{1e-6};
    Scalar xtol_rel{1e-6};
    std::int32_t max_eval{500};
    Scalar constraint_tol{1e-8};
    warm_start_mode warm_start{warm_start_mode::curvature};
    Scalar kkt_tol{1e-5};
};

template <typename Scalar>
struct argmin_timed_settings : argmin_settings<Scalar> {
    Scalar max_time{0};
};
```

| Field | Default | Meaning |
|-------|---------|---------|
| `ftol_rel` | `1e-6` | relative objective-change threshold (argmin's relative objective criterion) |
| `xtol_rel` | `1e-6` | relative step threshold (argmin's relative step criterion) |
| `max_eval` | `500` | cap on argmin iterations, not on function evaluations |
| `constraint_tol` | `1e-8` | argmin's constraint and feasibility tolerance; a value of `0` or below keeps argmin's own defaults. Ignored with `Constrained = false` |
| `warm_start` | `curvature` | what a solve keeps from the previous one, below |
| `kkt_tol` | `1e-5` | threshold on argmin's gradient criterion, the composite KKT residual (the largest of the stationarity and feasibility legs) in the infinity norm; absolute, not relative. The default is the backend's own, the `pgtol` default of Byrd, Lu, Nocedal and Zhu (1995) |
| `max_time` | `0` | `argmin_timed_settings` only: wall-clock budget per solve in seconds; `0` runs the step-budget driver |

| `warm_start_mode` | A solve after the first |
|-------------------|-------------------------|
| `curvature` | restarts the policy from `update.x0` and keeps its curvature state (argmin's `reset`) |
| `primal_only` | restarts from `update.x0` and clears every algorithm state (argmin's `reset_clear`) |
| `cold` | the same as `primal_only` |

Every mode starts from `update.x0`; the controllers pass their shifted warm start there.

The budget is chosen at compile time. The default instantiates only argmin's step-budget driver, so the solver reads no
clock, and its settings have no `max_time`: setting one is a compile error, not a value the solver ignores. A
wall-clock budget is opt-in, and is not real-time safe:

```cpp
using timed = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, true, Eigen::Dynamic, Eigen::Dynamic,
                                    ctrlpp::argmin_step_and_time_budget>;

timed::settings_type settings;
settings.max_time = 0.5;
```

Under `argmin_step_and_time_budget` the driver is chosen once, at the first solve: the timed driver when `max_time` is
positive, the step-budget driver otherwise. `nlopt_solver` differs here: its `max_time` is a run-time field on every
instance.

## Construction

```cpp
explicit argmin_solver();
explicit argmin_solver(settings_type settings);
```

Construction cannot fail and does not touch argmin; the argmin driver is built at the first `solve` or `step` after a
`setup`. The solver is copyable and movable. A move keeps the driver. A copy binds a fresh bridge to the same problem and
leaves the driver empty, so the copy's next solve builds a new one and starts without the source's curvature.

Inject the solver as the controller's `Solver` parameter. The two-argument `create` value-initializes it; the
three-argument form takes a configured instance:

```cpp
using solver = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, NV, MaxM>;

solver::settings_type settings;
settings.max_eval = 50;

auto created = ctrlpp::nmpc<double, NX, NU, NH, solver, plant>::create(plant{}, config, solver{settings});
```

## Rejections

`setup` refuses a problem its bridge cannot hold, reports the first refusal and changes nothing. It checks in this
order:

| Condition | `argmin_setup_error` |
|-----------|----------------------|
| `Constrained = false` and `n_constraints > 0` | `constraints_unsupported` |
| a fixed `NV` and `n_vars != NV` | `dimension_mismatch` |
| a fixed `MaxM` and more partitioned rows than `MaxM` | `constraint_cap_exceeded` |
| more equality rows than `n_vars`, which argmin's active-set QP cannot take | `equality_rows_exceed_dimension` |

The enumerators are declared in the order `constraints_unsupported`, `constraint_cap_exceeded`, `dimension_mismatch`,
`equality_rows_exceed_dimension`. A row whose `c_lower` equals its `c_upper` is an equality; any other row is
partitioned into one inequality per finite side, so a two-sided row counts twice. Bounds on the decision variables do
not count. A controller turns a refusal into the `argmin_setup_error` alternative of its `construction_error`.

## Member Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `setup` | `expected<void, argmin_setup_error> setup(const problem_type&)` | Admits the problem, binds the bridge to it and discards any earlier driver |
| `solve` | `nlp_result<Scalar> solve(const nlp_update<Scalar>&)` | Solves from `update.x0` within the settings' budget |
| `solve_into` | `nlp_result<Scalar> solve_into(const nlp_update<Scalar>&, Eigen::Ref<Eigen::VectorX<Scalar>> solution)` | As `solve`, but writes the primal into `solution`, which must already have the problem's dimension, and leaves the result's `x` empty |
| `step` | `nlp_result<Scalar> step(const nlp_update<Scalar>&, std::int32_t max_steps)` | Restarts from `update.x0` and runs at most `max_steps` iterations |

The problem is held by reference and must outlive the solver's use of it. The bridge reads `cost_cb` and `gradient_cb`,
and `constraints_cb` with `Constrained = true`; `gradient_cb` must be set. It reads `constraint_jacobian_cb` when one is
set, and otherwise takes the constraint Jacobian by central differences with the step `cbrt(eps) * max(1, |x_j|)`.
The constrained bridge hands argmin each row in argmin's convention, satisfied where non-negative: `c - c_lower` for an
equality or a finite lower side, `c_upper - c` for a finite upper side.

The result's `status` maps argmin's status: `converged`, `ftol_reached` and `xtol_reached` to `optimal`;
`max_iterations`, `budget_exhausted` and `maxeval_reached` to `max_iterations`; `time_limit_reached` to `time_limit`;
`stalled`, `roundoff_limited`, `objective_stalled`, `line_search_failed` and `trust_region_step_rejected` to
`solved_inaccurate`; `invalid_problem`, `diverged`, `aborted` and `running` to `error`. `stop_criterion` names the test,
as the [nmpc diagnostics](nmpc.md) list it. `primal_residual` is argmin's constraint violation. `solve_time` is the
timed driver's wall time, and `0` under the step-budget driver, which reads no clock.

`argmin_solver` satisfies `nlp_stepper` as well as `nlp_solver`; the controllers call only `solve` and `solve_into`.

## Usage Example

```cpp
#include <ctrlpp/nmpc.h>
#include <ctrlpp/mpc/argmin_solver.h>

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <variant>

struct pendulum
{
    Eigen::Vector2d operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const
    {
        constexpr double dt = 0.05;
        return {x[0] + dt * x[1], x[1] + dt * (9.81 * std::sin(x[0]) + u[0])};
    }
};

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NH = 20;
    constexpr int NV = static_cast<int>((NH + 1) * NX + NH * NU);
    constexpr int MaxM = static_cast<int>(NX * (NH + 1));

    using solver = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, NV, MaxM>;
    using controller = ctrlpp::nmpc<double, NX, NU, NH, solver, pendulum>;

    ctrlpp::nmpc_config<double, NX, NU> cfg{
        .horizon = NH,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Constant(0.01),
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-5.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(5.0)};

    auto created = controller::create(pendulum{}, cfg);
    if(!created)
    {
        if(std::holds_alternative<ctrlpp::argmin_setup_error>(created.error()))
            std::cerr << "argmin refused the posed problem\n";
        return 1;
    }

    Eigen::Vector2d x(1.0, 0.0);
    for(int k = 0; k < 60; ++k)
    {
        auto u = created->solve(x);
        if(!u)
            return 1;
        std::cout << k << "," << x[0] << "," << u->input[0] << "\n";
        x = pendulum{}(x, u->input);
    }
}
```

## Limitations

### argmin_auglag penalizes the inequality rows as well as the equalities

With its default inner solver, `argmin_lbfgsb`, which takes only bounds, argmin's augmented Lagrangian folds every
equality row and every inequality row into the augmented objective and passes only the variable bounds to L-BFGS-B. NLopt's `AUGLAG_EQ`, behind
`nlopt_algorithm::auglag_mma` and `auglag_ccsaq`, penalizes only the equalities and passes the inequalities to its
inner MMA or CCSAQ solver. argmin itself switches to that equality-only scheme when the inner solver tracks
inequalities in its own state, as its moving-asymptotes family does, but no ctrlpp tag names such an inner solver. On
the pendulum of the usage example (60 steps, default settings), `argmin_auglag<>` with a runtime dimension ended 56 of
the 60 solves on `max_eval` (`budget_exhausted`), where `argmin_nw_sqp` converged on all 60 in 1 to 23 iterations; both
brought the pendulum to about 0.145 rad.

### A bound-only policy drops the constraint rows

`argmin_lbfgsb` with `Constrained = true` sets up against a problem with constraint rows and then solves it without
them, reporting a constraint violation of `0`. On an NMPC this drops the initial-state and continuity rows, so the
solve converges to a decision vector that does not follow the dynamics: on `x+ = [x_1 + 0.1 x_2, x_2 + 0.1 u]` from
`[1, 0]` it returned the input `0` and a predicted `x_1` of `[0, 0]`. Use `argmin_lbfgsb` only on a problem with
bounds alone, or inside `argmin_auglag`.

### The adapter checks dimensions only at setup

`setup` checks the shapes the table above lists and nothing more: the lengths of `x_lower`, `x_upper`, `c_lower` and
`c_upper` are not compared with `n_vars` and `n_constraints`. `solve`, `solve_into` and `step` do not check the length
of `update.x0`, nor, for `solve_into`, that of `solution`. A `solve` before the first successful `setup` reads a problem
that was never bound and is undefined behavior; it does not report `solve_status::error`. A `setup` refused after an
earlier successful one leaves the earlier problem posed, and the next `solve` answers it. The controllers never reach
these cases: they set up at construction, refuse a controller whose setup was refused, and pass vectors of the posed
dimension.

## See Also

- [nmpc](nmpc.md)<br/> nonlinear MPC
- [nmhe](nmhe.md)<br/> nonlinear MHE
- [nlopt-solver](nlopt-solver.md)<br/> NLopt NLP solver backend
- [argmin-qp-solver](argmin-qp-solver.md)<br/> argmin's QP backend for the linear controllers
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
