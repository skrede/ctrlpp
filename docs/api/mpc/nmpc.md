# nmpc

> **Requires NLopt.** Enable with `-DCTRLPP_BUILD_NLOPT=ON` when configuring CMake.

Nonlinear Model Predictive Controller using NLP optimization with multiple shooting. Supports arbitrary nonlinear dynamics via the `dynamics_model` concept, custom stage and terminal cost overrides, nonlinear path and terminal constraints with optional soft relaxation, and automatic warm-starting via trajectory shifting.

## Header and Alias

The public `nmpc` name is the **compile-time-horizon** static controller (the
horizon `NH` is a template parameter) and is the default. Below a decision
dimension `NV` of 48 its solve is allocation-free at the solver's fixed-N floor
on the static path, in steady state and at an exact equilibrium; the default
`nmpc` row of the [RT-safety matrix](../../rt-safety-matrix.md) states the
measured coverage.
`nmpc_dynamic` is the **runtime-horizon** controller (horizon taken from
`nmpc_config::horizon`), the opt-in soft-RT path this NLopt example uses.

| Form | Header |
|------|--------|
| `nmpc<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>` (default, compile-time horizon) | `#include <ctrlpp/nmpc.h>` |
| `nmpc_dynamic<Scalar, NX, NU, Solver, Dynamics, NC, NTC>` (opt-in, runtime horizon) | `#include <ctrlpp/nmpc.h>` |

```cpp
// Default: compile-time horizon NH, allocation-free on the static path in steady state
// and at an exact equilibrium (measured coverage: the RT-safety matrix's nmpc row).
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NH,
          typename Solver,
          dynamics_model<Scalar, NX, NU> Dynamics,
          std::size_t NC = 0, std::size_t NTC = 0>
using nmpc = nmpc_static<Scalar, NX, NU, NH, Solver, Dynamics, NC, NTC>;

// Opt-in: runtime horizon from nmpc_config::horizon, soft-RT.
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
| `horizon` | `std::int32_t` | `1` | Prediction horizon length N |
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | State cost weight |
| `R` | `Matrix<Scalar, NU, NU>` | Identity | Input cost weight |
| `Qf` | `optional<Matrix<Scalar, NX, NX>>` | none | Terminal cost weight |
| `u_min` | `optional<Vector<Scalar, NU>>` | none | Element-wise lower input bound |
| `u_max` | `optional<Vector<Scalar, NU>>` | none | Element-wise upper input bound |
| `x_min` | `optional<Vector<Scalar, NX>>` | none | Element-wise lower state bound |
| `x_max` | `optional<Vector<Scalar, NX>>` | none | Element-wise upper state bound |
| `du_max` | `optional<Vector<Scalar, NU>>` | none | Maximum input rate of change per step |
| `stage_cost_cb` | `optional<function<Scalar(x, u)>>` | none | Custom stage cost override (replaces Q/R quadratic) |
| `terminal_cost_cb` | `optional<function<Scalar(x)>>` | none | Custom terminal cost override (replaces Qf quadratic) |
| `path_constraint_cb` | `optional<function<Vector<NC>(x, u)>>` | none | Path constraint g(x,u) <= 0 element-wise |
| `terminal_constraint_cb` | `optional<function<Vector<NTC>(x)>>` | none | Terminal constraint h(x_N) <= 0 element-wise |
| `soft_constraints` | `bool` | `true` | Soften constraints with slack variables |
| `path_penalty` | `Vector<Scalar, NC>` | `1e4` each | Per-constraint L1 penalty for path constraints |
| `terminal_penalty` | `Vector<Scalar, NTC>` | `1e4` each | Per-constraint L1 penalty for terminal constraints |

## Construction

### nmpc (compile-time horizon)

```cpp
using construction_error = nmpc_construction_error<SetupError>;
// = std::variant<controller_construction_error, nlp_formulation_error, SetupError>

static expected<nmpc, construction_error> create(Dynamics dynamics,
                                                 const nmpc_config<Scalar, NX, NU, NC, NTC>& config);

static expected<nmpc, construction_error> create(Dynamics dynamics,
                                                 const nmpc_config<Scalar, NX, NU, NC, NTC>& config,
                                                 Solver solver);
```

`create` is the only construction path for the compile-time-horizon controller, and so for the default `nmpc`. Its constructors are public only so `ctrlpp::expected` can build the controller in the caller's storage, and they take a passkey that only the class can make, so a configuration the controller cannot solve is a value the caller has to inspect, never a controller whose every solve fails. The horizon is the template parameter `NH`, so it carries a compile-time domain: `NH == 0` fails to compile, because a zero horizon leaves no input to apply and makes the warm-start shift and the input offset read outside the decision vector.

This is the slack-free cut: `config.horizon` must equal `NH`, and the configuration must not add slack decision variables, so the decision dimension equals the compile-time `NV = (NH+1)*NX + NH*NU`. **The shared `soft_constraints` default is `true`, so a compile-time-horizon controller with a path or terminal constraint needs `soft_constraints = false`** (or the runtime-horizon `nmpc_dynamic`, which poses the slack); `create` refuses it otherwise. A configuration without a path or terminal constraint adds no slack and is accepted with the default. Checked in order, each refusal returned before any controller exists:

| Check | Alternative | Value |
|-------|-------------|-------|
| the configuration rules of `nmpc_dynamic::create` below, from the bounds on, in the same order | `controller_construction_error` | `invalid_input_bounds` through `invalid_constraint_penalty` |
| `config.horizon != NH` | `nlp_formulation_error` | `horizon_mismatch` |
| soft constraints with a path or terminal constraint set, which would add slack | `nlp_formulation_error` | `slack_not_supported` |
| the solver refused the posed problem | `SetupError` | `argmin_setup_error`, `nlopt_setup_error`, or a custom backend's own enum |

The horizon rules of the runtime-horizon class do not apply: a non-positive or oversized `config.horizon` is a `horizon_mismatch` here. A penalty is inspected while its soft constraint is active, so an invalid penalty on a slack-adding configuration is reported as `invalid_constraint_penalty`, ahead of `slack_not_supported`. The formulation factory `detail::build_nmpc_problem_static` applies the same two formulation refusals, through the same check. The controller sets its solver up only against a problem that formulation posed: should posing refuse, `create` returns that refusal with the same alternative and value, never a controller whose solves report `setup_incomplete`.

```cpp
using controller = ctrlpp::nmpc<double, NX, NU, NH, solver, plant, NC>;

nmpc_config<double, NX, NU, NC> cfg;
cfg.horizon = NH;
cfg.path_constraint_cb = path;
cfg.soft_constraints = false;  // the default true would be refused with slack_not_supported

auto created = controller::create(plant{}, cfg);
if (!created)
{
    if (const auto* refusal = std::get_if<ctrlpp::nlp_formulation_error>(&created.error()))
        return handle_formulation(*refusal);
}
```

The controller is built once, in place, in the returned result: the two-argument form value-initializes the solver directly in the controller, and the three-argument form moves the caller's solver into it through an rvalue reference, so construction holds no second controller. The solver-taking form still takes its parameter by value, so a whole solver sits on the caller's stack as the argument until it is moved in. This matters for a solver that holds its workspace inline, such as the fixed-dimension argmin solver, which is about 41 KB at `NV = 17` on a 32-bit target; a solver with non-default settings has to come through the three-argument form and pays for it. Measured with `-fstack-usage` on the board demo's configuration (`NV = 17`, controller object 45,104 bytes on x86-64 and 44,104 on the Cortex-M7), built into static storage as the board image does, g++ 16.2.1 and arm-none-eabi-g++ 16.2.0 at `-O0`, `-Og`, `-Os` and `-O2`: no frame on the construction chain through `create` holds a controller, the out-of-line chain is between 1,400 and 3,904 bytes, and a demo built as a local takes a frame that exceeds the controller object by at most 96 bytes.

To keep the result in place, initialize the holding object directly from `create`'s prvalue (a member initializer or `auto created = controller::create(...)`); guaranteed elision then builds the controller in it. Moving the controller out of the result afterwards moves its inline solver workspace.

### nmpc_dynamic (runtime horizon)

```cpp
using construction_error = nmpc_construction_error<SetupError>;
// = std::variant<controller_construction_error, nlp_formulation_error, SetupError>

static expected<nmpc_dynamic, construction_error> create(Dynamics dynamics,
                                                         const nmpc_config<Scalar, NX, NU, NC, NTC>& config);

static expected<nmpc_dynamic, construction_error> create(Dynamics dynamics,
                                                         const nmpc_config<Scalar, NX, NU, NC, NTC>& config,
                                                         Solver solver);
```

`create` is the only construction path for the runtime-horizon controller. Its constructor is public only so `ctrlpp::expected` can build the controller in the caller's storage, and it takes a passkey that only the class can make, so a rejected configuration or a refused solver setup is a value the caller has to inspect and never a controller that quietly stands in for one. The two-argument form default-constructs the solver and chains into the three-argument form; both move the dynamics and the solver in before the NLP is posed.

`SetupError` is the error type of the bound solver's own `setup`, so the backend's cause of a setup failure reaches the caller unflattened. The error is one of three alternatives:

| Alternative | Carries |
|-------------|---------|
| `controller_construction_error` | a configuration refusal, in the order of the table below |
| `nlp_formulation_error` | never produced by `nmpc_dynamic`, whose runtime-horizon formulation has no refusal; it is there so both nonlinear MPC classes share one error shape |
| `SetupError` | the solver refused the posed problem: `argmin_setup_error` for `argmin_solver`, `nlopt_setup_error` for `nlopt_solver`, or a custom backend's own enum |

The configuration is checked first, then the controller is built in place and its solver is set up against the posed problem; a refused setup turns the result into the `SetupError` alternative, so no controller exists whose every solve would fail. The three alternative types must differ (a `static_assert` says so), because a variant whose alternatives repeat a type cannot be built from a value of that type.

```cpp
auto created = ctrlpp::nmpc_dynamic<double, NX, NU, ctrlpp::nlopt_solver<double>, Dynamics>::create(
    dynamics, cfg);
if (!created)
{
    if (const auto* refusal = std::get_if<ctrlpp::controller_construction_error>(&created.error()))
        return handle_configuration(*refusal);
    if (const auto* cause = std::get_if<ctrlpp::nlopt_setup_error>(&created.error()))
        return handle_setup(*cause);
}
auto controller = *std::move(created);
```

The configuration is validated once here, before any dimension product is formed and before any storage is reserved. Rejections, checked in order:

| Condition | `controller_construction_error` |
|-----------|----------------------------------|
| `config.horizon <= 0` | `non_positive_horizon` |
| `config.horizon` above the representable bound of the derived dimensions | `horizon_overflow` |
| `u_min` or `u_max` holds a NaN, `u_min` holds `+inf`, `u_max` holds `-inf`, or `u_min > u_max` in some entry | `invalid_input_bounds` |
| the same on `x_min` and `x_max` | `invalid_state_bounds` |
| `du_max` holds a NaN or a negative entry | `invalid_rate_bound` |
| `Q`, then `R`, then `Qf` when set, holds a NaN or an infinity | `non_finite_cost_weight` |
| that weight is asymmetric beyond rounding drift | `asymmetric_cost_weight` |
| that weight is not positive semidefinite | `indefinite_cost_weight` |
| `path_penalty` or `terminal_penalty` of an active soft constraint holds a non-finite or non-positive entry | `invalid_constraint_penalty` |

Each weight is checked completely (finite, symmetric, semidefinite) before the next one, so an indefinite `Q` is reported ahead of a non-finite `R`.

Equal bounds pin a variable and are accepted, as is an infinite bound in its open direction (`u_min = -inf`, `u_max = +inf`). A rate bound of zero, `-0.0` or `+inf` is accepted.

The weights follow the library's symmetry convention. An asymmetry within the rounding bound `N * eps * max|W_ij|` (never below `N * denorm_min`) of an `N x N` weight `W` is drift: it is symmetrized in the controller's own copy of the configuration, before the formulation captures it, so the posed cost and its analytic gradient `Q (x - x_ref)` use one symmetric weight. Without that, the gradient would not be the gradient of the cost, which depends only on the symmetric part. A larger asymmetry is structure and is refused, never silently symmetrized. Semidefiniteness is decided by the symmetric eigensolver with the floor `N * eps * max|W_ij|`, the order of its backward error. `R` is required positive semidefinite, not definite: no transcription needs it nonsingular. The weights are validated even when a custom `stage_cost_cb` or `terminal_cost_cb` replaces the quadratic cost, because they are configuration fields with valid defaults.

A penalty is inspected only while the soft constraint it weights is active, that is, when `soft_constraints` is true, the constraint dimension is nonzero and the matching constraint callable is set. Otherwise it weights nothing and is accepted whatever it holds.

`horizon` stays a signed `std::int32_t` deliberately: a mistaken negative value remains representable as negative and is therefore rejectable, whereas an unsigned field would turn the same mistake into an enormous allocation. The overflow bound is a representability condition on `std::int32_t`, not a chosen ceiling. At the worst-case configuration the horizon `N` scales the decision vector as `N*(NX + NU + NC) + NX + NTC` and the constraint rows as `N*(NX + 2*NU + NC) + NX + NTC`, so both stay representable exactly when `horizon <= (INT32_MAX - (NX + NTC)) / (NX + 2*NU + NC)`.

Two alternatives are deliberately not implemented. Validating at the first solve would surface a configuration error at the first control step, the worst possible moment. Clamping the horizon to one would turn a caller mistake into a silently different controller.

## Failure contract

`nmpc` shares the exact soft-constraint and failure contract documented for
[mpc](mpc.md). Every `solve` overload returns
`ctrlpp::expected<solve_output<Scalar, NU>, solver_error>`:

- The **success branch** holds a `solve_output` whose `input` field is the first
  control input and whose `status` field is a soft `solve_result_status`
  (`converged`, `solved_inaccurate`, or `budget_exhausted`). A budget-limited
  solve still returns its best iterate tagged `budget_exhausted`.
- The **error branch** holds a `solver_error` (`infeasible`, `invalid_problem`,
  `setup_incomplete`, or `invalid_backend_result`) and no input.

`invalid_backend_result` covers the case where the solver reports a status the
controller accepts and then returns a result the controller cannot use: a
decision vector shorter than the problem dimension, a non-finite value in the
consumed decision vector, a non-finite objective, solve time or primal residual,
or a negative iteration count. Both nonlinear controllers check all of these
before consuming the result, so the warm-start shift and the first-input slice,
whose offsets are derived from the horizon, cannot read past the end of the
solver's own storage, and no NaN or infinity becomes a command. The
runtime-horizon controller compares the length against the dimension it derived
at construction from its horizon and its state, input, path-slack and
terminal-slack contributions; the compile-time-horizon controller compares it
against `NV`. A longer-than-required result is accepted, and only its first `NV`
entries are consumed.

Nothing from a refused result is committed. The warm start of the next solve, the
previous input the rate constraints are posed against, `last_solution()` and
`trajectory()` all keep the last accepted solve, and `diagnostics().status` reads
`solve_status::invalid_backend_result`. On the compile-time-horizon controller, a
solver that writes into a caller buffer writes into a second, pre-sized candidate
buffer, and the controller exchanges the two buffers only when the result is
accepted, so the commit allocates nothing; a solver that returns its decision
vector by value leaves it in the returned result until it is accepted, and only
then are its first `NV` entries copied.

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
the first input and a soft status, or the error branch on failure. A solver setup
failure is a construction error returned by `create`, so a created controller is
always set up; only a copy whose re-posing of its problem or re-setup of its forked
solver was refused returns the error branch with `solver_error::setup_incomplete`
from every `solve` overload.

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

`stop_criterion` says which stopping test ended the last solve, for backends that report it (`argmin_solver` does; the others leave it `unreported`). `status` folds every successful stop into `solve_status::optimal`, which is enough for acting on the answer but not for saying how close the answer is to the exact optimum. Only a stop on `nlp_stop_criterion::stationarity` means the first-order optimality residual fell below its threshold, and a bound on the distance to the optimum can be built from that. A stop on the change in the objective (`objective_change`) or in the iterate (`step_size`) certifies no such distance.

| Value | The solve stopped because |
|-------|---------------------------|
| `unreported` | the backend does not say |
| `stationarity` | the KKT residual fell below `argmin_settings::kkt_tol` |
| `objective_change` | the relative change in the objective fell below `ftol_rel` |
| `step_size` | the relative change in the iterate fell below `xtol_rel` |
| `stalled` | the iterate or the objective stopped making progress, including a line search that accepted no step |
| `iteration_limit` | `max_eval` ran out |
| `time_limit` | `max_time` ran out (only under the wall-clock budget, below) |
| `other` | any other stop, including a failed solve |

### argmin backend settings

`argmin_solver` takes its `settings_type`. Under the default budget, `argmin_step_budget`, that is `argmin_settings<Scalar>`:

| Field | Default | Meaning |
|-------|---------|---------|
| `ftol_rel` | `1e-6` | relative objective-change threshold |
| `xtol_rel` | `1e-6` | relative step threshold |
| `max_eval` | `500` | iteration cap |
| `constraint_tol` | `1e-8` | constraint and feasibility tolerance |
| `warm_start` | `curvature` | what a re-solve keeps from the previous one |
| `kkt_tol` | `1e-5` | threshold on the composite KKT residual, the largest of the stationarity and feasibility legs in the infinity norm |

The third template parameter, `Constrained` (default `true`), selects the bridge. With `false` the bridge handles
variable bounds only, and its `setup` refuses a problem with general constraint rows with
`argmin_setup_error::constraints_unsupported` instead of solving it without them. The transcription always poses the
dynamics as equality rows, so `Constrained = false` cannot back an NMPC: `create` refuses it with the
`argmin_setup_error` alternative holding `constraints_unsupported`. Controllers use the default.

### argmin setup refusals

`argmin_solver`'s `setup` refuses a problem its bridge cannot hold and leaves the solver as it was. It checks in
this order and reports the first refusal:

| Refusal | The problem |
|---------|-------------|
| `constraints_unsupported` | poses general constraint rows on the bridge built with `Constrained = false` |
| `dimension_mismatch` | declares an `n_vars` other than the fixed decision dimension `NV` |
| `constraint_cap_exceeded` | poses more partitioned rows than the fixed cap `MaxM` |
| `equality_rows_exceed_dimension` | poses more equality rows than decision variables, which argmin's active-set QP cannot take |

The partitioned row count is the equality rows plus one per finite inequality side, so a two-sided row counts twice.
Bounds on the decision variables do not count. With horizon `N`, the transcription poses `NX * (N + 1)` dynamics
equalities, so `MaxM = NX * (NH + 1)` covers a controller with no other rows. Every other row has one finite side: a
rate limit (`du_max`) adds `2 * NU` rows per stage, a path constraint adds `NC` per stage, and a terminal constraint
adds `NTC`, so `2 * NU * N`, `NC * N` and `NTC` in total. An undersized `MaxM` is refused at construction: `create`
returns the `argmin_setup_error` alternative holding `constraint_cap_exceeded`.

The budget is the sixth template parameter, and it is chosen at compile time. The default instantiates only the
step-budget driver, so the solver reads no clock, and its settings have no `max_time`: setting one is a compile
error, not a value the solver ignores. A wall-clock budget is opt-in:

```cpp
using timed = ctrlpp::argmin_solver<double, ctrlpp::argmin_slsqp, true, Eigen::Dynamic, Eigen::Dynamic,
                                    ctrlpp::argmin_step_and_time_budget>;

timed::settings_type settings;  // argmin_timed_settings<double>: every field above, plus max_time
settings.max_time = 0.5;        // seconds; the default 0 runs the step-budget driver
```

This differs from `nlopt_solver`, whose `nlopt_settings::max_time` is a run-time field on every instance. The two
backends' settings otherwise carry the same tolerance and cap fields, but an argmin solver has a `max_time` to set
only when it is built with `argmin_step_and_time_budget`. A wall-clock budget is not real-time safe (see the
[RT-safety matrix](../../rt-safety-matrix.md)).

`kkt_tol` is absolute, not relative to the problem's scale. Its default is the backend's own, the `pgtol` default of Byrd, Lu, Nocedal and Zhu (1995), so leaving it unset changes nothing. On a problem whose state shrinks toward zero, an absolute threshold stops each solve once the residual is small in absolute terms. The input it returns is then accurate to an absolute amount, not a relative one. A caller who needs relative accuracy at small states tightens `kkt_tol`, and pays for it in iterations:

```cpp
using solver = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, NV, MaxM>;

solver::settings_type settings;
settings.kkt_tol = 1e-9;

auto controller = ctrlpp::nmpc<double, NX, NU, NH, solver, plant>::create(plant{}, config, solver{settings});
if(!controller)
    return handle_construction(controller.error());

auto u = controller->solve(x);
if(u && controller->diagnostics().stop_criterion == ctrlpp::nlp_stop_criterion::stationarity)
{
    // The KKT residual is below 1e-9 at this answer.
}
```

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
                                        pendulum_dynamics>::create(dynamics, cfg);
    if(!created)
    {
        // created.error() is the construction_error variant: a configuration
        // refusal or the solver's own nlopt_setup_error.
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
