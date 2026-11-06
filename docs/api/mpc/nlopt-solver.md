# nlopt_solver

> **Requires NLopt.** Enable with `-DCTRLPP_BUILD_NLOPT=ON` when configuring CMake.

NLopt-based NLP solver backend for nonlinear MPC and NMHE. Satisfies the `nlp_solver` concept and can be injected as the `Solver` template parameter of `nmpc_dynamic` and `nmhe`. It cannot back the compile-time-horizon `nmpc`, whose solver must set up against `nlp_problem_static` (see [argmin_solver](argmin-solver.md)). It offers seven NLopt algorithms, partitions the constraint rows into equalities and one-sided inequalities, and differentiates the constraints by central differences.

## Header

| Form | Header |
|------|--------|
| `nlopt_solver<Scalar>` | `#include <ctrlpp/mpc/nlopt_solver.h>` |
| `nlopt_algorithm`, `nlopt_settings<Scalar>` | `#include <ctrlpp/mpc/nlopt_types.h>` (also included by `nlopt_solver.h`) |

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
    && requires(S solver, const nlp_update<typename S::scalar_type>& upd) {
        { solver.solve(upd) } -> std::same_as<nlp_result<typename S::scalar_type>>;
    }
    && requires(S solver, const nlp_problem<typename S::scalar_type>& prob) {
        { solver.setup(prob).has_value() } -> std::convertible_to<bool>;
    };
```

There is one setup shape: `setup(problem)` returns a `ctrlpp::expected<void, E>` over the backend's own setup-error enum, through which the backend reports its setup failures. A backend with no setup failure mode writes a trivially succeeding fallible setup rather than an infallible one. `argmin_solver`'s setup refuses the problems its bridge cannot hold, with the `argmin_setup_error` enumerators listed on the [argmin_solver](argmin-solver.md) page.

The adapter reads `cost_cb`, `gradient_cb` and `constraints_cb` from the problem, and never `constraint_jacobian_cb`. A gradient-based algorithm (`slsqp`, `mma`, `ccsaq` and the two `auglag` variants) needs `gradient_cb`: without one the adapter leaves the gradient NLopt hands it unwritten. A row whose `c_lower` equals its `c_upper` is an equality; any other row contributes one inequality per finite side.

> **Note:** NLopt's C++ API throws by upstream design. `setup` and `solve` contain every exception that NLopt, the adapter or a problem callable raises, including one of a type not derived from `std::exception`, except a forced unwind on libstdc++ (below), so the containment lives in this adapter rather than in the controllers and estimators that consume it. In `setup`, a `std::bad_alloc` reports `nlopt_setup_error::allocation_failed` and anything else reports `nlopt_setup_error::backend_rejected`. In `solve`, any exception reports `solve_status::error`. The adapter itself still requires exception support and is exempt from the library's embedded no-exceptions floor.
>
> glibc implements thread cancellation as a forced unwind. On libstdc++ the adapter re-raises it, so a cancellation acting in `setup` or in the adapter's own evaluations (the cost re-evaluation after a round-off-limited stop, and the constraint-violation evaluation after `optimize`) unwinds the thread. A cancellation acting inside a callback that NLopt invokes during `optimize` is caught by NLopt's own C++ wrapper, which does not re-raise it, and glibc then aborts the process, so cancellation during `optimize` is not supported. On libc++ the adapter cannot recognize the forced unwind, so a cancellation there is contained like any other exception and glibc aborts. MSVC has no such cancellation.

## Supporting Types

### nlopt_algorithm

Defined in `<ctrlpp/mpc/nlopt_types.h>`.

```cpp
enum class nlopt_algorithm : std::uint8_t {
    slsqp,
    mma,
    cobyla,
    isres,
    auglag_mma,
    auglag_ccsaq,
    ccsaq
};
```

| Enumerator | NLopt algorithm | Notes |
|------------|-----------------|-------|
| `slsqp` | `LD_SLSQP` | gradient-based SQP; takes equalities and inequalities |
| `mma` | `LD_MMA` | gradient-based; refuses equality rows at setup |
| `cobyla` | `LN_COBYLA` | derivative-free |
| `isres` | `GN_ISRES` | global, derivative-free; NLopt requires every variable bound finite |
| `auglag_mma` | `AUGLAG_EQ` around `LD_MMA` | equalities in the outer penalty, inequalities and bounds to the inner solver |
| `auglag_ccsaq` | `AUGLAG_EQ` around `LD_CCSAQ` | as `auglag_mma` |
| `ccsaq` | `LD_CCSAQ` | gradient-based; refuses equality rows at setup |

For the two `auglag` variants the inner optimizer gets the same `ftol_rel`, `xtol_rel`, `max_eval` and `max_time` as the outer one.

### nlopt_setup_error

Defined in `<ctrlpp/mpc/nlp_types.h>`.

```cpp
enum class nlopt_setup_error : std::uint8_t {
    incompatible_equality_constraints,  // raw mma/ccsaq selected on a problem with equality constraints
    backend_rejected,                   // an error other than an allocation failure while building or configuring
    allocation_failed                   // an allocation failed while building or configuring
};
```

| Enumerator | Cause |
|------------|-------|
| `incompatible_equality_constraints` | Raw MMA and raw CCSAQ cannot handle equality constraints; the `auglag_mma` and `auglag_ccsaq` variants absorb them into the outer augmented-Lagrangian penalty and accept the same problem. |
| `backend_rejected` | NLopt or the adapter raised an error other than an allocation failure while building or configuring the optimizer: an argument NLopt rejects (such as a bound vector whose length differs from `n_vars`), or any other error, standard or not. No optimizer exists. |
| `allocation_failed` | An allocation failed while the optimizer was built or configured (`std::bad_alloc`, which NLopt also raises for its own out-of-memory). No optimizer exists, and a later `setup` may succeed once memory is available. |

### nlopt_settings

Defined in `<ctrlpp/mpc/nlopt_types.h>`.

```cpp
template <typename Scalar>
struct nlopt_settings {
    nlopt_algorithm algorithm{nlopt_algorithm::slsqp};
    Scalar ftol_rel{1e-6};
    Scalar xtol_rel{1e-6};
    std::int32_t max_eval{500};
    Scalar max_time{0};
    Scalar constraint_tol{1e-8};
};
```

| Field | Default | Meaning |
|-------|---------|---------|
| `algorithm` | `slsqp` | the NLopt algorithm, above |
| `ftol_rel` | `1e-6` | NLopt's relative objective tolerance (`set_ftol_rel`) |
| `xtol_rel` | `1e-6` | NLopt's relative step tolerance (`set_xtol_rel`) |
| `max_eval` | `500` | cap on objective evaluations (`set_maxeval`), not on iterations |
| `max_time` | `0` | wall-clock limit in seconds (`set_maxtime`); `0` sets none |
| `constraint_tol` | `1e-8` | the tolerance NLopt applies to every constraint row |

The iteration count the solver reports, and with it `diagnostics().iterations` on a controller, is NLopt's count of objective evaluations.

## Constructors

```cpp
explicit nlopt_solver(nlopt_settings<Scalar> settings = {});
```

Constructs the solver with the given settings. Defaults to SLSQP with standard tolerances.

The solver is copyable and movable. NLopt registers every callback with the solver's address, so a copy or a move re-runs `setup` against the same referenced problem to rebind the callbacks to the new object. A constructor or assignment cannot return an error, so when that re-setup fails (because the referenced problem changed since the source was set up) the new solver records it and its next `solve` reports `solve_status::error`. The same holds for a copy of a solver that was never set up. Every setup, including this re-setup, rebuilds the constraint partition from scratch, so a solver re-set up against a problem without constraint rows registers no constraint callbacks and evaluates no constraint, whatever problem it was set up on before.

## Methods

### setup

```cpp
ctrlpp::expected<void, nlopt_setup_error> setup(const nlp_problem<Scalar>& problem);
```

Fallible setup: configures the NLopt optimizer from an NLP problem definition and partitions constraints into equality and inequality groups automatically. Returns an empty `expected` on success, `nlopt_setup_error::incompatible_equality_constraints` if raw MMA or raw CCSAQ is selected on a problem with equality constraints, `nlopt_setup_error::allocation_failed` if an allocation fails while the optimizer is built or configured, and `nlopt_setup_error::backend_rejected` for any other error raised meanwhile, standard or not. After a failed setup the solver holds no usable optimizer, and `solve` reports `solve_status::error` until a later `setup` succeeds. This is the only setup API.

### solve

```cpp
nlp_result<Scalar> solve(const nlp_update<Scalar>& update);
```

Solves the NLP from the initial guess in `update.x0`. Returns the solution, the objective value, the status, the number of objective evaluations, the wall-clock time of the call, and the largest violation of any constraint row by its bound. `stop_criterion` stays `unreported`: NLopt's return code is folded into `status` (`SUCCESS`, `STOPVAL_REACHED`, `FTOL_REACHED` and `XTOL_REACHED` into `optimal`, `MAXEVAL_REACHED` into `max_iterations`, `MAXTIME_REACHED` into `time_limit`, `ROUNDOFF_LIMITED` into `solved_inaccurate`, any other into `error`).

NLopt requires the starting point to lie inside the variable bounds, so `solve` first projects each component of `update.x0` onto `[x_lower[i], x_upper[i]]`. A consistent box that excludes the previous iterate (for example an input box that excludes the controllers' initial zero warm start) is therefore solved rather than refused. Infinite bounds leave a component unchanged, and inverted or NaN bounds are left for NLopt to refuse, which `solve` reports as an error.

A round-off-limited stop reports `solve_status::solved_inaccurate` with the objective re-evaluated at the returned point. Every other exception raised by NLopt or by a problem callable during the solve, of any type and including one raised by that re-evaluation or by the constraint evaluation behind `primal_residual`, and including a start vector whose length differs from `n_vars`, reports `solve_status::error` with an empty `x` and zero objective, iteration count and violation; the controllers map that status to `solver_error::invalid_problem`.

## Usage Example

<!-- snippet: exceptions -->
```cpp
// gnuplot: plot "< ./nlopt_demo" using 1:3 with lines title "state"
#include <ctrlpp/nmpc.h>
#include <ctrlpp/mpc/nlopt_solver.h>

#include <Eigen/Dense>

#include <cmath>
#include <utility>
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

    auto created = ctrlpp::nmpc_dynamic<double, NX, NU, ctrlpp::nlopt_solver<double>, spring_mass>::create(
        dynamics, cfg, ctrlpp::nlopt_solver<double>{solver_settings});
    if(!created)
        return 1;
    auto controller = *std::move(created);

    Eigen::Vector2d x(2.0, 0.0);
    Eigen::Vector2d x_ref(0.0, 0.0);

    for(int k = 0; k < 80; ++k)
    {
        auto u_opt = controller.solve(x, x_ref);
        if(!u_opt)
            break;

        std::cout << "k=" << k << "  x=[" << x.transpose()
                  << "]  u=" << u_opt->input[0] << "\n";
        x = dynamics(x, u_opt->input);
    }
}
```

## Limitations

### The auglag variants penalize only the equality rows

`auglag_mma` and `auglag_ccsaq` run NLopt's `AUGLAG_EQ`: the outer augmented Lagrangian takes the equality rows (for an NMPC, the initial-state and continuity rows), while the inequality rows and the variable bounds go to the inner MMA or CCSAQ solver, which handles them directly. `argmin_auglag` is not the same method. With its default inner solver, `argmin_lbfgsb`, it folds the inequality rows into the augmented objective as well and leaves only the variable bounds to the inner solver (see [argmin_solver](argmin-solver.md)). The two names do not select interchangeable solvers.

### The adapter checks no dimension itself

`setup` reads `c_lower` and `c_upper` for `n_constraints` rows without checking their lengths, and passes `x_lower` and `x_upper` to NLopt when they are non-empty, which refuses a length other than `n_vars` (reported as `backend_rejected`). `solve` passes `update.x0` to NLopt, which refuses a length other than `n_vars` (reported as `solve_status::error`). A `solve` before any `setup`, or after a refused one, reports `solve_status::error`; a refused `setup` leaves no earlier problem posed.

### The constraints are differentiated by central differences

The adapter never calls `constraint_jacobian_cb`. Each of the three constraint groups (the equalities, the finite upper sides, the finite lower sides) forms its gradient by central differences over the whole constraint vector, with the step `cbrt(eps) * max(1, |x_j|)`, so one gradient request costs `2 * n_vars` evaluations of `constraints_cb` per group.

### Only SLSQP solved the shipped NMPC transcription at the default settings

The NMPC transcription poses the initial-state and continuity rows as equalities and leaves the initial state unbounded. Measured on the pendulum of the [nmpc](nmpc.md) usage example (horizon 20, input bound 5, 60 steps from 1 rad), with every other setting at its default: `slsqp` converged on all 60 solves with row violations of at most `1e-8` and brought the pendulum to 0.145 rad; `mma` and `ccsaq` are refused at setup with `incompatible_equality_constraints`; `isres` fails the first solve with `solver_error::invalid_problem`, because NLopt's ISRES refuses an infinite bound; `cobyla`, `auglag_mma` and `auglag_ccsaq` ended every solve on `max_eval` (`budget_exhausted`), with row violations between 0.08 and 4.5, and the pendulum ended between 1.65 and 2.87 rad instead of near zero.

## See Also

- [nmpc](nmpc.md)<br/> nonlinear MPC
- [nmhe](nmhe.md)<br/> nonlinear MHE
- [argmin-solver](argmin-solver.md)<br/> header-only argmin NLP solver backend
- [osqp-solver](osqp-solver.md)<br/> OSQP QP solver backend
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
- [background/nmpc](../../background/nmpc.md)<br/> NLP formulation and theory
