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
    && requires(S solver, const nlp_update<typename S::scalar_type>& upd) {
        { solver.solve(upd) } -> std::same_as<nlp_result<typename S::scalar_type>>;
    }
    && requires(S solver, const nlp_problem<typename S::scalar_type>& prob) {
        { solver.setup(prob).has_value() } -> std::convertible_to<bool>;
    };
```

There is one setup shape: `setup(problem)` returns a `ctrlpp::expected<void, E>` over the backend's own setup-error enum, through which the backend reports its setup failures. A backend with no setup failure mode writes a trivially succeeding fallible setup rather than an infallible one. `argmin_solver` built with `Constrained = false` handles variable bounds only, so its setup refuses a problem with general constraint rows with `argmin_setup_error::constraints_unsupported`.

> **Note:** NLopt's C++ API throws by upstream design. `setup` and `solve` contain every exception that NLopt, the adapter or a problem callable raises, including one of a type not derived from `std::exception`, so the containment lives in this adapter rather than in the controllers and estimators that consume it. In `setup`, a `std::bad_alloc` reports `nlopt_setup_error::allocation_failed` and anything else reports `nlopt_setup_error::backend_rejected`. In `solve`, any exception reports `solve_status::error`. The adapter itself still requires exception support and is exempt from the library's embedded no-exceptions floor.
>
> Because every exception is contained, the backend does not support thread cancellation during a `setup` or `solve`. Where the platform implements cancellation as a forced unwind, as glibc does, NLopt's own callback wrappers and the adapter's handlers catch that unwind and do not re-raise it, and glibc then aborts the process.

## Supporting Types

### nlopt_algorithm

```cpp
enum class nlopt_algorithm : std::uint8_t {
    slsqp,         // Sequential Least Squares Programming (gradient-based, supports equality)
    mma,           // Method of Moving Asymptotes (gradient-based, NO equality constraints)
    cobyla,        // Constrained Optimization BY Linear Approximations (derivative-free)
    isres,         // Improved Stochastic Ranking Evolution Strategy (global, derivative-free)
    auglag_mma,    // AUGLAG_EQ outer with LD_MMA inner (absorbs equality constraints)
    auglag_ccsaq,  // AUGLAG_EQ outer with LD_CCSAQ inner (absorbs equality constraints)
    ccsaq          // Conservative Convex Separable Approximation (NO equality constraints)
};
```

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

The solver is copyable and movable. NLopt registers every callback with the solver's address, so a copy or a move re-runs `setup` against the same referenced problem to rebind the callbacks to the new object. A constructor or assignment cannot return an error, so when that re-setup fails (because the referenced problem changed since the source was set up) the new solver records it and its next `solve` reports `solve_status::error`. The same holds for a copy of a solver that was never set up. Every setup, including this re-setup, rebuilds the constraint partition from scratch, so a solver re-set up against a problem without constraint rows registers no constraint callbacks and evaluates no constraint, whatever problem it was set up on before.

## Methods

### setup

```cpp
auto setup(const nlp_problem<Scalar>& problem)
    -> ctrlpp::expected<void, nlopt_setup_error>;
```

Fallible setup: configures the NLopt optimizer from an NLP problem definition and partitions constraints into equality and inequality groups automatically. Returns an empty `expected` on success, `nlopt_setup_error::incompatible_equality_constraints` if raw MMA or raw CCSAQ is selected on a problem with equality constraints, `nlopt_setup_error::allocation_failed` if an allocation fails while the optimizer is built or configured, and `nlopt_setup_error::backend_rejected` for any other error raised meanwhile, standard or not. After a failed setup the solver holds no usable optimizer, and `solve` reports `solve_status::error` until a later `setup` succeeds. This is the only setup API.

### solve

```cpp
auto solve(const nlp_update<Scalar>& update) -> nlp_result<Scalar>;
```

Solves the NLP from the initial guess in `update.x0`. Returns the solution, objective value, solver status, iteration count, solve time, and maximum constraint violation.

NLopt requires the starting point to lie inside the variable bounds, so `solve` first projects each component of `update.x0` onto `[x_lower[i], x_upper[i]]`. A consistent box that excludes the previous iterate (for example an input box that excludes the controllers' initial zero warm start) is therefore solved rather than refused. Infinite bounds leave a component unchanged, and inverted or NaN bounds are left for NLopt to refuse, which `solve` reports as an error.

A round-off-limited stop reports `solve_status::solved_inaccurate` with the objective re-evaluated at the returned point. Every other exception raised by NLopt or by a problem callable during the solve, of any type and including one raised by that re-evaluation or by the constraint evaluation behind `primal_residual`, and including a start vector whose length differs from `n_vars`, reports `solve_status::error` with an empty `x` and zero objective, iteration count and violation; the controllers map that status to `solver_error::invalid_problem`.

## Usage Example

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

## See Also

- [nmpc](nmpc.md)<br/> nonlinear MPC using NLopt
- [nmhe](nmhe.md)<br/> nonlinear MHE using NLopt
- [osqp-solver](osqp-solver.md)<br/> OSQP QP solver backend
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
- [background/mpc](../../background/mpc.md)<br/> NLP formulation and theory
