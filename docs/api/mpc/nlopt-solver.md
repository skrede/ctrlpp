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
    && (requires(S solver, const nlp_problem<typename S::scalar_type>& prob) {
            { solver.setup(prob) } -> std::same_as<void>;
        }
        || requires(S solver, const nlp_problem<typename S::scalar_type>& prob) {
            { solver.try_setup(prob).has_value() } -> std::convertible_to<bool>;
        });
```

A solver models the concept with either setup shape: the classic `void setup(problem)` or the fallible `try_setup(problem)` returning an `expected<void, E>`. `nlopt_solver` provides `try_setup` unconditionally and keeps `setup` as a throwing convenience wrapper when exceptions are enabled.

> **Note:** NLopt's C++ API throws by upstream design (`solve` maps those exceptions to `solve_status` internally), so this opt-in backend requires exception support and is exempt from the library's embedded no-exceptions floor.

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
    incompatible_equality_constraints  // raw mma/ccsaq selected on a problem with equality constraints
};
```

Raw MMA and raw CCSAQ cannot handle equality constraints; the `auglag_mma` and `auglag_ccsaq` variants absorb them into the outer augmented-Lagrangian penalty and accept the same problem.

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

### try_setup

```cpp
[[nodiscard]] auto try_setup(const nlp_problem<Scalar>& problem)
    -> ctrlpp::expected<void, nlopt_setup_error>;
```

Fallible setup: configures the NLopt optimizer from an NLP problem definition and partitions constraints into equality and inequality groups automatically. Returns an empty `expected` on success and `nlopt_setup_error::incompatible_equality_constraints` if raw MMA or raw CCSAQ is selected on a problem with equality constraints.

### setup

```cpp
void setup(const nlp_problem<Scalar>& problem);  // only when CTRLPP_HAS_EXCEPTIONS
```

Throwing convenience wrapper over `try_setup`, available only when exceptions are enabled. Throws `std::invalid_argument` if MMA or CCSAQ is selected with equality constraints.

### solve

```cpp
auto solve(const nlp_update<Scalar>& update) -> nlp_result<Scalar>;
```

Solves the NLP from the initial guess in `update.x0`. Returns the solution, objective value, solver status, iteration count, solve time, and maximum constraint violation. The `solve_status` in the result is unchanged by the setup API split; the per-iteration solve result stays a plain value enum.

## Usage Example

```cpp
// gnuplot: plot "< ./nlopt_demo" using 1:3 with lines title "state"
#include <ctrlpp/nmpc.h>
#include <ctrlpp/mpc/nlopt_solver.h>

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

- [nmpc](nmpc.md)<br/> nonlinear MPC using NLopt
- [nmhe](nmhe.md)<br/> nonlinear MHE using NLopt
- [osqp-solver](osqp-solver.md)<br/> OSQP QP solver backend
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
- [background/mpc](../../background/mpc.md)<br/> NLP formulation and theory
