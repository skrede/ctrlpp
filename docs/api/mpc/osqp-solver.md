# osqp_solver

> **Requires OSQP.** Enable with `-DCTRLPP_BUILD_OSQP=ON` when configuring CMake.

OSQP-based QP solver backend for linear MPC and MHE. Satisfies the `qp_solver` concept and can be injected as the `Solver` template parameter of `mpc` and `mhe`. Wraps the OSQP C library with RAII lifetime management, automatic warm-starting, and configurable tolerances.

## Header

| Form | Header |
|------|--------|
| `osqp_solver` | `#include <ctrlpp/mpc/osqp_solver.h>` |

```cpp
class osqp_solver;
```

## qp_solver Concept

`osqp_solver` satisfies the `qp_solver` concept defined in `<ctrlpp/mpc/qp_solver.h>`:

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

Any type satisfying this concept can replace `osqp_solver` as the solver backend. There is one setup shape: `setup(problem)` returns a `ctrlpp::expected<void, E>` over the backend's own setup-error enum, through which the backend reports its setup failures. A backend with no setup failure mode writes a trivially succeeding fallible setup rather than an infallible one. `argmin_solver`'s setup refuses the problems its bridge cannot hold, with the `argmin_setup_error` enumerators listed under [argmin setup refusals](nmpc.md#argmin-setup-refusals). `osqp_solver` provides `setup` in every build mode, including `-fno-exceptions` and `CTRLPP_NO_EXCEPTIONS`.

## Constructors

```cpp
explicit osqp_solver(double eps_abs = 1e-3,
                     double eps_rel = 1e-3,
                     std::int32_t max_iter = 4000,
                     bool verbose = false,
                     bool warm_starting = true,
                     bool polishing = true);
```

Constructs the solver with OSQP settings. Default values provide a good balance of speed and accuracy for typical MPC/MHE problems.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eps_abs` | `1e-3` | Absolute convergence tolerance |
| `eps_rel` | `1e-3` | Relative convergence tolerance |
| `max_iter` | `4000` | Maximum ADMM iterations. Valid range `[1, 2147483647]`; `setup` refuses 0 or a negative value with `osqp_setup_error::setup_failed` |
| `verbose` | `false` | Print solver progress |
| `warm_starting` | `true` | Enable primal/dual warm-starting |
| `polishing` | `true` | Enable solution polishing |

### Preset constructor

```cpp
explicit osqp_solver(qp_preset preset);
```

Constructs the solver from a backend-agnostic accuracy/speed preset instead of individual settings. The presets differ in a single knob, solution polishing:

| Preset | Polishing | When to use |
|--------|-----------|-------------|
| `qp_preset::accuracy` | on | The QP solution feeds a tolerance-sensitive consumer, or tight constraint feasibility is required. |
| `qp_preset::speed` | off | Warm-resolve linear MPC, where the unpolished iterate already meets the control tolerance and the per-step polish refinement is unnecessary. |

Without polishing, OSQP stops once the primal residual is within `eps_abs + eps_rel * max(‖Ax‖∞, ‖z‖∞)`, so a constraint can be exceeded by up to that amount. The regulation example's first input is -1.0016 against `|u| <= 1`. Polishing re-solves on the active set and normally lands on the bound, but when polishing fails OSQP keeps the unpolished iterate. Use `qp_preset::accuracy` together with a tighter `eps_abs` and `eps_rel` when a bound must hold to a tighter tolerance.

`qp_preset` (defined in `ctrlpp/mpc/qp_types.h`) is shared across every backend modeling the `qp_solver` concept, so a call site can pick the accuracy/speed tradeoff independently of the chosen solver:

```cpp
ctrlpp::osqp_solver solver{ctrlpp::qp_preset::speed};
```

## Supporting Types

### osqp_setup_error

Defined in `<ctrlpp/mpc/qp_types.h>` beside `solve_status` (which is unchanged; the per-iteration solve result stays a plain value enum).

```cpp
enum class osqp_setup_error : std::uint8_t {
    settings_allocation_failed,  // OSQPSettings_new returned null
    setup_failed                 // osqp_setup rejected the data or settings
};
```

## Methods

### setup

```cpp
auto setup(const qp_problem<double>& problem)
    -> ctrlpp::expected<void, osqp_setup_error>;
```

Fallible setup: initializes the OSQP workspace from a QP problem (cost matrices P, q and constraint matrices A, l, u). Returns an empty `expected` on success and an `osqp_setup_error` on failure. This is the only setup API and is available in every build, including `-fno-exceptions` and `CTRLPP_NO_EXCEPTIONS` builds.

An absent bound side is passed to OSQP in its own convention, `OSQP_INFTY` (1e30 in a double build). A lower bound at or below `-OSQP_INFTY`, including `-inf`, is passed as `-OSQP_INFTY`, and an upper bound at or above `+OSQP_INFTY`, including `+inf`, as `+OSQP_INFTY`. OSQP does not map an IEEE infinity itself: its duality gap would evaluate `inf * 0` to NaN, and every solve would run to `max_iter`. Every other value, NaN included, is passed unchanged. A lower bound of `+inf` or an upper bound of `-inf` poses an empty row rather than an absent side, so it is passed unchanged and OSQP's `l <= u` check refuses it with `osqp_setup_error::setup_failed`. A problem whose `l` or `u` length differs from the row count of A is refused the same way.

OSQP reads `OSQP_INFTY` as absent through its Ruiz equilibration, so on a row whose entries are far above 1e4 an absent side can read as a very large finite bound instead.

### solve

```cpp
auto solve(const qp_update<double>& update) -> qp_result<double>;
```

Solves the QP with updated cost vector, bounds, and optional warm-start vectors. Returns a `qp_result` containing the primal/dual solution, solver status, objective value, iteration count, solve time, and residuals. An empty vector in the update leaves that part of the problem as it was.

The bounds are mapped to `OSQP_INFTY` as at setup, into buffers sized at setup, so a solve allocates nothing for them. The result is `solve_status::error`, with no solve run, when:

- `l` or `u` has a length other than 0 or the posed row count;
- OSQP refuses the update, for example an empty row from an infinite initial state. OSQP keeps the previous q, l and u on a refused update, so solving anyway would answer the previous problem;
- no setup has succeeded yet.

## Type Aliases

```cpp
using scalar_type = double;
```

OSQP operates in double precision only.

## Usage Example

```cpp
// gnuplot: plot "< ./osqp_demo" using 1:3 with lines title "state"
#include <ctrlpp/mpc.h>
#include <ctrlpp/mpc/osqp_solver.h>
#include <ctrlpp/model/propagate.h>

#include <Eigen/Dense>

#include <utility>
#include <iostream>

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
        .horizon = 15,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity()};

    // Custom solver settings: tighter tolerance, fewer iterations
    ctrlpp::osqp_solver solver(1e-5, 1e-5, 2000);

    auto created = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(sys, cfg, std::move(solver));
    if(!created)
        return 1;
    auto& controller = *created;

    Eigen::Vector2d x(3.0, 0.0);

    for(int k = 0; k < 50; ++k)
    {
        auto u_opt = controller.solve(x);
        if(!u_opt)
            break;

        x = ctrlpp::propagate(sys, x, u_opt->input);
        std::cout << "k=" << k << "  x=[" << x.transpose() << "]\n";
    }
}
```

## See Also

- [mpc](mpc.md)<br/> linear MPC using OSQP
- [mhe](mhe.md)<br/> linear MHE using OSQP
- [argmin-qp-solver](argmin-qp-solver.md)<br/> header-only argmin QP backend with the same concept surface
- [nlopt-solver](nlopt-solver.md)<br/> NLopt NLP solver backend
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
- [background/mpc](../../background/mpc.md)<br/> QP formulation and theory
