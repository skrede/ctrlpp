# argmin_qp_solver

> **Requires argmin 0.3.5 or a later 0.3 release** (`CTRLPP_ARGMIN_MIN_VERSION`). Enable with `-DCTRLPP_BUILD_ARGMIN=ON` when configuring CMake. Every argmin at or above that version ships `argmin/qp/`, so the backend is active whenever argmin is enabled.

argmin-native QP solver backend for linear MPC and MHE, and the header-only alternative to [`osqp_solver`](osqp-solver.md). Satisfies the `qp_solver` concept and can be named as the `Solver` template parameter of `mpc` and `mhe`. A preconfigured instance is injected through `mpc`'s `create` (its `solver` argument); `mhe` constructs its solver from the `Solver` type. Wraps `argmin::sparse_admm_qp_solver`, a header-only C++ implementation of the same OSQP-class operator-splitting algorithm, returning results on a typed error channel rather than through C-pointer ownership.

Both backends solve the canonical OSQP-form QP `min ½·xᵀPx + qᵀx s.t. l ≤ Ax ≤ u` and are interchangeable at a call site. Choosing this backend drops the vendored OSQP C library (and its C-pointer-ownership failure modes) from the build in favor of a header-only dependency.

> **Not a real-time upgrade.** argmin's *sparse* QP variant is a host-tier solver that disclaims real-time safety in its own header (unbounded per-call heap). This backend is dependency relief plus a typed-error-channel upgrade, not an RT-safety upgrade; it is soft-RT exactly as `osqp_solver` is. See [rt-safety-matrix](../../rt-safety-matrix.md).

## Header

| Form | Header |
|------|--------|
| `argmin_qp_solver` | `#include <ctrlpp/mpc/argmin_qp_solver.h>` |

```cpp
class argmin_qp_solver;
```

The class is defined only when `CTRLPP_BUILD_ARGMIN=ON` **and** argmin's QP header (`argmin/qp/sparse_admm_qp.h`) is present. When present, the macro `CTRLPP_HAS_ARGMIN_QP` is defined. Every argmin at or above `CTRLPP_ARGMIN_MIN_VERSION` (0.3.5) ships the header, so no local checkout is required. An argmin below that version is not supported, and the header gate is not a fallback for one; see [Solver backend components](../../getting-started.md#solver-backend-components) for how an older installed argmin is refused.

## qp_solver Concept

`argmin_qp_solver` satisfies the `qp_solver` concept defined in `<ctrlpp/mpc/qp_solver.h>`. The concept accepts one setup shape, a fallible `setup(problem)` returning a `ctrlpp::expected<void, E>` over the backend's own setup-error enum, so this solver models the concept in every build mode, including `-fno-exceptions` and `CTRLPP_NO_EXCEPTIONS`. Any type satisfying the concept can replace it as the solver backend.

## Constructors

```cpp
explicit argmin_qp_solver(double eps_abs = 1e-3,
                          double eps_rel = 1e-3,
                          std::int32_t max_iter = 4000,
                          bool verbose = false,
                          bool warm_starting = true,
                          bool polishing = true);
```

Constructs the solver, binding argmin's stable QP contract: tolerances, iteration budget, warm start and polish. The signature mirrors `osqp_solver` so the two backends are drop-in interchangeable.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eps_abs` | `1e-3` | Absolute convergence tolerance. Must be finite and nonnegative, and not zero together with `eps_rel`, or `setup` refuses it with `argmin_qp_setup_error::invalid_options` |
| `eps_rel` | `1e-3` | Relative convergence tolerance. Must be finite and nonnegative, and not zero together with `eps_abs`, or `setup` refuses it with `argmin_qp_setup_error::invalid_options` |
| `max_iter` | `4000` | Maximum ADMM iterations. Valid range `[1, 2147483647]`; `setup` refuses 0 or a negative value with `argmin_qp_setup_error::invalid_options` |
| `verbose` | `false` | Accepted for signature parity; argmin's sparse solver has no progress print, so this is ignored |
| `warm_starting` | `true` | Enable primal/dual warm-starting |
| `polishing` | `true` | Enable solution polishing |

The operator-splitting knobs (rho / sigma / alpha / adaptive_rho) are deliberately not exposed: per argmin coordination they are the volatile surface argmin intends to demote to an opt-in sub-struct, so this policy never touches them and a future reshape of those knobs is a no-op here.

### Preset constructor

```cpp
explicit argmin_qp_solver(qp_preset preset);
```

Constructs the solver from the backend-agnostic accuracy/speed preset instead of individual settings. The presets differ in a single knob, solution polishing:

| Preset | Polishing | When to use |
|--------|-----------|-------------|
| `qp_preset::accuracy` | on | The QP solution feeds a tolerance-sensitive consumer, or tight constraint feasibility is required. |
| `qp_preset::speed` | off | Warm-resolve linear MPC, where the unpolished iterate already meets the control tolerance and the per-step polish refinement is unnecessary. |

`qp_preset` (defined in `ctrlpp/mpc/qp_types.h`) is shared with `osqp_solver` and every other `qp_solver`, so a call site can pick the accuracy/speed tradeoff independently of the chosen backend:

```cpp
ctrlpp::argmin_qp_solver solver{ctrlpp::qp_preset::speed};
```

## Supporting Types

### argmin_qp_setup_error

Defined in `<ctrlpp/mpc/argmin_qp_solver.h>`.

```cpp
enum class argmin_qp_setup_error : std::uint8_t {
    pose_failed,     // argmin rejected the problem data (dimension mismatch, non-finite
                     // data, invalid bounds, an asymmetric P, or an unposeable KKT system)
    invalid_options  // an option lies outside argmin's domain (a budget outside
                     // [1, 2147483647], or tolerances negative, non-finite or both zero)
};
```

## Methods

### setup

```cpp
auto setup(const qp_problem<double>& problem)
    -> ctrlpp::expected<void, argmin_qp_setup_error>;
```

Fallible setup: poses and factorizes the problem once from the cost matrices P, q and constraint matrices A, l, u. P is read as upper-triangular, the same convention `osqp_solver` uses and `build_cost_matrix` emits; setup expands it to the full symmetric matrix argmin's solver requires. Returns an empty `expected` on success and an `argmin_qp_setup_error` on failure: `invalid_options` when an option lies outside argmin's domain (the budget or tolerances above), checked before the problem is posed, and `pose_failed` when argmin rejects the problem data. Either way no factorization exists, so `solve` returns `solve_status::error` until a `setup` succeeds. This is the only setup API and is available in every build, including `-fno-exceptions` and `CTRLPP_NO_EXCEPTIONS`.

### solve

```cpp
auto solve(const qp_update<double>& update) -> qp_result<double>;
```

Vectors-only resolve reusing the frozen factorization: warm-starts from the retained iterate (or from `update.warm_x`/`update.warm_y` when supplied) and runs argmin's `resolve_into`. Returns a `qp_result` with the primal/dual solution, solver status, objective value, iteration count, and residuals. Calling `solve` before a successful `setup` returns `solve_status::error`. Note `qp_result::solve_time` is always `0`: argmin's solver takes no timing measurements.

## Type Aliases

```cpp
using scalar_type = double;
```

The sparse ADMM QP backend operates in double precision only.

## Usage Example

`argmin_qp_solver` is a drop-in replacement for `osqp_solver`: swap the `Solver` template argument. A preconfigured instance is injected through the `solver` argument of `mpc`'s `create`, as below; `mhe` takes no solver argument and constructs its solver from the `Solver` type.

```cpp
#include <ctrlpp/mpc.h>
#include <ctrlpp/mpc/argmin_qp_solver.h>
#include <ctrlpp/model/propagate.h>

#include <Eigen/Dense>

#include <cstdint>
#include <cstdlib>
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

    // Skip per-step polishing on the warm-resolve path; see qp_preset.
    auto created = ctrlpp::mpc<double, NX, NU, ctrlpp::argmin_qp_solver>::create(
        sys, cfg, ctrlpp::argmin_qp_solver{ctrlpp::qp_preset::speed});
    if(!created)
    {
        // created.error() carries the controller_construction_error.
        return EXIT_FAILURE;
    }
    auto controller = *std::move(created);

    Eigen::Vector2d x(3.0, 0.0);
    for(std::int32_t step = 0; step < 20; ++step)
    {
        auto u_opt = controller.solve(x);
        if(!u_opt)
        {
            // u_opt.error() carries the solver_error describing the failure.
            std::cerr << "MPC solve failed at step " << step << "\n";
            return EXIT_FAILURE;
        }
        x = ctrlpp::propagate(sys, x, u_opt->input);
    }
    std::cout << "x after 20 steps: [" << x.transpose() << "]\n";
}
```

## See Also

- [osqp_solver](osqp-solver.md)<br/> OSQP QP backend (vendored C library) with the same concept surface
- [mpc](mpc.md)<br/> linear MPC
- [mhe](mhe.md)<br/> linear MHE
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
- [background/mpc](../../background/mpc.md)<br/> QP formulation and theory
