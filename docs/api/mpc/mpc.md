# mpc

> **Requires OSQP.** Enable with `-DCTRLPP_BUILD_OSQP=ON` when configuring CMake.

Linear Model Predictive Controller using sparse QP optimization. Solves a receding-horizon regulation or tracking problem for a discrete-time linear state-space system with input, state, and rate constraints. Terminal cost defaults to the DARE solution when not provided explicitly. State constraints are softened by default to ensure solver feasibility.

## Header and Alias

| Form | Header |
|------|--------|
| `mpc<Scalar, NX, NU, Solver>` | `#include <ctrlpp/mpc.h>` |

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, qp_solver Solver>
class mpc;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type (`double`, `float`) |
| `NX` | `>= 1` | State dimension |
| `NU` | `>= 1` | Input dimension |
| `Solver` | satisfies `qp_solver` | QP solver backend. See [osqp-solver](osqp-solver.md). |

## mpc_config

Configuration struct `mpc_config<Scalar, NX, NU>` passed at construction.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `horizon` | `int` | `1` | Prediction horizon length N |
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | State cost weight |
| `R` | `Matrix<Scalar, NU, NU>` | Identity | Input cost weight |
| `Qf` | `optional<Matrix<Scalar, NX, NX>>` | DARE solution | Terminal cost weight. Computed from DARE when not provided; if that solve fails the state weight stands in and [`diagnostics()`](#diagnostics) reports it. |
| `u_min` | `optional<Vector<Scalar, NU>>` | none | Element-wise lower input bound |
| `u_max` | `optional<Vector<Scalar, NU>>` | none | Element-wise upper input bound |
| `x_min` | `optional<Vector<Scalar, NX>>` | none | Element-wise lower state bound |
| `x_max` | `optional<Vector<Scalar, NX>>` | none | Element-wise upper state bound |
| `du_max` | `optional<Vector<Scalar, NU>>` | none | Maximum input rate of change per step |
| `soft_penalty` | `Scalar` | `1e4` | L1 penalty for soft state constraints |
| `soft_state_penalty` | `optional<Vector<Scalar, NX>>` | none | Per-state soft constraint penalty |
| `terminal_constraint_set` | `optional<terminal_set<Scalar, NX>>` | none | Ellipsoidal or polyhedral terminal set; see [Terminal sets](#terminal-sets) for which one is invariant |
| `hard_state_constraints` | `bool` | `false` | When true, state constraints are hard (not softened) |

## Construction

```cpp
static auto create(const discrete_state_space<Scalar, NX, NU, NX>& system,
                       const mpc_config<Scalar, NX, NU>& config)
    -> expected<mpc, controller_construction_error>;

static auto create(const discrete_state_space<Scalar, NX, NU, NX>& system,
                       const mpc_config<Scalar, NX, NU>& config,
                       Solver solver)
    -> expected<mpc, controller_construction_error>;
```

`create` is the only construction path. There is no non-fallible constructor, so a rejected configuration is a value the caller has to inspect and never a controller that quietly stands in for one. It validates the configuration, then builds the QP matrices, computes the terminal cost (via DARE if `Qf` is not set), and initializes the solver. The two-argument form default-constructs the solver and chains into the three-argument form, which injects a caller-supplied, pre-configured solver (moved in before the initial QP is posed), letting you choose the accuracy/speed tradeoff via a preset:

```cpp
// Skip per-step polishing on the warm-resolve path -- see qp_preset.
auto created = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(
    sys, cfg, ctrlpp::osqp_solver{ctrlpp::qp_preset::speed});
if (!created)
    return handle(created.error());
auto controller = *std::move(created);
```

### Rejections

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
| `soft_penalty`, or a set `soft_state_penalty`, is non-finite or non-positive in some entry while the state bounds are soft | `invalid_constraint_penalty` |

Each weight is checked completely (finite, symmetric, semidefinite) before the next one, so an indefinite `Q` is reported ahead of a non-finite `R`.

Equal bounds pin a variable and are accepted, as is an infinite bound in its open direction (`u_min = -inf`, `u_max = +inf`). A rate bound of zero, `-0.0` or `+inf` is accepted.

The weights follow the library's symmetry convention. An asymmetry within the rounding bound `N * eps * max|W_ij|` (never below `N * denorm_min`) of an `N x N` weight `W` is drift: it is symmetrized in the controller's own copy of the configuration before the QP is assembled. The cost assembly passes only the upper triangle to the backend, so without that step the backend would see the reflection of the upper triangle rather than the symmetric part the configured cost denotes. A larger asymmetry is structure and is refused, never silently symmetrized. Semidefiniteness is decided by the symmetric eigensolver with the floor `N * eps * max|W_ij|`, the order of its backward error. `R` is required positive semidefinite, not definite: the QP needs only a semidefinite Hessian to be convex, and the one consumer that needs `R` nonsingular, the Riccati solve for a missing `Qf`, already reports its failure through [`diagnostics()`](#diagnostics).

The state bounds are soft when `x_min` or `x_max` is set and `hard_state_constraints` is false; only then are the penalties inspected. `soft_penalty` is checked even while a per-state `soft_state_penalty` replaces it, because it is a configuration field with a valid default and an invalid value there is a mistake either way.

`horizon` stays a signed `int` deliberately: a mistaken negative value remains representable as negative and is therefore rejectable, whereas an unsigned field would turn the same mistake into an enormous allocation. The overflow bound is a representability condition on `int`, not a chosen ceiling. At the worst-case configuration the horizon `N` scales the decision vector as `N*(2*NX + NU) + NX` and the constraint rows as `N*(2*NX + 2*NU) + NX + n_terminal`, so both stay representable exactly when `horizon <= (INT_MAX - (NX + n_terminal)) / (2*NX + 2*NU)`.

Two alternatives are deliberately not implemented. Validating at the first solve would surface a configuration error at the first control step, the worst possible moment. Clamping the horizon to one would turn a caller mistake into a silently different controller.

## Failure contract

Both `mpc` and `nmpc` share one soft-constraint and failure contract. Every
`solve` overload returns `ctrlpp::expected<solve_output<Scalar, NU>, solver_error>`:

- On the **success branch** the result holds a `solve_output` whose `input` field
  is the first control input to apply and whose `status` field is a soft
  `solve_result_status`:
  - `converged`: the solver met its convergence tolerances.
  - `solved_inaccurate`: a usable iterate was returned but tolerances were only
    partially met.
  - `budget_exhausted`: the iteration or time budget ran out, but the best iterate
    found is still returned so it can be commanded knowingly.
- On the **error branch** the result holds a `solver_error` and no input:
  - `infeasible`: the problem as posed has no feasible point.
  - `invalid_problem`: the problem is unbounded, non-convex, the solver reported an
    internal error, or a reference span was too short.
  - `setup_incomplete`: the one-time solver setup failed, so no solve can run.
  - `invalid_backend_result`: the solver reported a status the controller accepts
    and then returned a primal shorter than the decision dimension, or a dual
    shorter than the constraint count, of the problem posed at construction. The
    controller checks both reported lengths before it consumes either vector, so
    the fixed-width slices it takes out of the primal cannot read past the end of
    the solver's own storage, and an undersized dual is never handed back as the
    next warm start. This is deliberately distinct from `invalid_problem`: there
    the problem must be fixed, here the problem is well formed and the backend's
    answer is not. A longer-than-required result is accepted; only the condition
    that makes the reads legal is enforced.

The input is reached explicitly through `->input`; there is no implicit conversion
to `Vector<Scalar, NU>`, so the soft status can never be silently dropped. On the
error branch the controller does not update its internal previous-input record and
applies no hidden fallback, so a failed solve never warms the rate constraint from
a phantom input. Use `set_applied_input` to record the input the caller actually
commanded.

## Methods

### solve (regulation)

```cpp
auto solve(const Vector<Scalar, NX>& x0)
    -> ctrlpp::expected<solve_output<Scalar, NU>, solver_error>;
```

Solves the QP for regulating state to the origin. Returns the success branch with
the first input and a soft status, or the error branch on failure. If solver setup
failed at construction (reported through the solver's `setup`), every `solve`
overload returns the error branch with `solver_error::setup_incomplete`.

### solve (constant reference)

```cpp
auto solve(const Vector<Scalar, NX>& x0,
                         const Vector<Scalar, NY>& y_ref)
    -> ctrlpp::expected<solve_output<Scalar, NU>, solver_error>;
```

Solves the QP for tracking a constant output reference across the entire horizon.

### solve (trajectory reference)

```cpp
auto solve(const Vector<Scalar, NX>& x0,
                         std::span<const Vector<Scalar, NY>> y_ref)
    -> ctrlpp::expected<solve_output<Scalar, NU>, solver_error>;
```

Solves the QP for tracking a time-varying output reference trajectory. The span
must contain at least `horizon + 1` elements; an undersized span is rejected via
the error branch (`solver_error::invalid_problem`) rather than read out of bounds.

### set_applied_input

```cpp
void set_applied_input(const Vector<Scalar, NU>& u);
```

Records the control input the caller actually commanded. The internal previous
input (which anchors the rate constraint) is updated only on a successful solve;
this accessor lets the caller keep that reference consistent when it commands a
different input, for example after an error branch or external saturation.

### trajectory

```cpp
auto trajectory() const
    -> ctrlpp::expected<std::pair<std::vector<Vector<Scalar, NX>>,
                                  std::vector<Vector<Scalar, NU>>>,
                        solver_error>;
```

Returns the full predicted state and input trajectories from the last solve. It is
guarded: before the first valid solve it returns the error branch
(`solver_error::setup_incomplete`) rather than stale or default-initialized data.

### diagnostics

```cpp
auto diagnostics() const -> mpc_diagnostics<Scalar>;
```

Returns solver diagnostics from the last solve, including status, iteration count, solve time, cost, and residuals.

It also carries one disposition that is not about a solve at all:

| Field | Meaning |
|-------|---------|
| `used_state_weight_terminal_cost` | No terminal weight was configured, and the Riccati solve that would have supplied one produced no solution, so `Q` (mapped to state space) stands in for the terminal cost |

When `Qf` is not set, the infinite-horizon Riccati solution is what the terminal cost is meant to be. That solve can fail, for instance on a rank-deficient `A`. The controller then substitutes the state weight, which is a usable cost but a different problem from the configured one: every solve after it optimizes that different problem, and a stability argument resting on the terminal cost no longer holds. Nothing about a solve's success distinguishes the two, so the substitution is reported here instead.

The flag is latched when the terminal cost is computed, which happens at construction rather than per solve. It is therefore readable before the first solve, and reads the same on every solve afterwards.

```cpp
auto controller = *ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(system, config);

if (controller.diagnostics().used_state_weight_terminal_cost)
{
    // Not the terminal cost the configuration implies. Either supply Qf
    // explicitly, or accept that the horizon is doing all the work.
}
```

## Terminal sets

```cpp
#include <ctrlpp/mpc/invariant_set.h>

template <typename Scalar, std::size_t NX, std::size_t NU>
auto terminal_ingredients(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B,
                          const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R,
                          const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max,
                          const Vector<Scalar, NX>& x_min, const Vector<Scalar, NX>& x_max)
    -> expected<terminal_ingredients_result<Scalar, NX, NU>, terminal_set_error>;
```

`terminal_ingredients` solves the DARE for `P`, forms the LQR gain `K = -(R + B^T P B)^{-1} B^T P A`, and returns three fields. An infinite bound side is absent, so pass `+/-infinity` to leave a state or input unbounded on that side.

| Field | Type | What it is |
|-------|------|------------|
| `Qf` | `Matrix<Scalar, NX, NX>` | The DARE solution `P`, the terminal cost |
| `set` | `ellipsoidal_set<Scalar, NX>` | The ellipsoid `{x : x^T P x <= alpha}` with the largest `alpha` inside the input and state faces. The QP encodes its inscribed box in `P`'s eigenbasis, which is not positively invariant under `u = K x` and gives no recursive-feasibility guarantee |
| `invariant_set` | `expected<polytopic_set<Scalar, NX>, terminal_set_error>` | The maximal output admissible set of `A + B K` inside `{x_min <= x <= x_max, u_min <= K x <= u_max}`: the invariant terminal set of the LQR law, as a minimal set of unit rows |

The function's own error channel carries `dare_failed`, `input_zero_not_interior` and `empty_terminal_set`, the failures of `Qf` and `set`. The output admissible set has its own outcome in its field, so its refusal never fails the call:

| `invariant_set` error | Cause |
|-----------------------|-------|
| `origin_not_interior` | A finite state or input bound side does not strictly straddle zero, or a bound is NaN |
| `closed_loop_not_stable` | `A + B K` is not strictly stable (the Stein equation `A_K^T X A_K - X = -I` has no positive definite solution), or it or `K` is not finite |
| `unbounded_terminal_set` | The finitely bounded outputs do not observe `A + B K`, so the set is unbounded |
| `not_converged` | No fixed point within the step bound derived from the Lyapunov contraction of `A + B K` |
| `dimension_unsupported` | `NX > 4` |

### compute_output_admissible_set

```cpp
#include <ctrlpp/mpc/output_admissible_set.h>

template <typename Scalar, std::size_t NX, std::size_t NU>
auto compute_output_admissible_set(const Matrix<Scalar, NX, NX>& A_K, const Matrix<Scalar, NU, NX>& K,
                                   const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max,
                                   const Vector<Scalar, NX>& x_min, const Vector<Scalar, NX>& x_max,
                                   Scalar rounding_factor = Scalar{1})
    -> expected<polytopic_set<Scalar, NX>, terminal_set_error>;
```

Builds the maximal output admissible set `O_inf = {x : C A_K^k x in Y for all k >= 0}` of `x+ = A_K x`, with `C = [I; K]` and `Y` the finite bound sides, for any closed loop (the `invariant_set` field calls it with `A_K = A + B K`). Its hypotheses are checked in this order: every finite bound side strictly straddles zero (else `origin_not_interior`), `A_K` is strictly stable and `A_K` and `K` are finite (else `closed_loop_not_stable`), and the finitely bounded outputs observe `A_K` (else `unbounded_terminal_set`). `not_converged` is reported if the construction passes its derived step bound. `NX <= 4` is a compile-time requirement.

The construction stacks the rows `C A_K^k` until they bound the set, then appends, step by step, the rows that the current set's vertices do not imply, prunes redundant rows every step, and stops at the first step that appends nothing. There is no linear program: every implication and facet test runs on enumerated vertices, which costs `C(m, NX)` small solves per step for `m` rows. Each comparison of a row with a vertex carries its own rounding bound (the residual as an inner product, plus the vertex's forward error from the residual of its defining rows); `rounding_factor` scales the order-one prefactor of every such bound and must be positive. A row declared implied within its bound can leave the result larger than the exact set by at most that bound; a looser feasibility test only keeps more rows. On the double integrator with `dt = 0.1`, `Q = I`, `R = 1`, `|u| <= 1` and `|x_i| <= 10`, the set has 30 rows, the last needed at step 14.

### compute_polytopic_invariant_set

`compute_polytopic_invariant_set(A, B, state_constraints, input_constraints, max_iterations, convergence_tol)` iterates the robust pre-image to the set of states that stay inside the state polytope under every input in the input polytope. That is a robust control-invariant set, a valid object for robust MPC; it is not the terminal set of the LQR law. Use `invariant_set` for that.

### Using the invariant terminal set

```cpp
auto ingredients = ctrlpp::terminal_ingredients<double, NX, NU>(A, B, Q, R, u_min, u_max, x_min, x_max);
if(!ingredients || !ingredients->invariant_set)
    return EXIT_FAILURE;

ctrlpp::mpc_config<double, NX, NU> cfg{
    .horizon                 = 30,
    .Q                       = Q,
    .R                       = R,
    .Qf                      = ingredients->Qf,
    .u_min                   = u_min,
    .u_max                   = u_max,
    .x_min                   = x_min,
    .x_max                   = x_max,
    .terminal_constraint_set = ctrlpp::terminal_set<double, NX>{*ingredients->invariant_set},
    .hard_state_constraints  = true};
```

Recursive feasibility and stability of the origin follow from this configuration only with hard state constraints, a nominal model and regulation to the origin; see [background/mpc](../../background/mpc.md). The full closed loop is in `examples/mpc/ctrlpp_mpc_05_terminal_set.cpp`.

## Usage Example

```cpp
// gnuplot: plot "< ./mpc_regulation" using 1:3 with lines title "position"
#include <ctrlpp/mpc.h>
#include <ctrlpp/mpc/osqp_solver.h>
#include <ctrlpp/model/propagate.h>

#include <Eigen/Dense>

#include <cstdlib>
#include <iostream>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr double dt = 0.1;

    // Double integrator: x = [position, velocity]
    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        .A = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        .B = (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(),
        .C = Eigen::Matrix2d::Identity(),
        .D = Eigen::Matrix<double, 2, 1>::Zero()};

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 20,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-1.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(1.0)};

    auto created = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(sys, cfg);
    if(!created)
    {
        // created.error() carries the controller_construction_error.
        return 1;
    }
    auto controller = *std::move(created);

    Eigen::Vector2d x(5.0, 0.0);

    for(double t = 0.0; t < 5.0; t += dt)
    {
        auto u_opt = controller.solve(x);
        if(!u_opt)
        {
            // u_opt.error() carries the solver_error describing the failure.
            std::cerr << "MPC solve failed at t=" << t << "\n";
            return EXIT_FAILURE;
        }

        auto diag = controller.diagnostics();
        std::cout << "t=" << t << "  x=[" << x.transpose()
                  << "]  u=" << u_opt->input[0]
                  << "  cost=" << diag.cost << "\n";

        x = ctrlpp::propagate(sys, x, u_opt->input);
    }
}
```

## See Also

- [nmpc](nmpc.md)<br/> nonlinear MPC
- [osqp-solver](osqp-solver.md)<br/> OSQP QP solver backend
- [mhe](mhe.md)<br/> linear moving horizon estimation
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
- [guides/intro/your-first-mpc](../../guides/intro/your-first-mpc.md)<br/> introductory MPC tutorial
- [background/mpc](../../background/mpc.md)<br/> MPC theory and background
