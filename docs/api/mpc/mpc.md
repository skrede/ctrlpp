# mpc

> **Requires a QP backend.** `osqp_solver` needs `-DCTRLPP_BUILD_OSQP=ON` and `argmin_qp_solver` needs `-DCTRLPP_BUILD_ARGMIN=ON` when configuring CMake.

Linear Model Predictive Controller using sparse QP optimization. Solves a receding-horizon regulation or tracking problem for a discrete-time linear state-space system with input, state, and rate constraints. Terminal cost defaults to the DARE solution when not provided explicitly. State bounds are softened by default with a quadratic penalty on a slack, so the state bounds alone cannot make a solve infeasible; the input, rate and terminal-set rows stay hard (see [Soft state constraints](#soft-state-constraints)).

## Header and Alias

| Form | Header |
|------|--------|
| `mpc<Scalar, NX, NU, Solver, NY>` | `#include <ctrlpp/mpc.h>` |
| `mpc_config<Scalar, NX, NU, NY>` | `#include <ctrlpp/mpc/mpc_config.h>` (also included by `ctrlpp/mpc.h`) |
| `mpc_diagnostics<Scalar>` | `#include <ctrlpp/mpc/diagnostics.h>` (also included by `ctrlpp/mpc.h`) |
| `solve_output`, `solver_error`, `controller_construction_error`, `qp_problem`, `qp_update`, `qp_result` | `#include <ctrlpp/mpc/qp_types.h>` (also included by `ctrlpp/mpc.h`) |
| the `qp_solver` concept | `#include <ctrlpp/mpc/qp_solver.h>` (also included by `ctrlpp/mpc.h`) |
| `mpc_construction_error<SetupError>` | `#include <ctrlpp/mpc/terminal_cost.h>` (also included by `ctrlpp/mpc.h`) |

`ctrlpp/mpc.h` also includes two headers whose functions live in `ctrlpp::detail` and are not called directly: `ctrlpp/mpc/qp_formulation.h` assembles the sparse QP described under [Soft state constraints](#soft-state-constraints), and `ctrlpp/mpc/config_validation.h` holds the admission rules listed under [Rejections](#rejections), which `nmpc` shares.

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU, qp_solver Solver, std::size_t NY = NX>
class mpc;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar`, and the solver's `scalar_type` | Numeric type. Both shipped backends are `double` only |
| `NX` | `>= 1` | State dimension |
| `NU` | `>= 1` | Input dimension |
| `Solver` | satisfies `qp_solver` | QP solver backend. See [osqp-solver](osqp-solver.md) and [argmin-qp-solver](argmin-qp-solver.md). |
| `NY` | defaults to `NX` | Output dimension: the weights `Q` and `Qf` and the tracking references are in output space |

## mpc_config

Configuration struct `mpc_config<Scalar, NX, NU, NY = NX>` passed at construction.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `horizon` | `std::int32_t` | `1` | Prediction horizon length N |
| `Q` | `Matrix<Scalar, NY, NY>` | Identity | Output cost weight; the controller weighs the state by `C^T Q C` |
| `R` | `Matrix<Scalar, NU, NU>` | Identity | Input cost weight, about `u = 0` |
| `Qf` | `optional<Matrix<Scalar, NY, NY>>` | DARE solution | Terminal output weight, posed on the state as `C^T Qf C`. When not provided, the terminal state weight is the DARE solution for `(A, B, C^T Q C, R)`; `create` refuses the configuration if that solve declines. |
| `u_min` | `optional<Vector<Scalar, NU>>` | none | Element-wise lower input bound on `u_0 ... u_{N-1}` |
| `u_max` | `optional<Vector<Scalar, NU>>` | none | Element-wise upper input bound on `u_0 ... u_{N-1}` |
| `x_min` | `optional<Vector<Scalar, NX>>` | none | Element-wise lower state bound on `x_1 ... x_N` |
| `x_max` | `optional<Vector<Scalar, NX>>` | none | Element-wise upper state bound on `x_1 ... x_N` |
| `du_max` | `optional<Vector<Scalar, NU>>` | none | Symmetric rate bound `|u_k - u_{k-1}| <= du_max`; for `k = 0`, `u_{-1}` is the previous input (see [set_applied_input](#set_applied_input)) |
| `soft_penalty` | `Scalar` | `1e4` | Quadratic weight on the state-bound slacks (see [Soft state constraints](#soft-state-constraints)) |
| `soft_state_penalty` | `optional<Vector<Scalar, NX>>` | none | Per-state slack weights; replaces `soft_penalty` when set |
| `terminal_constraint_set` | `optional<terminal_set<Scalar, NX>>` | none | Ellipsoidal or polyhedral set for `x_N`, always hard, validated at construction (see [Rejections](#rejections)); see [Terminal sets](#terminal-sets) for which one is invariant |
| `hard_state_constraints` | `bool` | `false` | When true, state constraints are hard (not softened) |

### Soft state constraints

The state bounds are softened when `x_min` or `x_max` is set and `hard_state_constraints` is false. The QP then carries one free slack `s_{k,i}` per state component and stage `k = 1 ... N`, couples it into the bound row as `x_min_i <= x_{k,i} - s_{k,i} <= x_max_i`, and puts the weight `rho_i` (`soft_state_penalty(i)` when set, else `soft_penalty`) on the slack's diagonal of the cost Hessian, with no linear term. The QP objective is `z^T P z / 2 + q^T z`, so a slack costs `rho_i s^2 / 2`: a quadratic penalty on the distance from `x_{k,i}` to its bounds.

The penalty is not exact. Where the hard problem is feasible and its bound carries the multiplier `mu`, the softened bound is still violated by about `mu / rho_i`, however large `rho_i` is. On `x+ = x + u`, `N = 1`, `Q = R = Qf = 1`, `x_max = 1`, tracking `y_ref = 10` from `x_0 = 0`, the hard problem has `mu = 8`, and the soft solution violates the bound by `8 / (2 + rho)`: 0.078 at `rho = 1e2` and 8.0e-4 at `rho = 1e4`. Set `hard_state_constraints = true` where a bound must hold.

## Construction

```cpp
using construction_error = mpc_construction_error<SetupError>;
// = std::variant<controller_construction_error, dare_error, SetupError>

static expected<mpc, construction_error> create(const discrete_state_space<Scalar, NX, NU, NY>& system,
                                                const mpc_config<Scalar, NX, NU, NY>& config);

static expected<mpc, construction_error> create(const discrete_state_space<Scalar, NX, NU, NY>& system,
                                                const mpc_config<Scalar, NX, NU, NY>& config,
                                                Solver solver);
```

`create` is the only construction path. Its constructor is public only so `ctrlpp::expected` can build the controller in the caller's storage, and it takes a passkey that only the class can make, so a rejected configuration or a refused solver setup is a value the caller has to inspect and never a controller that quietly stands in for one. It validates the configuration, computes the terminal cost (via DARE if `Qf` is not set), then builds the controller in place and sets its solver up against the posed QP. The two-argument form default-constructs the solver and chains into the three-argument form, which injects a caller-supplied, pre-configured solver (moved in before the initial QP is posed), letting you choose the accuracy/speed tradeoff via a preset.

`SetupError` is the error type of the bound solver's own `setup`, so the backend's cause of a setup failure reaches the caller unflattened. The error is one of three alternatives:

| Alternative | Carries |
|-------------|---------|
| `controller_construction_error` | a configuration refusal, in the order of the table below |
| `dare_error` | with `Qf` unset, the Riccati solve's own refusal of the terminal cost it supplies |
| `SetupError` | the solver refused the posed problem: `osqp_setup_error` for `osqp_solver`, `argmin_qp_setup_error` for `argmin_qp_solver`, or a custom backend's own enum |

The three alternative types must differ (a `static_assert` says so), because a variant whose alternatives repeat a type cannot be built from a value of that type.

```cpp
// Skip per-step polishing on the warm-resolve path -- see qp_preset.
auto created = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(
    sys, cfg, ctrlpp::osqp_solver{ctrlpp::qp_preset::speed});
if (!created)
{
    if (const auto* cause = std::get_if<ctrlpp::osqp_setup_error>(&created.error()))
        return handle_setup_refusal(*cause);
    if (const auto* cause = std::get_if<ctrlpp::dare_error>(&created.error()))
        return handle_riccati_refusal(*cause);
    return handle(std::get<ctrlpp::controller_construction_error>(created.error()));
}
auto controller = *std::move(created);
```

### Rejections

The configuration is validated once here, before any dimension product is formed and before any storage is reserved. Rejections, checked in order; every refusal but the last is returned before any controller exists:

| Condition | Alternative | Refusal |
|-----------|-------------|---------|
| `config.horizon <= 0` | `controller_construction_error` | `non_positive_horizon` |
| `config.horizon` above the representable bound of the derived dimensions | `controller_construction_error` | `horizon_overflow` |
| `u_min` or `u_max` holds a NaN, `u_min` holds `+inf`, `u_max` holds `-inf`, or `u_min > u_max` in some entry | `controller_construction_error` | `invalid_input_bounds` |
| the same on `x_min` and `x_max` | `controller_construction_error` | `invalid_state_bounds` |
| `du_max` holds a NaN or a negative entry | `controller_construction_error` | `invalid_rate_bound` |
| `Q`, then `R`, then `Qf` when set, holds a NaN or an infinity | `controller_construction_error` | `non_finite_cost_weight` |
| that weight is asymmetric beyond rounding drift | `controller_construction_error` | `asymmetric_cost_weight` |
| that weight is not positive semidefinite | `controller_construction_error` | `indefinite_cost_weight` |
| `soft_penalty`, or a set `soft_state_penalty`, is non-finite or non-positive in some entry while the state bounds are soft | `controller_construction_error` | `invalid_constraint_penalty` |
| a polytopic `terminal_constraint_set` has `H.rows() != h.size()` | `controller_construction_error` | `mismatched_terminal_set` |
| that polytope's `H` holds a NaN or an infinity, or its `h` holds a NaN | `controller_construction_error` | `non_finite_terminal_set` |
| an ellipsoidal `terminal_constraint_set` has a `P` that holds a NaN or an infinity | `controller_construction_error` | `non_finite_terminal_set` |
| that `P` is asymmetric beyond rounding drift | `controller_construction_error` | `asymmetric_terminal_set` |
| that `P` is not positive semidefinite | `controller_construction_error` | `indefinite_terminal_set` |
| that set's `alpha` is a NaN or an infinity | `controller_construction_error` | `non_finite_terminal_set` |
| that `alpha` is negative | `controller_construction_error` | `negative_terminal_set_level` |
| `system.A`, `system.B` or `system.C` holds a NaN or an infinity, whether or not `Qf` is set. `system.D` is not used by the controller and is not checked | `controller_construction_error` | `non_finite_system` |
| `C^T Q C`, then `C^T Qf C` when `Qf` is set, is not all finite: the weight mapped to state space overflows although every configured entry is finite | `controller_construction_error` | `non_finite_state_space_weight` |
| `Qf` is not set and `R` is singular to the Riccati solve that supplies the terminal cost, checked after every rule above | `controller_construction_error` | `singular_input_weight` |
| `Qf` is not set and the Riccati solve declines, checked after `singular_input_weight` | `dare_error` | the solve's own enumerator |
| the solver's `setup` refuses the posed problem, checked after every rule above | `SetupError` | `osqp_setup_error`, `argmin_qp_setup_error`, or a custom backend's own enum |

Each weight is checked completely (finite, symmetric, semidefinite) before the next one, so an indefinite `Q` is reported ahead of a non-finite `R`.

The terminal set is checked after the penalties and before the plant. An ellipsoid `{x : x^T P x <= alpha}` needs `P` positive semidefinite, not definite, under the weights' own rules: drift within their rounding bound is symmetrized in the controller's copy, so the posed rows use one symmetric `P`. The QP holds the ellipsoid as its inscribed box in `P`'s eigenbasis, `|v_i^T x| <= sqrt(alpha / (NX lambda_i))`. An eigenvalue within the admission's semidefiniteness floor of zero is a null direction of `P`, and the box is left open along it, which is exact for every `alpha >= 0`; `alpha = 0` pins the other directions to zero. In a polytope, an entry of `h` that is `+inf` leaves its row vacuous. An entry of `-inf` poses an empty set, which `create` does not test for: both shipped backends refuse its bounds at setup, so `create` returns their setup error (`osqp_setup_error::setup_failed`, `argmin_qp_setup_error::pose_failed`). Emptiness of finite rows is not tested either; it shows up as an infeasible solve.

Equal bounds pin a variable and are accepted, as is an infinite bound in its open direction (`u_min = -inf`, `u_max = +inf`). A rate bound of zero, `-0.0` or `+inf` is accepted.

The weights follow the library's symmetry convention. An asymmetry within the rounding bound `N * eps * max|W_ij|` (never below `N * denorm_min`) of an `N x N` weight `W` is drift: it is symmetrized in the controller's own copy of the configuration before the QP is assembled. The cost assembly passes only the upper triangle to the backend, so without that step the backend would see the reflection of the upper triangle rather than the symmetric part the configured cost denotes. A larger asymmetry is structure and is refused, never silently symmetrized. Semidefiniteness is decided by the symmetric eigensolver with the floor `N * eps * max|W_ij|`, the order of its backward error. `R` is required positive semidefinite in general: the QP needs only a semidefinite Hessian to be convex, so a singular `R` is accepted when `Qf` is set. With `Qf` unset, the Riccati solve that supplies the terminal cost inverts `R`, so `create` refuses a singular `R` there with `singular_input_weight` rather than build a controller with a different terminal cost. The verdict is the solve's own: the reciprocal-pivot rank test at `NU * eps` on the equilibrated weighting `R / s`, where `s` is the smallest power of two at or above `max(max|Q_ij|, max|R_ij|)`, with `Q` mapped to state space. It is a rank test, not a zero test, so `R = diag(1, 0)` is refused, and so is a weighting that is positive at its own scale but underflows to zero at the equilibrated one (`Q = 1e300 I`, `R = 1e-300`). It does not depend on `A`.

With `Qf` unset, no other terminal cost is ever posed in place of the Riccati solution: a solve that declines refuses the configuration with the solve's own `dare_error`. Setting `Qf` poses a terminal cost of the caller's own choosing; `Qf = Q` poses the stage weight. A rank-deficient `A`, such as an input delay carried as states or a deadbeat or FIR model, is refused as `singular_a`, because the symplectic method the Riccati solve uses inverts `A`. That is a limit of the method and not of the plant, so such a plant needs an explicit `Qf`. Of the remaining enumerators, `singular_r` never arrives, because `singular_input_weight` is checked first with the same test on the same operand; `non_finite_input` never arrives, because a non-finite plant, `R` or state-space weight is refused first; `schur_failed` has no known input that reaches it; and `non_psd_solution` arrives through the rounding allowance the semidefiniteness check grants `Q`. On a stable mode with eigenvalue `a` that `B` cannot reach, the Riccati solution carries that mode's weight amplified by `1 / (1 - a^2)`, while the solve's definiteness band scales only with `P`'s own magnitude, so a `Q` accepted at its floor can give a `P` the solve refuses. `A = diag(0.1, 0.8125)`, `B = [1; 0]`, `Q = diag(1, -2 eps)`, `R = 1` is such a pose: `p22 = -2 eps / (1 - 0.8125^2) = -512 eps / 87`, about `-5.89 eps`, lies beyond the band of about `4.02 eps` ([the definiteness check](../../guides/patterns/numerical-behavior.md#what-the-library-does-guard-against)). Both checks are kept: the refusal is the exact-arithmetic verdict on the `Q` that was admitted, and a band wide enough to accept that `P` would pose an indefinite terminal weight. Branch on `has_value()` and read the enumerator as a diagnostic; see [how far the enumerator can be relied on](../control/dare.md#how-far-the-enumerator-can-be-relied-on).

The state bounds are soft when `x_min` or `x_max` is set and `hard_state_constraints` is false; only then are the penalties inspected. `soft_penalty` is checked even while a per-state `soft_state_penalty` replaces it, because it is a configuration field with a valid default and an invalid value there is a mistake either way.

`horizon` stays a signed `std::int32_t` deliberately: a mistaken negative value remains representable as negative and is therefore rejectable, whereas an unsigned field would turn the same mistake into an enormous allocation. The overflow bound is a representability condition on `std::int32_t`, not a chosen ceiling. At the worst-case configuration the horizon `N` scales the decision vector as `N*(2*NX + NU) + NX` and the constraint rows as `N*(2*NX + 2*NU) + NX + n_terminal`, so both stay representable exactly when `horizon <= (INT32_MAX - (NX + n_terminal)) / (2*NX + 2*NU)`.

Two alternatives are deliberately not implemented. Validating at the first solve would surface a configuration error at the first control step, the worst possible moment. Clamping the horizon to one would turn a caller mistake into a silently different controller.

## Failure contract

Both `mpc` and `nmpc` share one failure contract. Every
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
  - `setup_incomplete`: never returned by `mpc`'s `solve`, because a refused
    solver setup refuses `create`. Among the solves, only an `nmpc` copy whose
    re-setup failed returns it.
  - `invalid_backend_result`: the solver reported a status the controller accepts
    and then returned a primal shorter than the decision dimension, or a dual
    shorter than the constraint count, of the problem posed at construction, or
    a non-finite entry in the consumed primal or dual prefix, a non-finite
    objective, solve time or residual, or a negative iteration count. The
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
ctrlpp::expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0);
```

Solves the QP for regulating state to the origin. Returns the success branch with
the first input and a soft status, or the error branch on failure. A controller
always holds a solver that accepted its problem at setup: a refused setup refuses
`create` with the `SetupError` alternative instead.

### solve (constant reference)

```cpp
ctrlpp::expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0,
                                                               const Vector<Scalar, NY>& y_ref);
```

Solves the QP for tracking a constant output reference across the entire horizon.
The stage cost becomes `(C x_k - y_ref)^T Q (C x_k - y_ref) / 2` plus the input
cost `u_k^T R u_k / 2`, still about `u = 0`. With `Qf` set the terminal cost is
`(C x_N - y_ref)^T Qf (C x_N - y_ref) / 2`; with `Qf` unset it is the DARE weight
`P` about `C^T (C C^T)^{-1} y_ref`, the minimum-norm state whose output is
`y_ref`. No steady-state target `(x_ss, u_ss)` is computed, so a reference whose
equilibrium needs a nonzero input settles with an offset; see
[Limitations](#limitations).

### solve (trajectory reference)

```cpp
ctrlpp::expected<solve_output<Scalar, NU>, solver_error> solve(const Vector<Scalar, NX>& x0,
                                                               std::span<const Vector<Scalar, NY>> y_ref);
```

Solves the QP for tracking a time-varying output reference trajectory, with
`y_ref[k]` in the stage cost of `x_k` and `y_ref[N]` in the terminal cost, each
weighted as for the constant reference. The span must contain at least
`horizon + 1` elements; an undersized span is rejected via the error branch
(`solver_error::invalid_problem`) rather than read out of bounds.

### set_applied_input

```cpp
void set_applied_input(const Vector<Scalar, NU>& u);
```

Records the control input the caller actually commanded. The internal previous
input (which anchors the rate constraint) starts at zero and is updated only on a
successful solve;
this accessor lets the caller keep that reference consistent when it commands a
different input, for example after an error branch or external saturation.

### trajectory

```cpp
ctrlpp::expected<std::pair<std::vector<Vector<Scalar, NX>>,
                           std::vector<Vector<Scalar, NU>>>,
                 solver_error> trajectory() const;
```

Returns the full predicted state and input trajectories from the last solve. It is
guarded: before the first valid solve it returns the error branch
(`solver_error::setup_incomplete`) rather than stale or default-initialized data.

### diagnostics

```cpp
mpc_diagnostics<Scalar> diagnostics() const;
```

Returns the diagnostics of the last solve, defined in `ctrlpp/mpc/diagnostics.h`. They are written on every solve that reaches the backend, failed ones included; before the first solve every field holds its default (`status = solve_status::error`, zeros).

| Field | What `mpc` writes |
|-------|-------------------|
| `status` | The backend's `solve_status`, or `invalid_backend_result` when the controller discards the answer (see [Failure contract](#failure-contract)) |
| `iterations` | The backend's iteration count |
| `solve_time` | The backend's solve time; always `0` on `argmin_qp_solver` |
| `cost` | The backend's QP objective, which omits the constant reference terms, so it can be negative while tracking |
| `primal_residual`, `dual_residual` | The backend's residuals |
| `max_constraint_violation` | The largest violation of the posed QP's rows at the returned primal `x`, `max_i max(l_i - (A x)_i, (A x)_i - u_i, 0)`, over the initial-state and dynamics equalities and the state, input, rate and terminal-set rows. Under soft state constraints the state rows are the relaxed rows `x_k - s_k`, so an overshoot the slack absorbs is not counted. A converged solve can report a positive value, since the backend stops within its tolerance. Written only when the answer is accepted, `0` otherwise |
| `max_path_constraint_violation` | The largest violation of the stage constraints at the returned primal: the state bounds on `x_1` to `x_N` as stated, so under soft state constraints an overshoot the slack absorbs is counted here, and the input and rate rows. Written only when the answer is accepted, `0` otherwise |
| `max_terminal_constraint_violation` | The largest violation of the terminal-set rows at the returned `x_N`, as the controller poses them: `H x_N <= h` for a polytope, and for an ellipsoid the inscribed box in `P`'s eigenbasis that stands in for it. `0` without a terminal set. Written only when the answer is accepted, `0` otherwise |
| `total_slack` | The sum of the absolute values of the state slack entries over the horizon, under soft state constraints; `0` under hard ones or without state bounds. Written only when the answer is accepted, `0` otherwise |
| `stop_criterion` | Always `nlp_stop_criterion::unreported` |

## Terminal sets

<!-- snippet: fragment -->
```cpp
#include <ctrlpp/mpc/invariant_set.h>

template <typename Scalar, std::size_t NX, std::size_t NU>
expected<terminal_ingredients_result<Scalar, NX, NU>, terminal_ingredients_error>
terminal_ingredients(const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NX, NU>& B,
                     const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R,
                     const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max,
                     const Vector<Scalar, NX>& x_min, const Vector<Scalar, NX>& x_max);

// in ctrlpp/mpc/terminal_set.h
using terminal_ingredients_error = std::variant<terminal_set_error, dare_error>;
```

`terminal_ingredients` solves the DARE for `P`, forms the LQR gain `K = -(R + B^T P B)^{-1} B^T P A`, and returns three fields. An infinite bound side is absent, so pass `+/-infinity` to leave a state or input unbounded on that side.

| Field | Type | What it is |
|-------|------|------------|
| `Qf` | `Matrix<Scalar, NX, NX>` | The DARE solution `P`, the terminal cost |
| `set` | `ellipsoidal_set<Scalar, NX>` | The ellipsoid `{x : x^T P x <= alpha}` with the largest `alpha` inside the input and state faces. The QP encodes its inscribed box in `P`'s eigenbasis, which is not positively invariant under `u = K x` and gives no recursive-feasibility guarantee. Only the input bounds are checked to straddle zero: a state bound that excludes the origin still caps `alpha` at a positive level, and the ellipsoid then lies outside that bound |
| `invariant_set` | `expected<polytopic_set<Scalar, NX>, terminal_set_error>` | The maximal output admissible set of `A + B K` inside `{x_min <= x <= x_max, u_min <= K x <= u_max}`: the invariant terminal set of the LQR law, as a minimal set of unit rows |

The function's own error channel carries the failures of `Qf` and `set`, as one of two alternatives. The `terminal_set_error` alternative is `singular_input_weight`, checked first with the same rank test as `mpc::create` and independently of `A`, or the ellipsoid's refusal, `input_zero_not_interior` or `empty_terminal_set`. The `dare_error` alternative is a Riccati solve that declines, carrying the solve's own enumerator, as `mpc::create` does with `Qf` unset. The output admissible set has its own outcome in its field, so its refusal never fails the call:

| `invariant_set` error | Cause |
|-----------------------|-------|
| `non_finite_input` | `A + B K` or `K` contains NaN or Inf, or a bound is NaN. An infinite bound side is absent, not corrupt |
| `origin_not_interior` | A finite state or input bound side does not strictly straddle zero |
| `closed_loop_not_stable` | `A + B K` is not strictly stable (the Stein equation `A_K^T X A_K - X = -I` has no positive definite solution) |
| `unbounded_terminal_set` | The rows `C A_K^k` of the first `NX` steps leave a direction in which the set recedes, so the polytope the construction starts from is unbounded. One-sided bounds that would close the set only at later steps are refused as well, conservatively |
| `not_converged` | No fixed point within the step bound derived from the Lyapunov contraction of `A + B K` |
| `certification_failed` | The constructed set failed its own check on its vertices (bounded, inside `Y`, mapped into itself by `A + B K`, each within the construction's rounding bounds), or a vertex enumeration returned fewer than `NX + 1` vertices, too few to describe a full-dimensional polytope. Plants with large gains can be refused this way, because the rounding error that accumulates in the output rows `C A_K^t` can exceed the certificate's bound |
| `dimension_unsupported` | `NX > 4` |

### terminal_set_error

Every terminal-set construction reports its refusals as one `terminal_set_error` (in `ctrlpp/mpc/terminal_set.h`):

| Enumerator | Meaning |
|------------|---------|
| `input_zero_not_interior` | Some input component has `u_min(i) >= 0` or `u_max(i) <= 0`, so `u = 0`, the LQR input at the origin, is not strictly interior to `[u_min, u_max]` and the input-face `alpha` bound would be meaningless |
| `empty_terminal_set` | Capping `alpha` against the input and state faces yields a non-positive or non-finite `alpha`: no consistent bounded ellipsoid exists. From `compute_polytopic_invariant_set`, the input polytope has no vertex |
| `halfplanes_truncated` | The polytopic invariant-set halfplane count hit its resource bound, so the representation would be a truncated, and therefore unsound, outer approximation |
| `not_converged` | The construction reached its iteration bound without a fixed point; no set is returned |
| `closed_loop_not_stable` | `A + B K` is not strictly stable (the Stein equation has no positive definite solution), so no bounded output admissible set exists |
| `unbounded_terminal_set` | The rows `C A_K^k` for `k < NX` leave a direction in which the set recedes, so the starting polytope is unbounded; one-sided bounds that would close the set only at later steps are refused as well, conservatively |
| `origin_not_interior` | A finite state or input bound does not strictly straddle zero |
| `dimension_unsupported` | The state dimension exceeds the four that vertex enumeration is restricted to |
| `certification_failed` | The constructed set failed its own check on its vertices, or a vertex enumeration returned fewer than `NX + 1` vertices; no set is returned |
| `non_finite_input` | `A_K` or `K` contains NaN or Inf, or a bound is NaN; an infinite bound is an absent side, not corrupt data |
| `singular_input_weight` | `R` is rank-deficient to the reciprocal-pivot test the Riccati solve applies to its equilibrated weighting, so the terminal cost and gain cannot be formed |

### compute_output_admissible_set

<!-- snippet: fragment -->
```cpp
#include <ctrlpp/mpc/output_admissible_set.h>

template <typename Scalar, std::size_t NX, std::size_t NU>
expected<polytopic_set<Scalar, NX>, terminal_set_error>
compute_output_admissible_set(const Matrix<Scalar, NX, NX>& A_K, const Matrix<Scalar, NU, NX>& K,
                              const Vector<Scalar, NU>& u_min, const Vector<Scalar, NU>& u_max,
                              const Vector<Scalar, NX>& x_min, const Vector<Scalar, NX>& x_max,
                              Scalar rounding_factor = Scalar{1},
                              std::int64_t max_steps = std::numeric_limits<std::int64_t>::max());
```

Builds the maximal output admissible set `O_inf = {x : C A_K^k x in Y for all k >= 0}` of `x+ = A_K x`, with `C = [I; K]` and `Y` the finite bound sides, for any closed loop (the `invariant_set` field calls it with `A_K = A + B K`). `A_K` must be the closed loop `A + B K` of the `K` passed: the input rows are formed from `K` and the dynamics from `A_K`, so a mismatched pair describes a law that is never applied. Its hypotheses are checked in this order: every entry of `A_K` and `K` is finite and no bound is NaN (else `non_finite_input`), every finite bound side strictly straddles zero (else `origin_not_interior`), `A_K` is strictly stable (else `closed_loop_not_stable`), and the rows `C A_K^k` of the first `NX` steps bound a polytope (else `unbounded_terminal_set`). `not_converged` is reported if the construction passes its derived step bound. The result is checked on its own vertices before it is returned, and an enumeration too small to describe a full-dimensional polytope is refused; either failure is `certification_failed`. `NX <= 4` is a compile-time requirement.

The finite determination of `O_inf` and its step bound need only a bounded starting polytope and the origin interior to `Y`; `Y` itself need not be bounded. One-sided bounds are therefore accepted whenever the rows of the first `NX` steps close the set, as two-sided input bounds do whenever `(K, A)` is observable. Whether one-sided rows that leave a direction open close the set at some later step is a positivity question about the sequence `C A_K^k d`, which has no finite test in general, so such a configuration is refused conservatively. The recession test compares each row with each candidate direction within that comparison's rounding bound plus the direction's own forward error, so its tolerance leans toward refusal; rows that are linearly dependent to working precision define no candidate direction.

The construction stacks the rows `C A_K^k` until they bound the set, then appends, step by step, the rows that the current set's vertices do not imply, prunes redundant rows every step, and stops at the first step that appends nothing. There is no linear program: every implication and facet test runs on enumerated vertices, which costs `C(m, NX)` small solves per step for `m` rows. Each comparison of a row with a vertex carries its own rounding bound (the residual as an inner product, plus the vertex's forward error from the residual of its defining rows); `rounding_factor` scales the order-one prefactor of every such bound and must be positive. A row declared implied within its bound can leave the result larger than the exact set by at most that bound; a looser feasibility test only keeps more rows. A row is kept as a facet when some `NX` of the vertices on it, projected onto its hyperplane, span `NX - 1` dimensions by more than their forward errors can move a singular value, so a vertex that is located only roughly, such as a near-duplicate at a degenerate corner, cannot remove a facet. On the double integrator with `dt = 0.1`, `Q = I`, `R = 1`, `|u| <= 1` and `|x_i| <= 10`, the set has 30 rows, the last needed at step 14.

`max_steps` bounds the construction's wall time. The derived step bound grows with the Lyapunov matrix of `A_K`, so a closed loop near the unit circle can need very many steps, each enumerating `C(m, NX)` subsets. A step index `t > max_steps` that still appends rows is refused with `not_converged`; the rows of the first `NX` steps are always stacked, so a value below `NX` refuses as soon as a later step is needed. The default leaves only the derived bound, so omitting the argument gives the same result. There is no row budget, because rows are pruned to facets every step, so the step budget already bounds the work of each step. On the case above, `max_steps = 14` returns the same 30 rows and `max_steps = 13` is refused. `terminal_ingredients` runs the construction with the default budget; a caller who needs to bound the time calls `compute_output_admissible_set` directly.

### compute_polytopic_invariant_set

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
expected<polytopic_set<Scalar, NX>, terminal_set_error>
compute_polytopic_invariant_set(const Matrix<Scalar, NX, NX>& A_sys, const Matrix<Scalar, NX, NU>& B_sys,
                                const polytopic_set<Scalar, NX>& state_constraints,
                                const polytopic_set<Scalar, NU>& input_constraints,
                                std::int32_t max_iterations = 100, Scalar convergence_tol = Scalar{1e-6});
```

Declared in `ctrlpp/mpc/invariant_set.h`. It iterates the robust pre-image `{x : A x + B u in S for every vertex u of the input polytope}`, intersected with the current set, starting from the state polytope: the result is the set of states that stay inside the state polytope under every input in the input polytope. That is a robust control-invariant set, a valid object for robust MPC; it is not the terminal set of the LQR law. Use `invariant_set` for that.

A row is dropped only when it is parallel to a kept row and no tighter, within `convergence_tol`; no linear program removes other redundant rows. The iteration stops when a step leaves the row count unchanged. It refuses with `empty_terminal_set` when the input polytope has no vertex, `halfplanes_truncated` when more than 500 rows survive the filter, and `not_converged` after `max_iterations` steps. `NX <= 4` is a compile-time requirement.

### Using the invariant terminal set

```cpp
auto ingredients = ctrlpp::terminal_ingredients<double, NX, NU>(A, B, Q, R, u_min, u_max, x_min, x_max);
if(!ingredients)
    return report(ingredients.error());
if(!ingredients->invariant_set)
    return report(ingredients->invariant_set.error());

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

The outer channel carries a `terminal_ingredients_error`, the `std::variant` of a `terminal_set_error` (a singular `R`, or the ellipsoid's refusal) and the Riccati solve's own `dare_error`; the `invariant_set` channel carries a `terminal_set_error`. Either way a refusal names its cause; `report` stands for the caller's handling, which in the example visits the variant, prints the enumerator's name through an exhaustive `switch` for each error type, and exits. `terminal_ingredients` takes the state weight, while `mpc_config`'s `Q` and `Qf` are output weights posed as `C^T Q C` and `C^T Qf C`; the two agree here because the plant has `C = I`. Recursive feasibility and stability of the origin follow from this configuration only with hard state constraints, a nominal model and regulation to the origin; see [background/mpc](../../background/mpc.md). The full closed loop is in `examples/mpc/ctrlpp_mpc_05_terminal_set.cpp`.

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
        (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(),
        Eigen::Matrix2d::Identity(),
        Eigen::Matrix<double, 2, 1>::Zero()};

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 20,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-1.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(1.0)};

    auto created = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(sys, cfg);
    if(!created)
    {
        // created.error() carries the construction_error: a configuration
        // refusal, the Riccati solve's refusal, or the solver's setup refusal.
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

## Limitations

### No offset-free formulation exists

The controller takes the state it is given and predicts with the nominal `(A, B)`. The MPC code has no disturbance model, no augmented disturbance state and no coupling to an estimator, so a constant disturbance or a model mismatch leaves a steady-state offset. Regulating `x+ = x + u` at `horizon = 20` with the default weights while the plant adds a constant `0.1` per step settles at `x = 0.162`, not 0. An offset-free scheme needs the caller's own disturbance estimate folded into `x0` or the reference.

### Tracking has no steady-state input target

The tracking overloads penalize the input about `u = 0` and compute no `(x_ss, u_ss)` with `C x_ss = y_ref`, so a reference whose equilibrium needs a nonzero input settles short of it, by an amount that grows with `R` relative to `Q`. On `x+ = 0.9 x + 0.1 u`, `C = 1`, `horizon = 20` and the default weights, tracking `y_ref = 1` settles at `y = 0.489`. The double integrator's equilibria need `u = 0`, so it tracks without this offset.

### The softened state bounds are not an exact penalty

The default softening is a quadratic penalty on a free slack, which violates an active bound by about the multiplier over the weight even when the hard problem is feasible; see [Soft state constraints](#soft-state-constraints). Only `hard_state_constraints = true` poses the bounds as hard constraints, which then hold to the solver's tolerance.

### Each terminal-set function computes a different set

`compute_output_admissible_set` builds the maximal output admissible set of a given closed loop, and `terminal_ingredients` fills its `invariant_set` field with that set for `A + B K` only when `NX <= 4`; for larger `NX` the field holds `dimension_unsupported`. `compute_polytopic_invariant_set` builds a robust control-invariant set (every input in the input polytope keeps the state inside), not the LQR law's invariant set, and also requires `NX <= 4`. The ellipsoid in `set` is encoded as an inscribed box that is not invariant under the LQR law. Only the maximal output admissible set, posed hard with hard state constraints, carries the recursive-feasibility argument; see [the terminal constraint set](../../background/mpc.md#terminal-constraint-set-mathcalx_f).

## See Also

- [nmpc](nmpc.md)<br/> nonlinear MPC
- [osqp-solver](osqp-solver.md)<br/> OSQP QP solver backend
- [argmin-qp-solver](argmin-qp-solver.md)<br/> header-only argmin QP backend
- [mhe](mhe.md)<br/> linear moving horizon estimation
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
- [guides/intro/your-first-mpc](../../guides/intro/your-first-mpc.md)<br/> introductory MPC tutorial
- [background/mpc](../../background/mpc.md)<br/> MPC theory and background
