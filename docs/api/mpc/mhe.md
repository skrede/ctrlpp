# mhe

> **Requires a QP backend.** The estimator itself needs only Eigen. The shipped QP backends are
> [osqp_solver](osqp-solver.md) (`-DCTRLPP_BUILD_OSQP=ON`) and [argmin_qp_solver](argmin-qp-solver.md)
> (`-DCTRLPP_BUILD_ARGMIN=ON`).

Moving horizon estimator that poses its window as one sparse quadratic program. Each window solve linearizes the
dynamics and measurement models once, at the window's middle sample, and estimates every state in the window from
an arrival cost, a process-noise cost and a measurement cost, subject to optional state bounds (hard or softened)
and measurement residual bounds. An embedded `ekf` runs alongside it: it supplies the arrival cost, produces the
estimate while the window fills, and takes over whenever a window solve fails. The estimator is exact for linear
models; for a nonlinear or affine model, read [Limitations](#limitations) first, or use [nmhe](nmhe.md).

## Headers

| Header | Provides |
|--------|----------|
| `#include <ctrlpp/mhe.h>` | `mhe` |
| `#include <ctrlpp/mhe/mhe_config.h>` | `mhe_config`, `nmhe_config`, `moving_horizon_configuration_error`, `moving_horizon_construction_error` |
| `#include <ctrlpp/mhe/mhe_diagnostics.h>` | `mhe_diagnostics` |
| `#include <ctrlpp/mhe/mhe_qp_formulation.h>` | the window QP builder (in `ctrlpp::detail`; see [Formulation](#formulation)) |

`ctrlpp/mhe.h` includes the other three.

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU,
          std::size_t NY, std::size_t N,
          typename Solver, typename Dynamics, typename Measurement>
    requires qp_solver<Solver>
          && dynamics_model<Dynamics, Scalar, NX, NU>
          && measurement_model<Measurement, Scalar, NX, NY>
class mhe;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` (through the embedded `ekf`), equal to `Solver::scalar_type` | Numeric type; both shipped QP backends are `double` |
| `NX` | `>= 1` | State dimension |
| `NU` | `>= 1` | Input dimension |
| `NY` | `>= 1` | Measurement dimension |
| `N` | `>= 1` | Window length: the window holds `N + 1` states and measurements and `N` inputs |
| `Solver` | satisfies `qp_solver`, default-constructible | QP backend: [osqp_solver](osqp-solver.md) or [argmin_qp_solver](argmin-qp-solver.md) |
| `Dynamics` | satisfies `dynamics_model<Scalar, NX, NU>` | Discrete-time dynamics callable |
| `Measurement` | satisfies `measurement_model<Scalar, NX, NY>` | Measurement callable |

`NX`, `NU` and `NY` are checked by the embedded `ekf`'s static assertions, `N` by the estimator's and the
configuration's. The estimator constructs its `Solver` with the type's default constructor; no solver instance or
setting can be passed in.

## mhe_config

A configuration aggregate `mhe_config<Scalar, NX, NU, NY, N>` passed to `create`; set only the fields that differ
from the defaults.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | Process noise covariance; its inverse weights the process-noise cost |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance; its inverse weights the measurement cost |
| `x0` | `Vector<Scalar, NX>` | Zero | Initial state of the embedded EKF and of the window |
| `P0` | `Matrix<Scalar, NX, NX>` | Identity | Initial covariance of the embedded EKF |
| `arrival_cost_weight` | `Scalar` | `1` | Factor on the arrival cost's inverse covariance |
| `x_min` | `optional<Vector<Scalar, NX>>` | none | Lower bound on every state in the window |
| `x_max` | `optional<Vector<Scalar, NX>>` | none | Upper bound on every state in the window |
| `residual_bound` | `optional<Vector<Scalar, NY>>` | none | Hard bound `\|z_k - H x_k\| <= residual_bound` on every measurement in the window |
| `soft_constraints` | `bool` | `true` | Soften the state bounds with non-negative slacks; the residual bound is never softened |
| `soft_penalty` | `Scalar` | `1e4` | Linear penalty per unit of each slack (an L1 penalty, since the slacks are non-negative) |
| `numerical_eps` | `Scalar` | `cbrt(eps)` | Finite-difference step for the Jacobians of a model that is not analytically differentiable, used by the window linearization and the embedded EKF |

Only one of `x_min` and `x_max` may be given; the missing side is unbounded.

## Creation

```cpp
static ctrlpp::expected<mhe, moving_horizon_construction_error> create(
    Dynamics dynamics,
    Measurement measurement,
    const mhe_config<Scalar, NX, NU, NY, N>& config);
```

Validates the configuration, then constructs the embedded EKF, the windows and the warm start. Every window starts
at `x0` and every arrival prior at `(x0, P0)`.

`moving_horizon_construction_error` is `std::variant<filter_error, moving_horizon_configuration_error>`. The checks
run in this order and the first failure is returned:

1. The embedded `ekf`'s own validation of `Q`, `R`, `x0` and `P0`, reported as its `filter_error`
   (see [ekf](../estimation/ekf.md)).
2. The `moving_horizon_configuration_error` checks:

| Enumerator | Condition |
|------------|-----------|
| `non_invertible_process_noise` | `Q` is singular under a full-pivot LU, or its inverse is not finite |
| `non_invertible_measurement_noise` | the same for `R` |
| `non_invertible_initial_covariance` | the same for `P0` |
| `non_positive_arrival_cost_weight` | `arrival_cost_weight` is not finite or not positive |
| `invalid_state_bounds` | a given bound is not finite, or some element of `x_min` exceeds `x_max` |
| `invalid_residual_bound` | `residual_bound` has a non-finite or negative element |
| `non_positive_soft_penalty` | `soft_penalty` is not finite or not positive (checked even when nothing is softened) |
| `non_positive_numerical_eps` | `numerical_eps` is not finite or not positive |
| `invalid_path_penalty` | `nmhe` only |
| `non_finite_path_constraint` | `nmhe` only |

The window length `N` is a template parameter on both the class and `mhe_config`, so `N == 0` fails at compile time.
Callers must branch on the returned `expected` before moving out the estimator.

## Formulation

At each window solve the estimator builds the QP `min 1/2 z' P z + q' z` subject to `l <= A z <= u` and hands it to
the backend's `setup`, then `solve`s it warm-started from the previous solution. The decision vector is
`z = [x_0 .. x_N]`, followed by the lower and upper slacks `[s_lo_0 .. s_lo_N, s_up_0 .. s_up_N]` when the state
bounds are softened. With `(A, B, H)` the single linearization described under [Limitations](#limitations), the
objective is

```
1/2 w |x_0 - xbar|^2_{Pbar^-1} + 1/2 sum_{k<N} |x_{k+1} - A x_k - B u_k|^2_{Q^-1} + 1/2 sum_{k<=N} |z_k - H x_k|^2_{R^-1}
  + soft_penalty * sum (s_lo + s_up) + 1/2 * 1e-6 * |s|^2
```

where `w` is `arrival_cost_weight` and `(xbar, Pbar)` is the embedded EKF's predicted state and covariance at the
window's first sample, captured before that sample's measurement correction so the measurement is not counted
twice. The dynamics enter only through the process-noise cost; there are no equality rows. The constraint rows are
`x_min <= x_k <= x_max` (hard), or `x_k + s_lo_k >= x_min`, `x_k - s_up_k <= x_max` and `s >= 0` (softened), and
`z_k - residual_bound <= H x_k <= z_k + residual_bound`, for every `k` in the window.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates the embedded EKF one step and shifts `u` into the input window.

### update

```cpp
ctrlpp::expected<void, ekf_update_error> update(const output_vector_t& z);
```

Incorporates a new measurement, or reports why it could not. The embedded EKF is corrected first and always runs:
its state and covariance are never overwritten by the window estimate. The first `N` accepted measurements after
construction (or after a refusal) use the EKF's estimate; from the `N + 1`-th on, when the window holds `N + 1`
measurements, every update solves the window QP. The QP is rebuilt and the backend's `setup` runs again on every
window solve.

A window solve falls back to the EKF's estimate, and says so in `diagnostics()`, when:

| Cause | `diagnostics().status` |
|-------|------------------------|
| The arrival covariance is singular under a full-pivot LU or its inverse is not finite | `solve_status::error` |
| The backend's `setup` returns an error | `solve_status::error` |
| The solve returns any status other than `optimal` or `solved_inaccurate` | `solve_status::error` |
| An accepted solve returns a primal shorter than the decision vector, a dual shorter than the constraint rows, a non-finite consumed value or diagnostic, or a negative iteration count | `solve_status::invalid_backend_result` |

The three causes reported as `solve_status::error` cannot be told apart from the diagnostics.

**Two channels, and which one carries what.** A measurement the embedded filter cannot use is a FAILURE: the step did not happen, and `update` returns the filter's own `ekf_update_error` verbatim rather than restating the same three conditions under a second name that could drift from the one that decides. Everything else is a DISPOSITION on a step that succeeded -- which of the two estimators produced the estimate, and how the solve went -- and that is what `diagnostics()` carries.

A refusal cannot undo the preceding prediction: that input was applied to the plant, so `state()` exposes the embedded filter's predicted prior. The missing measurement does invalidate the fixed-step horizon. The estimator clears its input, measurement, prior, and warm-start windows, sets `is_initialized()` to false, and uses the embedded EKF for the next N accepted measurements. Optimization resumes only after a coherent horizon has been refilled. `diagnostics()` continues to describe the last successful update until the next accepted measurement; the failed `update` return is what identifies the current state as an uncorrected prior.

### state

```cpp
const state_vector_t& state() const;
```

Returns the current state estimate, the last element of the window: the window solve's `x_N`, or the EKF's estimate
during fill-up and after a fallback.

### covariance

```cpp
const cov_matrix_t& covariance() const;
```

Returns the embedded EKF's covariance. No covariance of the window estimate is computed.

### innovation

```cpp
const output_vector_t& innovation() const;
```

After a window solve, the residual `z_N - h(x_N)` of the newest measurement against the window's final estimate.
During fill-up and after a fallback, the EKF's innovation, the newest measurement against its prediction.

### trajectory

```cpp
std::span<const state_vector_t> trajectory() const;
```

Returns the `N + 1` window states of the last window solve, oldest first. During fill-up the slots not yet written
hold `x0` and the others hold the EKF's estimates; a fallback replaces only the last element.

### arrival_state and arrival_covariance

```cpp
const state_vector_t& arrival_state() const;
const cov_matrix_t& arrival_covariance() const;
```

`arrival_state()` returns the first element of the window, the estimate of the window's oldest state (not the
arrival prior the cost is anchored to). `arrival_covariance()` returns the arrival prior's covariance, the EKF's
predicted covariance at the window's first sample.

### diagnostics

```cpp
const mhe_diagnostics<Scalar>& diagnostics() const;
```

Returns the `mhe_diagnostics<Scalar>` of the last accepted update. This is the DISPOSITION channel and it describes
a step that succeeded. A refused measurement never appears here, because no estimate was produced for the aggregate
to describe -- it is returned by `update` instead, and the aggregate goes on describing the last accepted step.

| Field | Type | After a window solve | When the EKF produced the estimate |
|-------|------|----------------------|------------------------------------|
| `status` | `solve_status` | the backend's status (`optimal` or `solved_inaccurate`) | `optimal` during fill-up; the fallback reason otherwise |
| `iterations` | `std::int32_t` | the backend's iteration count | `0` |
| `solve_time` | `Scalar` | the backend's solve time | `0` |
| `cost` | `Scalar` | the backend's objective | `0` |
| `primal_residual` | `Scalar` | the backend's primal residual | `0` |
| `dual_residual` | `Scalar` | the backend's dual residual | `0` |
| `max_constraint_violation` | `Scalar` | `max(primal_residual, dual_residual)`, not a measured constraint violation | `0` |
| `max_residual_bound_violation` | `Scalar` | never computed: always `0` | `0` |
| `total_slack` | `Scalar` | sum of the slacks when the state bounds are softened, else `0` | `0` |
| `used_ekf_fallback` | `bool` | `false` | `true` |

A default-constructed `mhe_diagnostics` has `status == solve_status::error`, every number zero and
`used_ekf_fallback == false`; that is what `diagnostics()` returns before the first update.

### is_initialized

```cpp
bool is_initialized() const;
```

Returns `true` once `N + 1` measurements have been accepted since construction or the last refusal, that is, once
updates solve the window.

## Usage Example

```cpp
// gnuplot: plot "< ./mhe_estimation" using 1:4 with lines title "estimate"
#include <ctrlpp/mhe.h>
#include <ctrlpp/mpc/osqp_solver.h>

#include <Eigen/Dense>

#include <iostream>
#include <random>
#include <utility>

struct constant_dynamics
{
    Eigen::Vector2d operator()(const Eigen::Vector2d& x,
                               const Eigen::Matrix<double, 1, 1>&) const
    {
        return x;
    }
};

struct position_measurement
{
    Eigen::Matrix<double, 1, 1> operator()(const Eigen::Vector2d& x) const
    {
        return Eigen::Matrix<double, 1, 1>{x[0]};
    }
};

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY = 1;
    constexpr std::size_t N = 10;

    ctrlpp::mhe_config<double, NX, NU, NY, N> cfg{
        .Q = Eigen::Matrix2d::Identity() * 0.01,
        .R = Eigen::Matrix<double, 1, 1>::Identity() * 1.0,
        .P0 = Eigen::Matrix2d::Identity() * 10.0};

    using estimator_type =
        ctrlpp::mhe<double, NX, NU, NY, N,
                    ctrlpp::osqp_solver,
                    constant_dynamics,
                    position_measurement>;
    auto estimator_result = estimator_type::create(
        constant_dynamics{}, position_measurement{}, cfg);
    if(!estimator_result)
    {
        std::cerr << "invalid MHE configuration\n";
        return 1;
    }
    auto estimator = std::move(*estimator_result);

    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 1.0);

    Eigen::Vector2d x_true(3.0, -1.0);
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < 30; ++k)
    {
        estimator.predict(u);

        Eigen::Matrix<double, 1, 1> z;
        z << x_true[0] + noise(rng);
        if(!estimator.update(z))
            continue; // predicted prior retained; horizon refill begins

        auto x_hat = estimator.state();
        std::cout << "k=" << k
                  << "  true=[" << x_true.transpose() << "]"
                  << "  est=[" << x_hat.transpose() << "]"
                  << "  init=" << estimator.is_initialized() << "\n";
    }
}
```

## Limitations

### The window shares one linearization point

Each window solve linearizes the dynamics and the measurement once, at the window's middle index: `A` and `B` are the
Jacobians of `Dynamics` at `(x_ref, u_ref)` and `H` the Jacobian of `Measurement` at `x_ref`, analytic when the
model provides them and finite differences with `numerical_eps` otherwise. `x_ref` is the element at index `N / 2`
of the trajectory the previous solve left (on the first solve, what the fill-up wrote there), and `u_ref` is the input at index
`N / 2` of the current window; the trajectory is not shifted before the linearization, so `x_ref` is one sample older
than `u_ref`. Every stage of the window then uses the same `(A, B, H)`. A strongly nonlinear model over a long window
is approximated by that one linearization, and the model is never re-linearized within a solve. [nmhe](nmhe.md)
solves the full nonlinear problem instead.

### The linearized model drops the affine offset

The QP models the window as `x_{k+1} = A x_k + B u_k` and `z_k = H x_k`: the constant terms of the linearization,
`f(x_ref, u_ref) - A x_ref - B u_ref` and `h(x_ref) - H x_ref`, are left out. The estimate is exact only when the
models are linear without an offset. Measured with `NX = NU = NY = 1`, `N = 10`, `Q = 0.01`, `R = 1`, `P0 = 1`,
noise-free measurements and the state started at its true value: with `z = x + 5` and a constant true state of `1`,
the window estimate is `5.67` on the first solve and `4.34` after 40 steps, while an `ekf` on the same models stays at
`1`; with the drift `x+ = x + 0.1` and `z = x`, the window estimate is `3.40` against a true `4.00` after 40 steps.

### The solver cannot be configured

`mhe` default-constructs its `Solver` and offers no way to pass a configured instance, so the backend runs at its
default settings; `osqp_solver`, for example, runs with polishing on. Other settings take a wrapper type, named as
`Solver`, whose default constructor builds the backend with them and whose `setup` and `solve` forward to it. With
the OSQP 1.0.0 library, a polish that finds no active set prints
`Polishing not needed - no active set detected at optimal point` to standard output even though `verbose` is off,
once per window solve; the usage example above prints it.

## See Also

- [nmhe](nmhe.md)<br/> nonlinear moving horizon estimation
- [osqp-solver](osqp-solver.md)<br/> OSQP QP solver backend
- [argmin-qp-solver](argmin-qp-solver.md)<br/> argmin QP solver backend
- [ekf](../estimation/ekf.md)<br/> the embedded extended Kalman filter
- [dynamics-model](../model/dynamics-model.md)<br/> dynamics model concept
- [measurement-model](../model/measurement-model.md)<br/> measurement model concept
- [background/mhe](../../background/mhe.md)<br/> MHE theory and background
