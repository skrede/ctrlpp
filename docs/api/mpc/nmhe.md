# nmhe

> **Requires an NLP backend.** The estimator itself needs only Eigen. The shipped NLP backends are
> [nlopt_solver](nlopt-solver.md) (`-DCTRLPP_BUILD_NLOPT=ON`) and [argmin_solver](argmin-solver.md)
> (`-DCTRLPP_BUILD_ARGMIN=ON`).

Nonlinear moving horizon estimator that poses its window as one nonlinear program. Every state in the window is a
decision variable, and the nonlinear dynamics and measurement models enter the cost directly, so no linearization
is made. The window is estimated from an arrival cost, a process-noise cost and a measurement cost, subject to
optional hard state bounds, hard measurement residual bounds and a per-stage nonlinear path constraint that can be
softened. An embedded `ekf` runs alongside it: it supplies the arrival cost, produces the estimate while the window
fills, and takes over whenever a window solve fails.

## Headers

| Header | Provides |
|--------|----------|
| `#include <ctrlpp/nmhe.h>` | `nmhe` |
| `#include <ctrlpp/mhe/mhe_config.h>` | `nmhe_config`, `moving_horizon_configuration_error`, `moving_horizon_construction_error` |
| `#include <ctrlpp/mhe/mhe_diagnostics.h>` | `mhe_diagnostics` |
| `#include <ctrlpp/mhe/mhe_nlp_formulation.h>` | `nmhe_formulation_state` and the window NLP builder (see [Formulation](#formulation)) |

`ctrlpp/nmhe.h` includes the other three.

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU,
          std::size_t NY, std::size_t N,
          nlp_solver Solver,
          typename Dynamics, typename Measurement,
          std::size_t NC = 0>
    requires dynamics_model<Dynamics, Scalar, NX, NU>
          && measurement_model<Measurement, Scalar, NX, NY>
class nmhe;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | `ctrlpp_floating_scalar` (through the embedded `ekf`), equal to `Solver::scalar_type` | Numeric type; `nlopt_solver` is `double` only |
| `NX` | `>= 1` | State dimension |
| `NU` | `>= 1` | Input dimension |
| `NY` | `>= 1` | Measurement dimension |
| `N` | `>= 1` | Window length: the window holds `N + 1` states and measurements and `N` inputs |
| `Solver` | satisfies `nlp_solver`, default-constructible | NLP backend: `nlopt_solver<double>` or `argmin_solver<Scalar, Policy>` |
| `Dynamics` | satisfies `dynamics_model<Scalar, NX, NU>` | Discrete-time dynamics callable |
| `Measurement` | satisfies `measurement_model<Scalar, NX, NY>` | Measurement callable |
| `NC` | `>= 0` | Number of path constraint rows per stage (default 0) |

The estimator constructs its `Solver` with the type's default constructor; no solver instance or setting can be
passed in. `nlopt_solver<double>` therefore runs SLSQP at its default settings. The window problem has inequality
rows only (the path and residual rows) and variable bounds. On the argmin side, a bound-only policy
(`argmin_lbfgsb`) solves without the rows, and `argmin_auglag` folds them into its augmented objective; see
[argmin-solver](argmin-solver.md). The SQP policies take them directly.

## nmhe_config

A configuration aggregate `nmhe_config<Scalar, NX, NU, NY, N, NC>` passed to `create`; set only the fields that
differ from the defaults.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | Process noise covariance; its inverse weights the process-noise cost |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance; its inverse weights the measurement cost |
| `x0` | `Vector<Scalar, NX>` | Zero | Initial state of the embedded EKF and of the window |
| `P0` | `Matrix<Scalar, NX, NX>` | Identity | Initial covariance of the embedded EKF |
| `arrival_cost_weight` | `Scalar` | `1` | Factor on the arrival cost's inverse covariance |
| `x_min` | `optional<Vector<Scalar, NX>>` | none | Hard lower bound on every state in the window |
| `x_max` | `optional<Vector<Scalar, NX>>` | none | Hard upper bound on every state in the window |
| `residual_bound` | `optional<Vector<Scalar, NY>>` | none | Hard bound `\|z_k - h(x_k)\| <= residual_bound` on every measurement in the window |
| `soft_constraints` | `bool` | `true` | Soften the path constraint with non-negative slacks; the bounds are never softened |
| `soft_penalty` | `Scalar` | `1e4` | Validated by `create`, otherwise unused by `nmhe`; the path slacks are weighted by `path_penalty` |
| `numerical_eps` | `Scalar` | `cbrt(eps)` | Finite-difference step of the embedded EKF's Jacobians for a model that is not analytically differentiable |
| `path_constraint_cb` | `optional<function<Vector<Scalar, NC>(const Vector<Scalar, NX>&)>>` | none | Path constraint `g(x_k) <= 0`, imposed on every state in the window; used only when `NC > 0` |
| `path_penalty` | `Vector<Scalar, NC>` | `1e4` each | Linear penalty per unit of each path slack |

## Creation

```cpp
static ctrlpp::expected<nmhe, moving_horizon_construction_error> create(
    Dynamics dynamics,
    Measurement measurement,
    const nmhe_config<Scalar, NX, NU, NY, N, NC>& config);
```

Validates the configuration, then constructs the embedded EKF, the windows and the warm start, builds the window
NLP and runs the backend's `setup` on it, once. A `setup` failure is not reported by `create`; see
[Limitations](#limitations).

The checks run in the same order as `mhe`'s and return the same errors, listed under [mhe](mhe.md#creation), followed
by two that only `nmhe` makes, when `NC > 0` and `path_constraint_cb` is set:

| Enumerator | Condition |
|------------|-----------|
| `invalid_path_penalty` | `path_penalty` has a non-finite or non-positive element (checked even when `soft_constraints` is false) |
| `non_finite_path_constraint` | `path_constraint_cb(x0)` returns a non-finite element |

The window length `N` is a template parameter on both the class and `nmhe_config`, so `N == 0` fails at compile
time. Callers must branch on the returned `expected` before moving out the estimator.

## Formulation

The decision vector is `[x_0 .. x_N]`, followed by the path slacks `[s_0 .. s_N]`, `NC` per stage, when the path
constraint is softened. The objective is

```
1/2 w |x_0 - xbar|^2_{Pbar^-1} + 1/2 sum_{k<N} |x_{k+1} - f(x_k, u_k)|^2_{Q^-1} + 1/2 sum_{k<=N} |z_k - h(x_k)|^2_{R^-1}
  + sum_{k<=N} path_penalty . s_k
```

where `w` is `arrival_cost_weight` and `(xbar, Pbar)` is the embedded EKF's predicted state and covariance at the
window's first sample, captured before that sample's measurement correction so the measurement is not counted
twice. The dynamics enter only through the process-noise residual `x_{k+1} - f(x_k, u_k)` in the cost; there are no
continuity equalities. The variable bounds are `x_min <= x_k <= x_max` and `s_k >= 0`. The constraint rows, every
one an upper bound of zero, are `g(x_k) - s_k` (or `g(x_k)` when not softened) and `z_k - h(x_k) - residual_bound`
and `-(z_k - h(x_k)) - residual_bound`, for every `k` in the window.

The cost gradient is a central difference of the whole cost, with the step `cbrt(eps) * max(1, |z_j|)` per
coordinate, and no constraint Jacobian is supplied, so the backend finite-differences the rows. The problem's
callbacks read the arrival prior and the input and measurement windows through a shared
`nmhe_formulation_state<Scalar, NX, NU, NY, N>`, which the estimator writes before each solve:

| Field | Type | Holds |
|-------|------|-------|
| `arrival_state` | `Vector<Scalar, NX>` | the arrival prior `xbar` |
| `arrival_P_inv` | `Matrix<Scalar, NX, NX>` | the inverse of `Pbar` (of `P0` before the first solve) |
| `u_window` | `std::array<Vector<Scalar, NU>, N>` | the window's inputs |
| `z_window` | `std::array<Vector<Scalar, NY>, N + 1>` | the window's measurements |

Each solve warm-starts from the previous solution shifted by one sample, with the EKF's estimate in the newest slot
and the newest path slacks at zero.

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
measurements, every update solves the window NLP.

A window solve falls back to the EKF's estimate, and says so in `diagnostics()`, when:

| Cause | `diagnostics().status` |
|-------|------------------------|
| The backend's `setup` failed when the estimator was created | `solve_status::error` |
| The arrival covariance is singular under a full-pivot LU or its inverse is not finite | `solve_status::error` |
| The solve returns any status other than `optimal` or `solved_inaccurate` | `solve_status::error` |
| An accepted solve returns a decision vector shorter than the NLP's, a non-finite consumed value or diagnostic, or a negative iteration count | `solve_status::invalid_backend_result` |

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

Returns the [`mhe_diagnostics<Scalar>`](mhe.md#diagnostics) of the last accepted update. This is the DISPOSITION
channel and it describes a step that succeeded. A refused measurement never appears here, because no estimate was
produced for the aggregate to describe -- it is returned by `update` instead, and the aggregate goes on describing
the last accepted step.

| Field | After a window solve | When the EKF produced the estimate |
|-------|----------------------|------------------------------------|
| `status` | the backend's status (`optimal` or `solved_inaccurate`) | `optimal` during fill-up; the fallback reason otherwise |
| `iterations`, `solve_time`, `cost`, `primal_residual` | the backend's values | `0` |
| `dual_residual`, `max_constraint_violation`, `max_residual_bound_violation` | never computed: always `0` | `0` |
| `total_slack` | sum of the path slacks when the path constraint is softened, else `0` | `0` |
| `used_ekf_fallback` | `false` | `true` |

### is_initialized

```cpp
bool is_initialized() const;
```

Returns `true` once `N + 1` measurements have been accepted since construction or the last refusal, that is, once
updates solve the window.

## Usage Example

<!-- snippet: exceptions -->
```cpp
// gnuplot: plot "< ./nmhe_estimation" using 1:4 with lines title "estimate"
#include <ctrlpp/nmhe.h>
#include <ctrlpp/mpc/nlopt_solver.h>

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <random>
#include <utility>

struct nonlinear_dynamics
{
    double dt{0.1};

    Eigen::Vector2d operator()(const Eigen::Vector2d& x,
                               const Eigen::Matrix<double, 1, 1>& u) const
    {
        return {x[0] + x[1] * dt,
                x[1] + (u[0] - 0.1 * x[1] * std::abs(x[1])) * dt};
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
    constexpr std::size_t N = 8;

    ctrlpp::nmhe_config<double, NX, NU, NY, N> cfg{
        .Q = Eigen::Matrix2d::Identity() * 0.1,
        .R = Eigen::Matrix<double, 1, 1>::Identity() * 0.5,
        .P0 = Eigen::Matrix2d::Identity() * 5.0};

    nonlinear_dynamics dynamics;
    position_measurement measurement;

    using estimator_type =
        ctrlpp::nmhe<double, NX, NU, NY, N,
                     ctrlpp::nlopt_solver<double>,
                     nonlinear_dynamics,
                     position_measurement>;
    auto estimator_result =
        estimator_type::create(dynamics, measurement, cfg);
    if(!estimator_result)
    {
        std::cerr << "invalid NMHE configuration\n";
        return 1;
    }
    auto estimator = std::move(*estimator_result);

    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 0.5);

    Eigen::Vector2d x_true(1.0, 0.5);
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < 40; ++k)
    {
        x_true = dynamics(x_true, u);
        estimator.predict(u);

        Eigen::Matrix<double, 1, 1> z;
        z << x_true[0] + noise(rng);
        if(!estimator.update(z))
            continue; // predicted prior retained; horizon refill begins

        auto x_hat = estimator.state();
        std::cout << "k=" << k
                  << "  true=[" << x_true.transpose() << "]"
                  << "  est=[" << x_hat.transpose() << "]\n";
    }
}
```

## Limitations

### create does not report a solver setup failure

The backend's `setup` runs once, inside the estimator's constructor, and its result is only recorded. `create` then
succeeds, and every window solve falls back to the EKF with `solve_status::error`, indistinguishable from a solve
that failed. Measured with a backend whose `setup` always returns an error: `create` succeeded, and every update from
the `N + 1`-th on reported `used_ekf_fallback` with `solve_status::error`.

### The solver cannot be configured

`nmhe` default-constructs its `Solver` and offers no way to pass a configured instance, so the backend runs at its
default settings: for `nlopt_solver<double>`, SLSQP with its default tolerances and evaluation budget. Other
settings take a wrapper type, named as `Solver`, whose default constructor builds the backend with them and whose
`setup` and `solve` forward to it.

### Only the path constraint is softened

`soft_constraints` adds slacks to the path constraint rows only. The state bounds are variable bounds and the
residual bound is a pair of rows without slack; both stay hard. `soft_penalty` is validated by `create` and
otherwise unused.

### The NLP derivatives are always finite differences

The cost gradient is a central difference of the whole cost, two cost evaluations per decision variable, and the
constraint Jacobian is left to the backend's finite differences. An analytically differentiable `Dynamics` or
`Measurement` speeds up only the embedded EKF; the window NLP never calls its Jacobians.

## See Also

- [mhe](mhe.md)<br/> linear moving horizon estimation
- [nlopt-solver](nlopt-solver.md)<br/> NLopt NLP solver backend
- [argmin-solver](argmin-solver.md)<br/> argmin NLP solver backend
- [ekf](../estimation/ekf.md)<br/> the embedded extended Kalman filter
- [dynamics-model](../model/dynamics-model.md)<br/> dynamics model concept
- [measurement-model](../model/measurement-model.md)<br/> measurement model concept
- [background/nmhe](../../background/nmhe.md)<br/> NMHE theory and background
- [background/mhe](../../background/mhe.md)<br/> MHE theory and background
