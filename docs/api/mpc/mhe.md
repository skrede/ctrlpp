# mhe

> **Requires OSQP.** Enable with `-DCTRLPP_BUILD_OSQP=ON` when configuring CMake.

Linear Moving Horizon Estimator using sparse QP optimization. Performs constrained state estimation over a sliding window of measurements, combining an arrival cost (propagated via an internal EKF) with process and measurement noise weighting. Supports box state constraints with optional soft relaxation and measurement residual bounds. Falls back to the internal EKF during the fill-up phase (fewer than N measurements received).

## Header and Alias

| Form | Header |
|------|--------|
| `mhe<Scalar, NX, NU, NY, N, Solver, Dynamics, Measurement>` | `#include <ctrlpp/mhe.h>` |

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
| `Scalar` | floating-point | Numeric type |
| `NX` | `>= 1` | State dimension |
| `NU` | `>= 1` | Input dimension |
| `NY` | `>= 1` | Measurement dimension |
| `N` | `>= 1` | Horizon length (number of measurement steps in window) |
| `Solver` | satisfies `qp_solver` | QP solver backend. See [osqp-solver](osqp-solver.md). |
| `Dynamics` | satisfies `dynamics_model<Scalar, NX, NU>` | Discrete-time dynamics callable |
| `Measurement` | satisfies `measurement_model<Scalar, NX, NY>` | Measurement callable |

## mhe_config

Configuration struct `mhe_config<Scalar, NX, NU, NY, N>` passed to `create`.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | Process noise covariance |
| `R` | `Matrix<Scalar, NY, NY>` | Identity | Measurement noise covariance |
| `x0` | `Vector<Scalar, NX>` | Zero | Initial state estimate |
| `P0` | `Matrix<Scalar, NX, NX>` | Identity | Initial covariance |
| `arrival_cost_weight` | `Scalar` | `1` | Weight on the arrival cost term |
| `x_min` | `optional<Vector<Scalar, NX>>` | none | Element-wise lower state bound |
| `x_max` | `optional<Vector<Scalar, NX>>` | none | Element-wise upper state bound |
| `residual_bound` | `optional<Vector<Scalar, NY>>` | none | Measurement residual bound |
| `soft_constraints` | `bool` | `true` | Soften state constraints with slack variables |
| `soft_penalty` | `Scalar` | `1e4` | L1 penalty for soft state constraints |
| `numerical_eps` | `Scalar` | `sqrt(eps)` | Perturbation for numerical Jacobians |

## Creation

```cpp
static auto create(
    Dynamics dynamics,
    Measurement measurement,
    const mhe_config<Scalar, NX, NU, NY, N>& config)
    -> ctrlpp::expected<mhe, moving_horizon_construction_error>;
```

Validates the embedded EKF and the moving-horizon formulation before constructing the estimator. On success, initializes the internal EKF for arrival-cost propagation and the measurement/input windows.

`moving_horizon_construction_error` is a variant. Its `filter_error` alternative reports invalid `Q`, `R`, `P0`, `x0`, or numerical-differentiation configuration for the embedded EKF. Its `moving_horizon_configuration_error` alternative reports non-invertible covariance weights, nonpositive arrival or soft-constraint weights, invalid state or residual bounds, and an invalid numerical-differentiation step.

The window length `N` is a template parameter on both the class and `mhe_config`, not a runtime field, so `N == 0` still fails at compile time. Runtime-valued configuration is checked by `create`; callers must branch on its `expected` before moving out the estimator.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates the internal EKF one step and records the input in the sliding window.

### update

```cpp
auto update(const output_vector_t& z)
    -> ctrlpp::expected<void, ekf_update_error>;
```

Incorporates a new measurement, or reports why it could not. During fill-up (fewer than N steps), delegates to the internal EKF. Once the window is full, solves the estimation problem to refine the state trajectory over the entire window. A solver setup failure (reported through the solver's `setup`) or a non-optimal solve falls back to the internal EKF; `diagnostics()` reports the fallback.

**Two channels, and which one carries what.** A measurement the embedded filter cannot use is a FAILURE: the step did not happen, and `update` returns the filter's own `ekf_update_error` verbatim rather than restating the same three conditions under a second name that could drift from the one that decides. Everything else is a DISPOSITION on a step that succeeded -- which of the two estimators produced the estimate, and how the solve went -- and that is what `diagnostics()` carries.

A refusal cannot undo the preceding prediction: that input was applied to the plant, so `state()` exposes the embedded filter's predicted prior. The missing measurement does invalidate the fixed-step horizon. The estimator clears its input, measurement, prior, and warm-start windows, sets `is_initialized()` to false, and uses the embedded EKF for the next N accepted measurements. Optimization resumes only after a coherent horizon has been refilled. `diagnostics()` continues to describe the last successful update until the next accepted measurement; the failed `update` return is what identifies the current state as an uncorrected prior.

An invalid solver result falls back the same way. A solver may report an accepted status and still return a primal shorter than the decision dimension, a dual shorter than the constraint count, a non-finite consumed value or diagnostic, or a negative iteration count. The estimator validates all of these before reading or storing the result, and engages the EKF fallback instead of writing the window when validation fails. The estimate is still produced -- by the fallback rather than by the window solve -- so this is a disposition and not a failure: `diagnostics().used_ekf_fallback` is `true` and `diagnostics().status` is `solve_status::invalid_backend_result`, which names this condition specifically rather than collapsing it into the general `solve_status::error` a non-optimal solve reports.

### state

```cpp
const state_vector_t& state() const;
```

Returns the current state estimate (last element of the window).

### covariance

```cpp
const cov_matrix_t& covariance() const;
```

Returns the current covariance estimate from the internal EKF.

### innovation

```cpp
const output_vector_t& innovation() const;
```

Returns the measurement innovation (residual) from the last update.

### trajectory

```cpp
std::span<const state_vector_t> trajectory() const;
```

Returns the smoothed state trajectory over the full estimation window (N+1 elements).

### diagnostics

```cpp
const mhe_diagnostics<Scalar>& diagnostics() const;
```

Returns solver diagnostics including status, cost, residuals, slack usage, and whether the EKF fallback was used. This is the DISPOSITION channel and it describes a step that succeeded: `solve_status::error` with `used_ekf_fallback` set for a setup failure or a non-optimal solve, and `solve_status::invalid_backend_result` for an incomplete or non-finite accepted solver result. A refused measurement never appears here, because no estimate was produced for the aggregate to describe -- it is returned by `update` instead, and the aggregate goes on describing the last accepted step.

### is_initialized

```cpp
bool is_initialized() const;
```

Returns `true` once the estimation window is full (at least N measurement steps received).

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
        .x0 = Eigen::Vector2d::Zero(),
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

## See Also

- [nmhe](nmhe.md)<br/> nonlinear moving horizon estimation
- [osqp-solver](osqp-solver.md)<br/> OSQP QP solver backend
- [dynamics-model](../model/dynamics-model.md)<br/> dynamics model concept
- [measurement-model](../model/measurement-model.md)<br/> measurement model concept
- [background/mhe](../../background/mhe.md)<br/> MHE theory and background
