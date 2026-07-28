# nmhe

> **Requires NLopt.** Enable with `-DCTRLPP_BUILD_NLOPT=ON` when configuring CMake.

Nonlinear Moving Horizon Estimator using NLP optimization with multiple shooting. Handles nonlinear dynamics and measurement models for constrained state estimation over a sliding window. Supports general nonlinear path constraints with soft relaxation and propagates arrival cost via an internal EKF. Falls back to the EKF during the fill-up phase.

## Header and Alias

| Form | Header |
|------|--------|
| `nmhe<Scalar, NX, NU, NY, N, Solver, Dynamics, Measurement, NC>` | `#include <ctrlpp/nmhe.h>` |

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
| `Scalar` | floating-point | Numeric type |
| `NX` | `>= 1` | State dimension |
| `NU` | `>= 1` | Input dimension |
| `NY` | `>= 1` | Measurement dimension |
| `N` | `>= 1` | Horizon length |
| `Solver` | satisfies `nlp_solver` | NLP solver backend. See [nlopt-solver](nlopt-solver.md). |
| `Dynamics` | satisfies `dynamics_model<Scalar, NX, NU>` | Discrete-time dynamics callable |
| `Measurement` | satisfies `measurement_model<Scalar, NX, NY>` | Measurement callable |
| `NC` | `>= 0` | Number of path constraints (default 0) |

## nmhe_config

Configuration struct `nmhe_config<Scalar, NX, NU, NY, N, NC>` passed at construction.

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
| `soft_constraints` | `bool` | `true` | Soften constraints with slack variables |
| `soft_penalty` | `Scalar` | `1e4` | L1 penalty for soft constraints |
| `numerical_eps` | `Scalar` | `sqrt(eps)` | Perturbation for numerical Jacobians |
| `path_constraint` | `optional<function<Vector<NC>(x)>>` | none | Nonlinear path constraint g(x) <= 0 |
| `path_penalty` | `Vector<Scalar, NC>` | `1e4` each | Per-constraint L1 penalty |

## Constructors

```cpp
nmhe(Dynamics dynamics, Measurement measurement,
     const nmhe_config<Scalar, NX, NU, NY, N, NC>& config);
```

Constructs the estimator from dynamics and measurement models plus configuration. Builds the NLP formulation and initializes the internal EKF for arrival cost propagation.

The window length `N` is a template parameter on both the class and `nmhe_config`, not a runtime field, so its domain is enforced at compile time: `N == 0` fails to compile on both. `N` sizes the fixed estimation window arrays, which the update rotates and reads the trailing element of, neither of which is defined for an empty window. There is therefore no fallible construction factory here: the horizon domain is closed before the program runs.

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

A refusal is a complete no-op: no window shifts, no counter increments, and the input the last `predict` supplied is not committed to the horizon either. `state()` therefore still returns the estimate the last ACCEPTED measurement produced, which is a perfectly plausible number that nothing about its value distinguishes from a fresh one -- and that is precisely why the refusal is returned instead of being encoded in a status flag. `diagnostics()` likewise still describes that last accepted step, so the two accessors always describe the same step.

An ill-shaped solver result falls back the same way. A solver may report an accepted status and still return a decision vector shorter than the NLP the estimator posed; the estimator compares the reported length against that dimension before the extraction reads the result, and on a violation engages the EKF fallback instead of writing the window. The estimate is still produced -- by the fallback rather than by the window solve -- so this is a disposition and not a failure: `diagnostics().used_ekf_fallback` is `true` and `diagnostics().status` is `solve_status::invalid_backend_result`, which names this condition specifically rather than collapsing it into the general `solve_status::error` a non-optimal solve reports.

### state

```cpp
const state_vector_t& state() const;
```

Returns the current state estimate (last element of the window).

### covariance

```cpp
const cov_matrix_t& covariance() const;
```

Returns the current covariance from the internal EKF.

### trajectory

```cpp
std::span<const state_vector_t> trajectory() const;
```

Returns the smoothed state trajectory over the estimation window (N+1 elements).

### diagnostics

```cpp
const mhe_diagnostics<Scalar>& diagnostics() const;
```

Returns solver diagnostics including constraint violation metrics and EKF fallback status. This is the DISPOSITION channel and it describes a step that succeeded: `solve_status::error` with `used_ekf_fallback` set for a setup failure or a non-optimal solve, and `solve_status::invalid_backend_result` for a solver result whose dimensions did not cover the posed problem. A refused measurement never appears here, because no estimate was produced for the aggregate to describe -- it is returned by `update` instead, and the aggregate goes on describing the last accepted step.

### is_initialized

```cpp
bool is_initialized() const;
```

Returns `true` once the estimation window is full.

## Usage Example

```cpp
// gnuplot: plot "< ./nmhe_estimation" using 1:4 with lines title "estimate"
#include <ctrlpp/nmhe.h>
#include <ctrlpp/mpc/nlopt_solver.h>

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <random>

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
        .x0 = Eigen::Vector2d::Zero(),
        .P0 = Eigen::Matrix2d::Identity() * 5.0};

    nonlinear_dynamics dynamics;
    position_measurement measurement;

    ctrlpp::nmhe<double, NX, NU, NY, N,
                 ctrlpp::nlopt_solver<double>,
                 nonlinear_dynamics,
                 position_measurement>
        estimator(dynamics, measurement, cfg);

    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 0.5);

    Eigen::Vector2d x_true(1.0, 0.5);
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < 40; ++k)
    {
        estimator.predict(u);

        Eigen::Matrix<double, 1, 1> z;
        z << x_true[0] + noise(rng);
        if(!estimator.update(z))
            continue; // the measurement was refused; the estimate is the last accepted one

        x_true = dynamics(x_true, u);

        auto x_hat = estimator.state();
        std::cout << "k=" << k
                  << "  true=[" << x_true.transpose() << "]"
                  << "  est=[" << x_hat.transpose() << "]\n";
    }
}
```

## See Also

- [mhe](mhe.md)<br/> linear moving horizon estimation
- [nlopt-solver](nlopt-solver.md)<br/> NLopt NLP solver backend
- [dynamics-model](../model/dynamics-model.md)<br/> dynamics model concept
- [measurement-model](../model/measurement-model.md)<br/> measurement model concept
- [background/mhe](../../background/mhe.md)<br/> MHE theory and background
