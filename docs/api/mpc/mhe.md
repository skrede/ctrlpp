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

Configuration struct `mhe_config<Scalar, NX, NU, NY, N>` passed at construction.

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

## Constructors

```cpp
mhe(Dynamics dynamics, Measurement measurement,
    const mhe_config<Scalar, NX, NU, NY, N>& config);
```

Constructs the estimator from dynamics and measurement models plus configuration. Initialises the internal EKF for arrival cost propagation and the measurement/input windows.

## Methods

### predict

```cpp
void predict(const input_vector_t& u);
```

Propagates the internal EKF one step and records the input in the sliding window.

### update

```cpp
void update(const output_vector_t& z);
```

Incorporates a new measurement. During fill-up (fewer than N steps), delegates to the internal EKF. Once the window is full, solves the MHE QP to refine the state trajectory over the entire window.

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

Returns solver diagnostics including status, cost, residuals, slack usage, and whether the EKF fallback was used.

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

    ctrlpp::mhe<double, NX, NU, NY, N,
                ctrlpp::osqp_solver,
                constant_dynamics,
                position_measurement>
        estimator(constant_dynamics{}, position_measurement{}, cfg);

    std::mt19937 rng(42);
    std::normal_distribution<double> noise(0.0, 1.0);

    Eigen::Vector2d x_true(3.0, -1.0);
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();

    for(int k = 0; k < 30; ++k)
    {
        estimator.predict(u);

        Eigen::Matrix<double, 1, 1> z;
        z << x_true[0] + noise(rng);
        estimator.update(z);

        auto x_hat = estimator.state();
        std::cout << "k=" << k
                  << "  true=[" << x_true.transpose() << "]"
                  << "  est=[" << x_hat.transpose() << "]"
                  << "  init=" << estimator.is_initialized() << "\n";
    }
}
```

## See Also

- [nmhe](nmhe.md) -- nonlinear moving horizon estimation
- [osqp-solver](osqp-solver.md) -- OSQP QP solver backend
- [dynamics-model](../model/dynamics-model.md) -- dynamics model concept
- [measurement-model](../model/measurement-model.md) -- measurement model concept
- [reference/mhe-theory](../reference/mhe-theory.md) -- MHE theory and background
