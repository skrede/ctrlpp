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

Constructs the estimator from dynamics and measurement models plus configuration. Builds the NLP formulation and initialises the internal EKF for arrival cost propagation.

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

Incorporates a new measurement. During fill-up, delegates to the internal EKF. Once the window is full, solves the NMHE NLP.

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

Returns solver diagnostics including constraint violation metrics and EKF fallback status.

### is_initialized

```cpp
bool is_initialized() const;
```

Returns `true` once the estimation window is full.

## Usage Example

```cpp
#include "ctrlpp/nmhe.h"
#include "ctrlpp/mpc/nlopt_solver.h"

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
        estimator.update(z);

        x_true = dynamics(x_true, u);

        auto x_hat = estimator.state();
        std::cout << "k=" << k
                  << "  true=[" << x_true.transpose() << "]"
                  << "  est=[" << x_hat.transpose() << "]\n";
    }
}
```

## See Also

- [mhe](mhe.md) -- linear moving horizon estimation
- [nlopt-solver](nlopt-solver.md) -- NLopt NLP solver backend
- [dynamics-model](../model/dynamics-model.md) -- dynamics model concept
- [measurement-model](../model/measurement-model.md) -- measurement model concept
- [reference/mhe-theory](../reference/mhe-theory.md) -- MHE theory and background
