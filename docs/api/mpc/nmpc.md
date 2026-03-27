# nmpc

> **Requires NLopt.** Enable with `-DCTRLPP_BUILD_NLOPT=ON` when configuring CMake.

Nonlinear Model Predictive Controller using NLP optimization with multiple shooting. Supports arbitrary nonlinear dynamics via the `dynamics_model` concept, custom stage and terminal cost overrides, nonlinear path and terminal constraints with optional soft relaxation, and automatic warm-starting via trajectory shifting.

## Header and Alias

| Form | Header |
|------|--------|
| `nmpc<Scalar, NX, NU, Solver, Dynamics, NC, NTC>` | `#include <ctrlpp/nmpc.h>` |

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU,
          nlp_solver Solver,
          dynamics_model<Scalar, NX, NU> Dynamics,
          std::size_t NC = 0, std::size_t NTC = 0>
class nmpc;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type (`double`, `float`) |
| `NX` | `>= 1` | State dimension |
| `NU` | `>= 1` | Input dimension |
| `Solver` | satisfies `nlp_solver` | NLP solver backend. See [nlopt-solver](nlopt-solver.md). |
| `Dynamics` | satisfies `dynamics_model<Scalar, NX, NU>` | Discrete-time dynamics callable. See [dynamics-model](../model/dynamics-model.md). |
| `NC` | `>= 0` | Number of path constraints (default 0) |
| `NTC` | `>= 0` | Number of terminal constraints (default 0) |

## nmpc_config

Configuration struct `nmpc_config<Scalar, NX, NU, NC, NTC>` passed at construction.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `horizon` | `int` | `1` | Prediction horizon length N |
| `Q` | `Matrix<Scalar, NX, NX>` | Identity | State cost weight |
| `R` | `Matrix<Scalar, NU, NU>` | Identity | Input cost weight |
| `Qf` | `optional<Matrix<Scalar, NX, NX>>` | none | Terminal cost weight |
| `u_min` | `optional<Vector<Scalar, NU>>` | none | Element-wise lower input bound |
| `u_max` | `optional<Vector<Scalar, NU>>` | none | Element-wise upper input bound |
| `x_min` | `optional<Vector<Scalar, NX>>` | none | Element-wise lower state bound |
| `x_max` | `optional<Vector<Scalar, NX>>` | none | Element-wise upper state bound |
| `du_max` | `optional<Vector<Scalar, NU>>` | none | Maximum input rate of change per step |
| `stage_cost` | `optional<function<Scalar(x, u)>>` | none | Custom stage cost override (replaces Q/R quadratic) |
| `terminal_cost` | `optional<function<Scalar(x)>>` | none | Custom terminal cost override (replaces Qf quadratic) |
| `path_constraint` | `optional<function<Vector<NC>(x, u)>>` | none | Path constraint g(x,u) <= 0 element-wise |
| `terminal_constraint` | `optional<function<Vector<NTC>(x)>>` | none | Terminal constraint h(x_N) <= 0 element-wise |
| `soft_constraints` | `bool` | `true` | Soften constraints with slack variables |
| `path_penalty` | `Vector<Scalar, NC>` | `1e4` each | Per-constraint L1 penalty for path constraints |
| `terminal_penalty` | `Vector<Scalar, NTC>` | `1e4` each | Per-constraint L1 penalty for terminal constraints |

## Constructors

```cpp
nmpc(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config);
```

Constructs the controller from a dynamics model and configuration. Builds the NLP formulation and initialises the solver.

## Methods

### solve (regulation)

```cpp
std::optional<Vector<Scalar, NU>> solve(const Vector<Scalar, NX>& x0);
```

Solves the NLP for regulating state to the origin.

### solve (constant reference)

```cpp
std::optional<Vector<Scalar, NU>> solve(const Vector<Scalar, NX>& x0,
                                        const Vector<Scalar, NX>& x_ref);
```

Solves the NLP for tracking a constant reference.

### solve (trajectory reference)

```cpp
std::optional<Vector<Scalar, NU>> solve(const Vector<Scalar, NX>& x0,
                                        std::span<const Vector<Scalar, NX>> x_ref);
```

Solves the NLP for tracking a time-varying reference trajectory.

### trajectory

```cpp
std::pair<std::vector<Vector<Scalar, NX>>,
          std::vector<Vector<Scalar, NU>>> trajectory() const;
```

Returns the full predicted state and input trajectories from the last solve.

### diagnostics

```cpp
mpc_diagnostics<Scalar> diagnostics() const;
```

Returns solver diagnostics including constraint violation metrics and total slack.

## Usage Example

```cpp
// gnuplot: plot "< ./nmpc_pendulum" using 1:2 with lines title "theta"
#include <ctrlpp/nmpc.h>
#include <ctrlpp/mpc/nlopt_solver.h>

#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iostream>

// Simple pendulum dynamics: x = [theta, omega]
struct pendulum_dynamics
{
    double dt{0.05};
    double g{9.81};
    double l{1.0};

    Eigen::Vector2d operator()(const Eigen::Vector2d& x,
                               const Eigen::Matrix<double, 1, 1>& u) const
    {
        double theta = x[0];
        double omega = x[1];
        double alpha = -g / l * std::sin(theta) + u[0];
        return {theta + omega * dt, omega + alpha * dt};
    }
};

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;

    pendulum_dynamics dynamics{.dt = 0.05};

    ctrlpp::nmpc_config<double, NX, NU> cfg{
        .horizon = 30,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity() * 0.01,
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-5.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(5.0)};

    ctrlpp::nmpc<double, NX, NU, ctrlpp::nlopt_solver<double>, pendulum_dynamics>
        controller(dynamics, cfg);

    Eigen::Vector2d x(1.0, 0.0);  // Start at 1 radian
    Eigen::Vector2d x_ref(0.0, 0.0);  // Swing up to vertical

    for(int k = 0; k < 100; ++k)
    {
        auto u_opt = controller.solve(x, x_ref);
        if(!u_opt)
        {
            std::cerr << "NMPC solve failed at step " << k << "\n";
            return EXIT_FAILURE;
        }

        std::cout << "k=" << k
                  << "  theta=" << x[0]
                  << "  omega=" << x[1]
                  << "  u=" << (*u_opt)[0] << "\n";

        x = dynamics(x, *u_opt);
    }
}
```

## See Also

- [mpc](mpc.md) -- linear MPC
- [nlopt-solver](nlopt-solver.md) -- NLopt NLP solver backend
- [dynamics-model](../model/dynamics-model.md) -- dynamics model concept
- [constraint-model](../model/constraint-model.md) -- constraint model concepts
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md) -- solver injection guide
- [background/mpc](../../background/mpc.md) -- MPC theory and background
