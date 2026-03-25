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
| `Qf` | `optional<Matrix<Scalar, NX, NX>>` | DARE solution | Terminal cost weight. Computed from DARE when not provided. |
| `u_min` | `optional<Vector<Scalar, NU>>` | none | Element-wise lower input bound |
| `u_max` | `optional<Vector<Scalar, NU>>` | none | Element-wise upper input bound |
| `x_min` | `optional<Vector<Scalar, NX>>` | none | Element-wise lower state bound |
| `x_max` | `optional<Vector<Scalar, NX>>` | none | Element-wise upper state bound |
| `du_max` | `optional<Vector<Scalar, NU>>` | none | Maximum input rate of change per step |
| `soft_penalty` | `Scalar` | `1e4` | L1 penalty for soft state constraints |
| `soft_state_penalty` | `optional<Vector<Scalar, NX>>` | none | Per-state soft constraint penalty |
| `terminal_constraint_set` | `optional<terminal_set<Scalar, NX>>` | none | Ellipsoidal or polyhedral terminal set |
| `hard_state_constraints` | `bool` | `false` | When true, state constraints are hard (not softened) |

## Constructors

```cpp
mpc(const discrete_state_space<Scalar, NX, NU, NX>& system,
    const mpc_config<Scalar, NX, NU>& config);
```

Constructs the controller from a discrete-time state-space model and configuration. Builds the QP matrices, computes terminal cost (via DARE if `Qf` is not set), and initialises the solver.

## Methods

### solve (regulation)

```cpp
[[nodiscard]] auto solve(const Vector<Scalar, NX>& x0)
    -> std::optional<Vector<Scalar, NU>>;
```

Solves the QP for regulating state to the origin. Returns the first optimal input or `std::nullopt` if the solver fails.

### solve (constant reference)

```cpp
[[nodiscard]] auto solve(const Vector<Scalar, NX>& x0,
                         const Vector<Scalar, NX>& x_ref)
    -> std::optional<Vector<Scalar, NU>>;
```

Solves the QP for tracking a constant reference across the entire horizon.

### solve (trajectory reference)

```cpp
[[nodiscard]] auto solve(const Vector<Scalar, NX>& x0,
                         std::span<const Vector<Scalar, NX>> x_ref)
    -> std::optional<Vector<Scalar, NU>>;
```

Solves the QP for tracking a time-varying reference trajectory. The span must contain at least `horizon + 1` elements.

### trajectory

```cpp
[[nodiscard]] auto trajectory() const
    -> std::pair<std::vector<Vector<Scalar, NX>>,
                 std::vector<Vector<Scalar, NU>>>;
```

Returns the full predicted state and input trajectories from the last solve.

### diagnostics

```cpp
[[nodiscard]] auto diagnostics() const -> mpc_diagnostics<Scalar>;
```

Returns solver diagnostics from the last solve, including status, iteration count, solve time, cost, and residuals.

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

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> controller(sys, cfg);

    Eigen::Vector2d x(5.0, 0.0);

    for(double t = 0.0; t < 5.0; t += dt)
    {
        auto u_opt = controller.solve(x);
        if(!u_opt)
        {
            std::cerr << "MPC solve failed at t=" << t << "\n";
            return EXIT_FAILURE;
        }

        auto diag = controller.diagnostics();
        std::cout << "t=" << t << "  x=[" << x.transpose()
                  << "]  u=" << (*u_opt)[0]
                  << "  cost=" << diag.cost << "\n";

        x = ctrlpp::propagate(sys, x, *u_opt);
    }
}
```

## See Also

- [nmpc](nmpc.md) -- nonlinear MPC
- [osqp-solver](osqp-solver.md) -- OSQP QP solver backend
- [mhe](mhe.md) -- linear moving horizon estimation
- [guides/mpc/solver-injection](../guides/mpc/solver-injection.md) -- solver injection guide
- [guides/intro/your-first-mpc](../guides/intro/your-first-mpc.md) -- introductory MPC tutorial
- [reference/mpc-theory](../reference/mpc-theory.md) -- MPC theory and background
