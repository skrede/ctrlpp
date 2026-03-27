# Your First MPC Controller

> **Requires OSQP.** Enable with `-DCTRLPP_BUILD_OSQP=ON`.

This tutorial builds a model predictive controller (MPC) for a double
integrator -- a system with position and velocity states driven by a force
input. MPC plans an optimal trajectory over a finite horizon while respecting
state and input constraints.

**Prerequisites:** ctrlpp installed per [Getting Started](../../getting-started.md),
OSQP enabled.

## The System

A double integrator discretised with sampling time `dt`:

```
x = [position, velocity]
u = [force]

A = [1  dt ]    B = [0.5*dt^2]
    [0  1  ]        [dt      ]
```

## Complete Program

```cpp
// Usage: ./your_first_mpc | gnuplot -p -e "set datafile separator ','; plot '-' skip 1 using 1:2 with lines title 'position', '' using 1:3 with lines title 'velocity', '' using 1:4 with lines title 'control'"
#include <ctrlpp/mpc.h>
#include <ctrlpp/mpc/osqp_solver.h>
#include <ctrlpp/model/propagate.h>

#include <Eigen/Dense>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr double dt = 0.1;

    // Double integrator state-space
    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        .A = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        .B = (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(),
        .C = Eigen::Matrix2d::Identity(),
        .D = Eigen::Matrix<double, 2, 1>::Zero()};

    // MPC configuration: horizon, cost weights, constraints
    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 20,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity(),
        .Qf = std::nullopt,
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-1.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(1.0),
        .x_min = Eigen::Vector2d(
            -std::numeric_limits<double>::infinity(), -2.0),
        .x_max = Eigen::Vector2d(
            std::numeric_limits<double>::infinity(), 2.0)};

    // Create MPC with OSQP solver backend
    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> controller(sys, cfg);

    // Initial state: position=5, velocity=0 (regulate to origin)
    Eigen::Vector2d x(5.0, 0.0);
    constexpr double duration = 10.0;

    std::cout << "time,position,velocity,control\n";

    for (double t = 0.0; t < duration; t += dt)
    {
        auto u_opt = controller.solve(x);
        if (!u_opt)
        {
            std::cerr << "MPC solve failed at t=" << t << "\n";
            return EXIT_FAILURE;
        }

        Eigen::Matrix<double, 1, 1> u = *u_opt;

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << x[0] << "," << x[1] << ","
                  << u[0] << "\n";

        x = ctrlpp::propagate(sys, x, u);
    }
}
```

## How MPC Works

At every time step:

1. **Solve** -- `controller.solve(x)` formulates a quadratic program (QP) over
   the prediction horizon. It finds the optimal sequence of inputs that
   minimises the cost while satisfying constraints.

2. **Apply first input** -- only the first input from the optimal sequence is
   applied. This is the *receding horizon* principle.

3. **Propagate** -- the system advances one step, and the whole process repeats.

## Key Configuration

| Field     | Purpose                                         |
| --------- | ----------------------------------------------- |
| `horizon` | Number of prediction steps                      |
| `Q`       | State cost weight (penalise deviation from zero)|
| `R`       | Input cost weight (penalise control effort)     |
| `u_min/max` | Input constraints (actuator limits)           |
| `x_min/max` | State constraints (e.g., velocity limits)     |

## Next Steps

- [MPC API Reference](../../api/mpc/mpc.md) -- full interface, diagnostics, warm
  starting
- [Solver Injection Guide](../mpc/solver-injection.md) -- swapping QP/NLP
  solver backends
- [Observer-Controller Composition](../estimation/observer-controller.md) --
  MPC with state estimation
- [MPC Theory](../../background/mpc.md) -- QP formulation, stability,
  terminal constraints
