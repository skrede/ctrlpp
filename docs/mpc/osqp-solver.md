# osqp_solver

> **Requires OSQP.** Enable with `-DCTRLPP_BUILD_OSQP=ON` when configuring CMake.

OSQP-based QP solver backend for linear MPC and MHE. Satisfies the `qp_solver` concept and can be injected as the `Solver` template parameter of `mpc` and `mhe`. Wraps the OSQP C library with RAII lifetime management, automatic warm-starting, and configurable tolerances.

## Header

| Form | Header |
|------|--------|
| `osqp_solver` | `#include <ctrlpp/mpc/osqp_solver.h>` |

```cpp
class osqp_solver;
```

## qp_solver Concept

`osqp_solver` satisfies the `qp_solver` concept defined in `<ctrlpp/mpc/qp_solver.h>`:

```cpp
template <typename S>
concept qp_solver = requires { typename S::scalar_type; }
    && requires(S solver,
                const qp_problem<typename S::scalar_type>& prob,
                const qp_update<typename S::scalar_type>& upd) {
        { solver.setup(prob) } -> std::same_as<void>;
        { solver.solve(upd) } -> std::same_as<qp_result<typename S::scalar_type>>;
    };
```

Any type satisfying this concept can replace `osqp_solver` as the solver backend.

## Constructors

```cpp
explicit osqp_solver(double eps_abs = 1e-3,
                     double eps_rel = 1e-3,
                     int max_iter = 4000,
                     bool verbose = false,
                     bool warm_starting = true,
                     bool polishing = true);
```

Constructs the solver with OSQP settings. Default values provide a good balance of speed and accuracy for typical MPC/MHE problems.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `eps_abs` | `1e-3` | Absolute convergence tolerance |
| `eps_rel` | `1e-3` | Relative convergence tolerance |
| `max_iter` | `4000` | Maximum ADMM iterations |
| `verbose` | `false` | Print solver progress |
| `warm_starting` | `true` | Enable primal/dual warm-starting |
| `polishing` | `true` | Enable solution polishing |

## Methods

### setup

```cpp
void setup(const qp_problem<double>& problem);
```

Initialises the OSQP workspace from a QP problem (cost matrices P, q and constraint matrices A, l, u). Throws `std::runtime_error` on allocation or setup failure.

### solve

```cpp
auto solve(const qp_update<double>& update) -> qp_result<double>;
```

Solves the QP with updated cost vector, bounds, and optional warm-start vectors. Returns a `qp_result` containing the primal/dual solution, solver status, objective value, iteration count, solve time, and residuals.

## Type Aliases

```cpp
using scalar_type = double;
```

OSQP operates in double precision only.

## Usage Example

```cpp
#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/propagate.h"

#include <Eigen/Dense>

#include <iostream>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr double dt = 0.1;

    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        .A = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        .B = (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(),
        .C = Eigen::Matrix2d::Identity(),
        .D = Eigen::Matrix<double, 2, 1>::Zero()};

    // Custom solver settings: tighter tolerance, fewer iterations
    ctrlpp::osqp_solver solver(1e-5, 1e-5, 2000);

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 15,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Identity()};

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> controller(sys, cfg);

    Eigen::Vector2d x(3.0, 0.0);

    for(int k = 0; k < 50; ++k)
    {
        auto u_opt = controller.solve(x);
        if(!u_opt)
            break;

        x = ctrlpp::propagate(sys, x, *u_opt);
        std::cout << "k=" << k << "  x=[" << x.transpose() << "]\n";
    }
}
```

## See Also

- [mpc](mpc.md) -- linear MPC using OSQP
- [mhe](mhe.md) -- linear MHE using OSQP
- [nlopt-solver](nlopt-solver.md) -- NLopt NLP solver backend
- [guides/mpc/solver-injection](../guides/mpc/solver-injection.md) -- solver injection guide
