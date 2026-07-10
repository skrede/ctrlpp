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

Constructs the controller from a discrete-time state-space model and configuration. Builds the QP matrices, computes terminal cost (via DARE if `Qf` is not set), and initializes the solver.

## Failure contract

Both `mpc` and `nmpc` share one soft-constraint and failure contract. Every
`solve` overload returns `ctrlpp::expected<solve_output<Scalar, NU>, solver_error>`:

- On the **success branch** the result holds a `solve_output` whose `input` field
  is the first control input to apply and whose `status` field is a soft
  `solve_result_status`:
  - `converged`: the solver met its convergence tolerances.
  - `solved_inaccurate`: a usable iterate was returned but tolerances were only
    partially met.
  - `budget_exhausted`: the iteration or time budget ran out, but the best iterate
    found is still returned so it can be commanded knowingly.
- On the **error branch** the result holds a `solver_error` and no input:
  - `infeasible`: the problem as posed has no feasible point.
  - `invalid_problem`: the problem is unbounded, non-convex, the solver reported an
    internal error, or a reference span was too short.
  - `setup_incomplete`: the one-time solver setup failed, so no solve can run.

The input is reached explicitly through `->input`; there is no implicit conversion
to `Vector<Scalar, NU>`, so the soft status can never be silently dropped. On the
error branch the controller does not update its internal previous-input record and
applies no hidden fallback, so a failed solve never warms the rate constraint from
a phantom input. Use `set_applied_input` to record the input the caller actually
commanded.

## Methods

### solve (regulation)

```cpp
[[nodiscard]] auto solve(const Vector<Scalar, NX>& x0)
    -> ctrlpp::expected<solve_output<Scalar, NU>, solver_error>;
```

Solves the QP for regulating state to the origin. Returns the success branch with
the first input and a soft status, or the error branch on failure. If solver setup
failed at construction (reported through the solver's `try_setup`), every `solve`
overload returns the error branch with `solver_error::setup_incomplete`.

### solve (constant reference)

```cpp
[[nodiscard]] auto solve(const Vector<Scalar, NX>& x0,
                         const Vector<Scalar, NY>& y_ref)
    -> ctrlpp::expected<solve_output<Scalar, NU>, solver_error>;
```

Solves the QP for tracking a constant output reference across the entire horizon.

### solve (trajectory reference)

```cpp
[[nodiscard]] auto solve(const Vector<Scalar, NX>& x0,
                         std::span<const Vector<Scalar, NY>> y_ref)
    -> ctrlpp::expected<solve_output<Scalar, NU>, solver_error>;
```

Solves the QP for tracking a time-varying output reference trajectory. The span
must contain at least `horizon + 1` elements; an undersized span is rejected via
the error branch (`solver_error::invalid_problem`) rather than read out of bounds.

### set_applied_input

```cpp
void set_applied_input(const Vector<Scalar, NU>& u);
```

Records the control input the caller actually commanded. The internal previous
input (which anchors the rate constraint) is updated only on a successful solve;
this accessor lets the caller keep that reference consistent when it commands a
different input, for example after an error branch or external saturation.

### trajectory

```cpp
[[nodiscard]] auto trajectory() const
    -> ctrlpp::expected<std::pair<std::vector<Vector<Scalar, NX>>,
                                  std::vector<Vector<Scalar, NU>>>,
                        solver_error>;
```

Returns the full predicted state and input trajectories from the last solve. It is
guarded: before the first valid solve it returns the error branch
(`solver_error::setup_incomplete`) rather than stale or default-initialized data.

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
            // u_opt.error() carries the solver_error describing the failure.
            std::cerr << "MPC solve failed at t=" << t << "\n";
            return EXIT_FAILURE;
        }

        auto diag = controller.diagnostics();
        std::cout << "t=" << t << "  x=[" << x.transpose()
                  << "]  u=" << u_opt->input[0]
                  << "  cost=" << diag.cost << "\n";

        x = ctrlpp::propagate(sys, x, u_opt->input);
    }
}
```

## See Also

- [nmpc](nmpc.md)<br/> nonlinear MPC
- [osqp-solver](osqp-solver.md)<br/> OSQP QP solver backend
- [mhe](mhe.md)<br/> linear moving horizon estimation
- [guides/mpc/solver-injection](../../guides/mpc/solver-injection.md)<br/> solver injection guide
- [guides/intro/your-first-mpc](../../guides/intro/your-first-mpc.md)<br/> introductory MPC tutorial
- [background/mpc](../../background/mpc.md)<br/> MPC theory and background
