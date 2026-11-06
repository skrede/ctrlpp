# Composition Patterns in ctrlpp

ctrlpp uses three composition mechanisms to build complex controllers from
simple building blocks. Each mechanism solves a different problem, and they
can be combined freely.

## 1. Policy Composition (PID)

Template parameter packs add behaviors at compile time. A policy that has
settings adds its own configuration struct to the controller's `pid_config`,
reached through `cfg.policy<P>()`. You compose exactly the features you need;
unused policies cost nothing.

```cpp
using Pid = ctrlpp::pid<double, 1,
                         ctrlpp::anti_windup<ctrlpp::back_calc>,
                         ctrlpp::deriv_filter,
                         ctrlpp::rate_limit>;
```

The compiler resolves all policies at instantiation, with no virtual dispatch and
no runtime branching on which policies are present.

**Details:** [PID Composition Guide](../pid/composition.md)

## 2. Concept-Based Injection (Solver Backends)

C++20 concepts define what a solver backend must provide. `mpc` constrains its
solver parameter with `qp_solver`, and `nmpc_dynamic` constrains its own with
`nlp_solver`:

```cpp
auto linear_mpc = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(sys, mpc_cfg);
auto nonlinear_mpc = ctrlpp::nmpc_dynamic<double, NX, NU, ctrlpp::nlopt_solver<double>,
                                          Dynamics>::create(dyn, nmpc_cfg);
```

Swap `osqp_solver` for your own type and MPC works identically. The concept
check catches mismatches at compile time. Regardless of the backend, both `mpc`
and `nmpc` expose one shared failure contract: `solve()` returns
`ctrlpp::expected<solve_output<Scalar, NU>, solver_error>`, so a backend that
fails surfaces the error on the same channel the caller already handles.

**Details:** [Solver Injection Guide](../mpc/solver-injection.md)

## 3. Observer-Controller Composition

Observers and controllers are separate objects connected at the application
level. Every ctrlpp estimator satisfies the `ObserverPolicy` concept, which
requires `predict(u)`, `update(z)` and `state()`, so the loop below is written
the same way whichever estimator it holds. No controller takes an observer as a
parameter; the application passes the estimate on:

```cpp
kf.predict(u);
if(!kf.update(z))
{
    std::cerr << "Kalman filter rejected the measurement\n";
    return EXIT_FAILURE;
}
auto u_next = pid.compute(setpoint, kf.state(), dt);
if(!u_next)
{
    std::cerr << "PID rejected the cycle\n";
    return EXIT_FAILURE;
}
```

This works with `kalman_filter`, `ekf`, `ukf`, `mekf`, `manifold_ukf`,
`luenberger_observer`, `complementary_filter`, `particle_filter`, `mhe`, `nmhe`,
or any custom type satisfying `ObserverPolicy`.

**Details:** [Observer-Controller Guide](../estimation/observer-controller.md)

## Combining All Three

A real system might use all three mechanisms simultaneously:

- **PID with policies** for inner-loop velocity control
- **MPC with OSQP solver** for outer-loop trajectory tracking
- **EKF observer** providing state estimates to both controllers

Each mechanism is orthogonal: policies configure a single controller, solver
injection selects an optimization backend, and observer-controller composition
connects estimators to controllers. They compose without interference.

## Next Steps

- [PID Composition](../pid/composition.md)<br/> policy details and examples
- [Solver Injection](../mpc/solver-injection.md)<br/> concept definitions and
  custom solvers
- [Observer-Controller](../estimation/observer-controller.md)<br/> estimation
  and control loop integration
