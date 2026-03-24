# Composition Patterns in ctrlpp

ctrlpp uses three composition mechanisms to build complex controllers from
simple building blocks. Each mechanism solves a different problem, and they
can be combined freely.

## 1. Policy Composition (PID)

Template parameter packs add behaviours at compile time. Each policy extends
the config struct with its own fields. You compose exactly the features you
need; unused policies cost nothing.

```cpp
using Pid = ctrlpp::pid<double, 1, 1, 1,
                         ctrlpp::anti_windup<ctrlpp::back_calc>,
                         ctrlpp::deriv_filter,
                         ctrlpp::rate_limit>;
```

The compiler resolves all policies at instantiation -- no virtual dispatch,
no runtime branching.

**Details:** [PID Composition Guide](../pid/composition.md)

## 2. Concept-Based Injection (Solver Backends)

C++23 concepts define what a solver backend must provide. MPC and NMPC accept
any type satisfying `qp_solver` or `nlp_solver` as a template parameter:

```cpp
ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> linear_mpc(sys, cfg);
ctrlpp::nmpc<double, NX, NU, ctrlpp::nlopt_solver> nonlinear_mpc(dyn, cfg);
```

Swap `osqp_solver` for your own type and MPC works identically. The concept
check catches mismatches at compile time.

**Details:** [Solver Injection Guide](../mpc/solver-injection.md)

## 3. Observer-Controller Composition

Observers and controllers are separate objects connected at the application
level. The `ObserverPolicy` concept ensures any observer can feed any
controller:

```cpp
kf.predict(u);
kf.update(z);
auto u_next = controller.compute(kf.state());
```

This works with `kalman_filter`, `ekf`, `ukf`, `particle_filter`, or any
custom type satisfying `ObserverPolicy`.

**Details:** [Observer-Controller Guide](../estimation/observer-controller.md)

## Combining All Three

A real system might use all three mechanisms simultaneously:

- **PID with policies** for inner-loop velocity control
- **MPC with OSQP solver** for outer-loop trajectory tracking
- **EKF observer** providing state estimates to both controllers

Each mechanism is orthogonal: policies configure a single controller, solver
injection selects an optimisation backend, and observer-controller composition
connects estimators to controllers. They compose without interference.

## Next Steps

- [PID Composition](../pid/composition.md) -- policy details and examples
- [Solver Injection](../mpc/solver-injection.md) -- concept definitions and
  custom solvers
- [Observer-Controller](../estimation/observer-controller.md) -- estimation
  and control loop integration
