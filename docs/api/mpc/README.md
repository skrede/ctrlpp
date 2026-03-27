# MPC and MHE

Model predictive control (MPC) and moving horizon estimation (MHE) types. Linear
variants solve quadratic programs via OSQP; nonlinear variants solve nonlinear
programs via NLopt. Solver backends are injected via C++23 concepts -- you can
swap in your own QP or NLP solver without changing controller code.

## Types

### Controllers

- [mpc](mpc.md) -- Linear model predictive control (sparse QP, OSQP backend)
- [nmpc](nmpc.md) -- Nonlinear model predictive control (multiple shooting, NLopt backend)

### Estimators

- [mhe](mhe.md) -- Linear moving horizon estimation (OSQP backend)
- [nmhe](nmhe.md) -- Nonlinear moving horizon estimation (NLopt backend)

### Solver Backends

- [osqp_solver](osqp-solver.md) -- OSQP quadratic program solver wrapper
- [nlopt_solver](nlopt-solver.md) -- NLopt nonlinear program solver wrapper

## When to use

Pick **mpc** for linear systems where you need constraint handling (input bounds,
state bounds, terminal constraints). It formulates a sparse QP solved by OSQP.

Pick **nmpc** for nonlinear systems. It uses multiple shooting with continuity
constraints, solved by NLopt. Supports path constraints `g(x,u) <= 0` and terminal
constraints `h(x_N) <= 0` with soft L1 penalty.

Pick **mhe** or **nmhe** when you need constrained state estimation -- the MHE
dual of MPC. Useful when states have physical bounds (e.g., non-negative
concentrations).

## Theory

For the mathematical background, see
[MPC Theory](../reference/mpc-theory.md) and
[MHE Theory](../reference/mhe-theory.md).
