# MPC and MHE

Model predictive control (MPC) and moving horizon estimation (MHE) types. The linear
variants solve quadratic programs through OSQP or argmin's header-only sparse ADMM QP
backend; the nonlinear variants solve nonlinear programs through NLopt or argmin's
header-only NLP backend. Solver backends are injected as a template parameter checked
by C++20 concepts, so you can supply your own QP or NLP solver without changing
controller code.

## Types

### Controllers

- [mpc](mpc.md)<br/> Linear model predictive control (sparse QP; OSQP or argmin backend)
- [nmpc](nmpc.md)<br/> Nonlinear model predictive control (multiple shooting; argmin or NLopt backend)

### Estimators

- [mhe](mhe.md)<br/> Linear moving horizon estimation (QP; OSQP or argmin backend)
- [nmhe](nmhe.md)<br/> Nonlinear moving horizon estimation (NLP; NLopt or argmin backend)

### Solver Backends

- [osqp_solver](osqp-solver.md)<br/> OSQP quadratic program solver wrapper (the OSQP C library)
- [argmin_qp_solver](argmin-qp-solver.md)<br/> argmin sparse ADMM QP backend (header-only; same concept as `osqp_solver`, different setup error)
- [nlopt_solver](nlopt-solver.md)<br/> NLopt nonlinear program solver wrapper (runtime-horizon `nmpc_dynamic` and `nmhe`)
- [argmin_solver](argmin-solver.md)<br/> argmin NLP backend (header-only; the only shipped backend for the compile-time-horizon `nmpc`)

## When to use

Pick **mpc** for linear systems where you need constraint handling (input bounds,
state bounds, terminal constraints). It formulates a sparse QP solved by the injected
QP backend.

Pick **nmpc** for nonlinear systems. It uses direct multiple shooting with continuity
constraints, solved by the injected NLP backend. Supports path constraints `g(x,u) <= 0`
and terminal constraints `h(x_N) <= 0`, softened by default with an L1-penalized slack.

Pick **mhe** or **nmhe** when you need constrained state estimation &mdash; the MHE
dual of MPC. Useful when states have physical bounds (e.g., non-negative
concentrations).

## Theory

For the mathematical background, see
[MPC Theory](../../background/mpc.md),
[NMPC Theory](../../background/nmpc.md) and
[MHE Theory](../../background/mhe.md).
