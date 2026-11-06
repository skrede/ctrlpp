# Control

Classical and optimal control types for feedback and feedforward control design.
PID handles the majority of single-loop regulation problems, while LQR provides
optimal state-feedback for multi-variable systems. The Riccati solvers are the
design tools the LQR gain functions build on; pole placement is an alternative
design method.

## Types

- [PID](pid/README.md)<br/> Policy-based PID controller with compile-time feature composition
- [lqr](lqr.md)<br/> Linear-quadratic regulator (infinite, finite, time-varying, integral action; discrete and continuous gains)
- [dare](dare.md)<br/> Discrete algebraic Riccati equation solver (real Schur method on the symplectic matrix)
- [care](dare.md#care)<br/> Continuous algebraic Riccati equation solver (matrix sign function by default, or real Schur on the Hamiltonian)
- [place](place.md)<br/> Single-input pole placement via Ackermann's formula
- [mrac](mrac.md)<br/> Model reference adaptive controller with dead-zone, sigma-modification, and e-modification robustification
- [l1](l1.md)<br/> L1 adaptive controller with state predictor, projection-based adaptation, and low-pass filtered control

`<ctrlpp/control.h>` includes every public control header, directly or through
the headers it includes. `<ctrlpp/control/solver_policy.h>` declares the `SolverPolicy`
concept, satisfied by a type with a nested `solver_tag` type; no function or class
in ctrlpp constrains on it. The Riccati solvers take a conditioning policy instead
(see [dare](dare.md#conditioning-policy)).

## When to Use

Pick **PID** when you have a single-loop SISO or MIMO regulation problem and want
to compose exactly the features you need (anti-windup, derivative filtering,
feed-forward) without paying for what you don't.

Pick **LQR** when you have a state-space model and want optimal full-state feedback
that minimizes a quadratic cost. The discrete gain functions call `dare`
internally; `lqr_gain_continuous` solves the continuous equation itself.

Pick **DARE** or **CARE** directly when you need the solution of the discrete or
continuous algebraic Riccati equation outside of LQR (e.g., for terminal cost
computation in MPC).

Pick **place** when you need direct pole assignment of a single-input system rather
than cost-based tuning.

Pick **MRAC** when you have an unknown or uncertain plant and want to track a
reference model adaptively. MRAC adapts online without identifying the plant,
given the sign structure of its input gain; with several inputs it also assumes the
reference model's input matrix is the plant's up to a positive scalar (see
[mrac](mrac.md#limitations)). Robustification policies are selected at compile
time.

Pick **L1** when you want the L1 architecture's separation of adaptation speed from
robustness: its low-pass filter, not the adaptation gain, sets the control
bandwidth. ctrlpp does not check the condition on the filter under which the L1
theory bounds the transient (see [l1](l1.md#limitations)).

## Theory

- [PID Theory](../../background/pid.md)<br/> PID mathematical background
- [LQR Theory](../../background/lqr.md)<br/> LQR and Riccati equation background
- [Adaptive Control Theory](../../background/adaptive-control.md)<br/> MRAC and L1 mathematical background
