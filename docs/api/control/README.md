# Control

Classical and optimal control types for feedback and feedforward control design.
PID handles the majority of single-loop regulation problems, while LQR provides
optimal state-feedback for multi-variable systems. DARE and place are the
underlying design tools that LQR builds on.

## Types

- [PID](pid/README.md) -- Policy-based PID controller with compile-time feature composition
- [lqr](lqr.md) -- Linear-quadratic regulator (infinite, finite, time-varying, integral action)
- [dare](dare.md) -- Discrete algebraic Riccati equation solver (complex Schur method)
- [place](place.md) -- Pole placement via Ackermann's formula

## When to Use

Pick **PID** when you have a single-loop SISO or MIMO regulation problem and want
to compose exactly the features you need (anti-windup, derivative filtering,
feed-forward) without paying for what you don't.

Pick **LQR** when you have a state-space model and want optimal full-state feedback
that minimises a quadratic cost. LQR calls DARE internally.

Pick **DARE** directly when you need the solution to the discrete algebraic Riccati
equation outside of LQR (e.g., for terminal cost computation in MPC).

Pick **place** when you need direct pole assignment rather than cost-based tuning.

## Theory

- [PID Theory](../reference/pid-theory.md) -- PID mathematical background
- [Optimal Control Theory](../reference/optimal-control-theory.md) -- LQR, DARE, Riccati equations
