# Control

Classical and optimal control types for feedback and feedforward control design.
PID handles the majority of single-loop regulation problems, while LQR provides
optimal state-feedback for multi-variable systems. DARE and place are the
underlying design tools that LQR builds on.

## Types

- [PID](pid/README.md)<br/> Policy-based PID controller with compile-time feature composition
- [lqr](lqr.md)<br/> Linear-quadratic regulator (infinite, finite, time-varying, integral action)
- [dare](dare.md)<br/> Discrete algebraic Riccati equation solver (complex Schur method)
- [place](place.md)<br/> Pole placement via Ackermann's formula
- [mrac](mrac.md)<br/> Model reference adaptive controller with dead-zone, sigma-modification, and e-modification robustification
- [l1](l1.md)<br/> L1 adaptive controller with state predictor, projection-based adaptation, and low-pass filtered control

## When to Use

Pick **PID** when you have a single-loop SISO or MIMO regulation problem and want
to compose exactly the features you need (anti-windup, derivative filtering,
feed-forward) without paying for what you don't.

Pick **LQR** when you have a state-space model and want optimal full-state feedback
that minimises a quadratic cost. LQR calls DARE internally.

Pick **DARE** directly when you need the solution to the discrete algebraic Riccati
equation outside of LQR (e.g., for terminal cost computation in MPC).

Pick **place** when you need direct pole assignment rather than cost-based tuning.

Pick **MRAC** when you have an unknown or uncertain plant and want to track a
reference model adaptively. MRAC adapts online without requiring a plant model,
supporting both SISO and MIMO configurations with compile-time robustification
policy selection.

Pick **L1** when you need guaranteed bounded transient performance independent of the
adaptation gain. L1 decouples adaptation speed from robustness through its low-pass
filter &mdash; suitable when fast adaptation and predictable transient behaviour are both
required.

## Theory

- [PID Theory](../../background/pid.md)<br/> PID mathematical background
- [Adaptive Control Theory](../../background/adaptive-control.md)<br/> MRAC and L1 mathematical background

