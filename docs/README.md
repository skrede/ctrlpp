# ctrlpp Documentation

Guides and API reference for the ctrlpp C++20 control library.

## Getting Started

- [Getting Started](getting-started.md) -- Install ctrlpp and run your first PID controller

## Guides

### Introduction

- [Your First PID](guides/intro/your-first-pid.md) -- A minimal PID controller from scratch
- [Your First Estimator](guides/intro/your-first-estimator.md) -- Add an observer to your control loop
- [Your First MPC](guides/intro/your-first-mpc.md) -- Model predictive control in under 30 lines

### Deep Dives

- [PID Composition](guides/pid/composition.md) -- Policy-based composition for anti-windup, filtering, and more
- [Cascade Control](guides/pid/cascade.md) -- Inner/outer loop cascade PID
- [Solver Injection](guides/mpc/solver-injection.md) -- Concept-based QP/NLP solver backends
- [Observer-Controller Patterns](guides/estimation/observer-controller.md) -- Composing observers with controllers
- [Sysid Workflow](guides/sysid/workflow.md) -- Identify, model, and control
- [Composition Patterns](guides/patterns/composition.md) -- Cross-cutting design patterns

## API Reference

### Control

- [PID](control/pid/README.md) -- Policy-based PID controller
- [lqr](control/lqr.md) -- Linear-quadratic regulator
- [dare](control/dare.md) -- Discrete algebraic Riccati equation solver
- [place](control/place.md) -- Pole placement

### Estimation

- [kalman](estimation/kalman.md) -- Linear Kalman filter
- [luenberger](estimation/luenberger.md) -- Luenberger observer
- [ekf](estimation/ekf.md) -- Extended Kalman filter
- [ukf](estimation/ukf.md) -- Unscented Kalman filter
- [particle_filter](estimation/particle-filter.md) -- Bootstrap SIR particle filter
- [mekf](estimation/mekf.md) -- Multiplicative extended Kalman filter (SO(3))
- [manifold_ukf](estimation/manifold-ukf.md) -- Manifold unscented Kalman filter (SO(3))
- [complementary_filter](estimation/complementary-filter.md) -- Mahony complementary filter
- [observer_policy](estimation/observer-policy.md) -- Observer concept interface

### MPC and MHE

- [mpc](mpc/mpc.md) -- Linear model predictive control (OSQP)
- [nmpc](mpc/nmpc.md) -- Nonlinear model predictive control (NLopt)
- [mhe](mpc/mhe.md) -- Linear moving horizon estimation (OSQP)
- [nmhe](mpc/nmhe.md) -- Nonlinear moving horizon estimation (NLopt)
- [osqp_solver](mpc/osqp-solver.md) -- OSQP QP solver backend
- [nlopt_solver](mpc/nlopt-solver.md) -- NLopt NLP solver backend

### Signal Processing

- [biquad](dsp/biquad.md) -- IIR second-order section filter
- [fir](dsp/fir.md) -- Finite impulse response filter
- [discrete_filter](dsp/discrete-filter.md) -- Cascaded biquad filter

### System Identification

- [rls](sysid/rls.md) -- Recursive least squares
- [batch_arx](sysid/batch-arx.md) -- Batch ARX identification (QR)
- [recursive_arx](sysid/recursive-arx.md) -- Recursive ARX identification
- [n4sid](sysid/n4sid.md) -- Subspace identification (BDCSVD)
- [fit_metrics](sysid/fit-metrics.md) -- Goodness-of-fit metrics (NRMSE, VAF)
- [sysid_result](sysid/sysid-result.md) -- Identification result container

### Lie Groups

- [so3](lie/so3.md) -- SO(3) quaternion rotation utilities

### Model Utilities

- [state_space](model/state-space.md) -- Linear state-space model
- [transfer_function](model/transfer-function.md) -- Transfer function representation
- [discretise](model/discretise.md) -- Continuous-to-discrete conversion
- [conversion](model/conversion.md) -- Transfer function / state-space conversion
- [analysis](model/analysis.md) -- Stability and controllability analysis
- [propagate](model/propagate.md) -- State propagation utilities
- [dynamics_model](model/dynamics-model.md) -- Dynamics model concept
- [measurement_model](model/measurement-model.md) -- Measurement model concept
- [differentiable_dynamics](model/differentiable-dynamics.md) -- Differentiable dynamics concept
- [differentiable_measurement](model/differentiable-measurement.md) -- Differentiable measurement concept
- [constraint_model](model/constraint-model.md) -- Constraint model concept

## Reference

- [Reference](reference/README.md) -- Theory and mathematical background for all algorithms
