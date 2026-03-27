# ctrlpp Documentation

Guides and API reference for the ctrlpp C++23 control library.

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

- [PID](api/control/pid/README.md) -- Policy-based PID controller
- [lqr](api/control/lqr.md) -- Linear-quadratic regulator
- [dare](api/control/dare.md) -- Discrete algebraic Riccati equation solver
- [place](api/control/place.md) -- Pole placement

### Estimation

- [kalman](api/estimation/kalman.md) -- Linear Kalman filter
- [luenberger](api/estimation/luenberger.md) -- Luenberger observer
- [ekf](api/estimation/ekf.md) -- Extended Kalman filter
- [ukf](api/estimation/ukf.md) -- Unscented Kalman filter
- [particle_filter](api/estimation/particle-filter.md) -- Bootstrap SIR particle filter
- [mekf](api/estimation/mekf.md) -- Multiplicative extended Kalman filter (SO(3))
- [manifold_ukf](api/estimation/manifold-ukf.md) -- Manifold unscented Kalman filter (SO(3))
- [complementary_filter](api/estimation/complementary-filter.md) -- Mahony complementary filter
- [observer_policy](api/estimation/observer-policy.md) -- Observer concept interface

### MPC and MHE

- [mpc](api/mpc/mpc.md) -- Linear model predictive control (OSQP)
- [nmpc](api/mpc/nmpc.md) -- Nonlinear model predictive control (NLopt)
- [mhe](api/mpc/mhe.md) -- Linear moving horizon estimation (OSQP)
- [nmhe](api/mpc/nmhe.md) -- Nonlinear moving horizon estimation (NLopt)
- [osqp_solver](api/mpc/osqp-solver.md) -- OSQP QP solver backend
- [nlopt_solver](api/mpc/nlopt-solver.md) -- NLopt NLP solver backend

### Signal Processing

- [biquad](api/dsp/biquad.md) -- IIR second-order section filter
- [fir](api/dsp/fir.md) -- Finite impulse response filter
- [discrete_filter](api/dsp/discrete-filter.md) -- Discrete filter concept

### System Identification

- [rls](api/sysid/rls.md) -- Recursive least squares
- [batch_arx](api/sysid/batch-arx.md) -- Batch ARX identification (QR)
- [recursive_arx](api/sysid/recursive-arx.md) -- Recursive ARX identification
- [n4sid](api/sysid/n4sid.md) -- Subspace identification (BDCSVD)
- [fit_metrics](api/sysid/fit-metrics.md) -- Goodness-of-fit metrics (NRMSE, VAF)
- [sysid_result](api/sysid/sysid-result.md) -- Identification result container

### Lie Groups

- [so3](api/lie/so3.md) -- SO(3) quaternion rotation utilities

### Model Utilities

- [state_space](api/model/state-space.md) -- Linear state-space model
- [transfer_function](api/model/transfer-function.md) -- Transfer function representation
- [discretise](api/model/discretise.md) -- Continuous-to-discrete conversion
- [conversion](api/model/conversion.md) -- Transfer function / state-space conversion
- [analysis](api/model/analysis.md) -- Stability and controllability analysis
- [propagate](api/model/propagate.md) -- State propagation utilities
- [dynamics_model](api/model/dynamics-model.md) -- Dynamics model concept
- [measurement_model](api/model/measurement-model.md) -- Measurement model concept
- [differentiable_dynamics](api/model/differentiable-dynamics.md) -- Differentiable dynamics concept
- [differentiable_measurement](api/model/differentiable-measurement.md) -- Differentiable measurement concept
- [constraint_model](api/model/constraint-model.md) -- Constraint model concept

### Trajectory (coming soon)

Trajectory generation primitives will be documented in [api/traj/](api/traj/) once implemented.

## Background Theory

Standalone theory and mathematical background for the algorithms in ctrlpp.

- [PID Theory](background/pid.md) -- Parallel form, derivative filter, anti-windup clamping and back-calculation
- [Kalman Theory](background/kalman.md) -- Linear Kalman filter predict/update equations, optimality, innovation
- [EKF/UKF Theory](background/ekf-ukf.md) -- EKF Jacobian linearisation, UKF sigma-point generation and weights
- [Particle Filter Theory](background/particle-filter.md) -- Importance sampling, weight update, systematic resampling, ESS
- [Attitude Estimation Theory](background/attitude-estimation.md) -- Quaternion kinematics, MEKF error-state, manifold UKF
- [MPC Theory](background/mpc.md) -- QP/NLP optimisation formulation, terminal cost and constraints, stability
- [MHE Theory](background/mhe.md) -- Moving horizon cost function, arrival cost approximation, duality with MPC
- [System Identification Theory](background/sysid.md) -- ARX regression model, RLS recursive update, N4SID Hankel matrix
- [DSP Theory](background/dsp.md) -- Biquad transfer function, bilinear transform, FIR convolution, cascading
