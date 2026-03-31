# ctrlpp Documentation

Guides and API reference for the ctrlpp C++23 control library.

## Getting Started

- [Getting Started](getting-started.md)<br/> Install ctrlpp and run your first PID controller

## Guides

### Introduction

- [Your First PID](guides/intro/your-first-pid.md)<br/> A minimal PID controller from scratch
- [Your First Estimator](guides/intro/your-first-estimator.md)<br/> Add an observer to your control loop
- [Your First MPC](guides/intro/your-first-mpc.md)<br/> Model predictive control in under 30 lines

### Detailed tutorials

- [PID Composition](guides/pid/composition.md)<br/> Policy-based composition for anti-windup, filtering, and more
- [Cascade Control](guides/pid/cascade.md)<br/> Inner/outer loop cascade PID
- [Solver Injection](guides/mpc/solver-injection.md)<br/> Concept-based QP/NLP solver backends
- [Observer-Controller Patterns](guides/estimation/observer-controller.md)<br/> Composing observers with controllers
- [Sysid Workflow](guides/sysid/workflow.md)<br/> Identify, model, and control
- [Composition Patterns](guides/patterns/composition.md)<br/> Cross-cutting design patterns
- [Point-to-Point Motion](guides/trajectory/point-to-point.md)<br/> Choosing the right trajectory profile
- [Multi-Waypoint Paths](guides/trajectory/multi-waypoint.md)<br/> Splines and B-splines for waypoint sequences
- [Real-Time Replanning](guides/trajectory/real-time.md)<br/> Online trajectory planners
- [Multi-Axis Coordination](guides/trajectory/multi-axis.md)<br/> Synchronizing multiple axes

## API Reference

### Control

- [PID](api/control/pid/README.md)<br/> Policy-based PID controller
- [lqr](api/control/lqr.md)<br/> Linear-quadratic regulator
- [dare](api/control/dare.md)<br/> Discrete algebraic Riccati equation solver
- [place](api/control/place.md)<br/> Pole placement

### Estimation

- [kalman](api/estimation/kalman.md)<br/> Linear Kalman filter
- [luenberger](api/estimation/luenberger.md)<br/> Luenberger observer
- [ekf](api/estimation/ekf.md)<br/> Extended Kalman filter
- [ukf](api/estimation/ukf.md)<br/> Unscented Kalman filter
- [particle_filter](api/estimation/particle-filter.md)<br/> Bootstrap SIR particle filter
- [mekf](api/estimation/mekf.md)<br/> Multiplicative extended Kalman filter (SO(3))
- [manifold_ukf](api/estimation/manifold-ukf.md)<br/> Manifold unscented Kalman filter (SO(3))
- [complementary_filter](api/estimation/complementary-filter.md)<br/> Mahony complementary filter
- [observer_policy](api/estimation/observer-policy.md)<br/> Observer concept interface

### MPC and MHE

- [mpc](api/mpc/mpc.md)<br/> Linear model predictive control (OSQP)
- [nmpc](api/mpc/nmpc.md)<br/> Nonlinear model predictive control (NLopt)
- [mhe](api/mpc/mhe.md)<br/> Linear moving horizon estimation (OSQP)
- [nmhe](api/mpc/nmhe.md)<br/> Nonlinear moving horizon estimation (NLopt)
- [osqp_solver](api/mpc/osqp-solver.md)<br/> OSQP QP solver backend
- [nlopt_solver](api/mpc/nlopt-solver.md)<br/> NLopt NLP solver backend

### Signal Processing

- [biquad](api/dsp/biquad.md)<br/> IIR second-order section filter
- [fir](api/dsp/fir.md)<br/>Finite impulse response filter
- [discrete_filter](api/dsp/discrete-filter.md)<br/>Discrete filter concept

### System Identification

- [rls](api/sysid/rls.md)<br/>Recursive least squares
- [batch_arx](api/sysid/batch-arx.md)<br/>Batch ARX identification (QR)
- [recursive_arx](api/sysid/recursive-arx.md)<br/>Recursive ARX identification
- [n4sid](api/sysid/n4sid.md)<br/>Subspace identification (BDCSVD)
- [fit_metrics](api/sysid/fit-metrics.md)<br/>Goodness-of-fit metrics (NRMSE, VAF)
- [sysid_result](api/sysid/sysid-result.md)<br/>Identification result container

### Lie Groups

- [so3](api/lie/so3.md)<br/>SO(3) quaternion rotation utilities

### Model Utilities

- [state_space](api/model/state-space.md)<br/>Linear state-space model
- [transfer_function](api/model/transfer-function.md)<br/>Transfer function representation
- [discretise](api/model/discretise.md)<br/>Continuous-to-discrete conversion
- [conversion](api/model/conversion.md)<br/>Transfer function / state-space conversion
- [analysis](api/model/analysis.md)<br/>Stability and controllability analysis
- [propagate](api/model/propagate.md)<br/>State propagation utilities
- [dynamics_model](api/model/dynamics-model.md)<br/>Dynamics model concept
- [measurement_model](api/model/measurement-model.md)<br/>Measurement model concept
- [differentiable_dynamics](api/model/differentiable-dynamics.md)<br/>Differentiable dynamics concept
- [differentiable_measurement](api/model/differentiable-measurement.md)<br/>Differentiable measurement concept
- [constraint_model](api/model/constraint-model.md)<br/>Constraint model concept

### Trajectory

- [Trajectory API Overview](api/trajectory/README.md)<br/>Complete trajectory module index
- [trajectory](api/trajectory/trajectory.md)<br/>Trajectory convenience header and concept
- [trajectory_segment](api/trajectory/trajectory-segment.md)<br/>Trajectory segment concept definition
- [trajectory_types](api/trajectory/trajectory-types.md)<br/>trajectory_point, Vector type aliases
- [path_segment](api/trajectory/path-segment.md)<br/>Path segment concept
- [piecewise_path](api/trajectory/piecewise-path.md)<br/>Piecewise path composition
- [piecewise_trajectory](api/trajectory/piecewise-trajectory.md)<br/>Piecewise trajectory composition
- [time_scaling](api/trajectory/time-scaling.md)<br/>Time scaling functions
- [trapezoidal_trajectory](api/trajectory/trapezoidal-trajectory.md)<br/>Trapezoidal velocity profile
- [double_s_trajectory](api/trajectory/double-s-trajectory.md)<br/>Double-S (jerk-limited) velocity profile
- [modified_trap_trajectory](api/trajectory/modified-trap-trajectory.md)<br/>Modified trapezoidal velocity profile
- [modified_sin_trajectory](api/trajectory/modified-sin-trajectory.md)<br/>Modified sinusoidal velocity profile
- [cubic_path](api/trajectory/cubic-path.md)<br/>Cubic polynomial path
- [cubic_trajectory](api/trajectory/cubic-trajectory.md)<br/>Cubic polynomial trajectory
- [quintic_path](api/trajectory/quintic-path.md)<br/>Quintic polynomial path
- [quintic_trajectory](api/trajectory/quintic-trajectory.md)<br/>Quintic polynomial trajectory
- [septic_path](api/trajectory/septic-path.md)<br/>Septic polynomial path
- [septic_trajectory](api/trajectory/septic-trajectory.md)<br/>Septic polynomial trajectory
- [harmonic_path](api/trajectory/harmonic-path.md)<br/>Harmonic path segment
- [cycloidal_path](api/trajectory/cycloidal-path.md)<br/>Cycloidal path segment
- [cubic_spline](api/trajectory/cubic-spline.md)<br/>Cubic spline interpolation (natural, clamped, periodic)
- [smoothing_spline](api/trajectory/smoothing-spline.md)<br/>Smoothing spline approximation
- [bspline_trajectory](api/trajectory/bspline-trajectory.md)<br/>B-spline trajectory with compile-time degree
- [online_planner_2nd](api/trajectory/online-planner-2nd.md)<br/>2nd-order online trajectory planner
- [online_planner_3rd](api/trajectory/online-planner-3rd.md)<br/>3rd-order online trajectory planner
- [synchronize](api/trajectory/synchronize.md)<br/>Multi-axis trajectory synchronization

## Validation

- [Validation Status](validation.md)<br/>Testing levels and cross-validation against GNU Octave

## Benchmarks

- [Benchmarks](benchmarks.md)<br/>Internal performance and comparison benchmarks

## Background Theory

Standalone theory and mathematical background for the algorithms in ctrlpp.

- [PID Theory](background/pid.md)<br/>Parallel form, derivative filter, anti-windup clamping and back-calculation
- [Kalman Theory](background/kalman.md)<br/>Linear Kalman filter predict/update equations, optimality, innovation
- [EKF Theory](background/ekf.md)<br/>Extended Kalman filter: Jacobian linearization, prediction and update equations, numerical stability
- [UKF Theory](background/ukf.md)<br/>Unscented Kalman filter: sigma-point generation, weights, unscented transform
- [Particle Filter Theory](background/particle-filter.md)<br/>Importance sampling, weight update, systematic resampling, ESS
- [Attitude Estimation Theory](background/attitude-estimation.md)<br/>Quaternion kinematics, MEKF error-state, manifold UKF
- [MPC Theory](background/mpc.md)<br/>QP/NLP optimization formulation, terminal cost and constraints, stability
- [MHE Theory](background/mhe.md)<br/>Moving horizon cost function, arrival cost approximation, duality with MPC
- [System Identification Theory](background/sysid.md)<br/>ARX regression model, RLS recursive update, N4SID Hankel matrix
- [DSP Theory](background/dsp.md)<br/>Biquad transfer function, bilinear transform, FIR convolution, cascading
- [Trajectory Generation Theory](background/trajectory-generation.md)<br/>Polynomial trajectories, velocity profiles, splines, online planners, synchronization
- [LQR Theory](background/lqr.md)<br/>Linear quadratic regulator: DARE, cost function, optimal gain
- [Observers Theory](background/observers.md)<br/>State observer design: Luenberger, pole placement, duality
- [NMPC Theory](background/nmpc.md)<br/>Nonlinear MPC: NLP formulation, direct multiple shooting, terminal constraints
- [NMHE Theory](background/nmhe.md)<br/>Nonlinear MHE: arrival cost, nonlinear optimization, EKF approximation
- [SO(3) Theory](background/so3.md)<br/>Special orthogonal group: rotation representations, quaternion algebra, exponential map
