# Trajectory API Reference

Trajectory generation primitives for point-to-point motion, multi-point interpolation, and multi-axis coordination. All trajectory types are domain-agnostic (scalar or vector-valued, no spatial awareness).

## Infrastructure

Core types, concepts, and building blocks used by all trajectory primitives.

- [trajectory](trajectory.md)<br/> Trajectory convenience header and concept
- [trajectory-segment](trajectory-segment.md)<br/> `trajectory_segment` concept definition
- [trajectory-types](trajectory-types.md)<br/> `trajectory_point`, `Vector` type aliases
- [path-segment](path-segment.md)<br/> `path_segment` concept for geometric paths
- [piecewise-path](piecewise-path.md)<br/> Piecewise path composition from path segments
- [piecewise-trajectory](piecewise-trajectory.md)<br/> Piecewise trajectory composition from trajectory segments
- [time-scaling](time-scaling.md)<br/> Time scaling functions for path-to-trajectory conversion

## Elementary Profiles

Single-segment point-to-point motion primitives with analytical velocity profiles.

- [trapezoidal-trajectory](trapezoidal-trajectory.md)<br/> Trapezoidal (bang-coast-bang) velocity profile
- [double-s-trajectory](double-s-trajectory.md)<br/> Double-S (jerk-limited) velocity profile
- [modified-trap-trajectory](modified-trap-trajectory.md)<br/> Modified trapezoidal velocity profile
- [modified-sin-trajectory](modified-sin-trajectory.md)<br/> Modified sinusoidal velocity profile

## Polynomial Paths and Trajectories

Polynomial motion primitives defined by boundary conditions on position and derivatives.

- [cubic-path](cubic-path.md)<br/> Cubic polynomial path (position + velocity BCs)
- [cubic-trajectory](cubic-trajectory.md)<br/>Cubic polynomial trajectory (time-parametrized)
- [quintic-path](quintic-path.md)<br/>Quintic polynomial path (position + velocity + acceleration BCs)
- [quintic-trajectory](quintic-trajectory.md)<br/>Quintic polynomial trajectory
- [septic-path](septic-path.md)<br/>Septic polynomial path (up to jerk BCs)
- [septic-trajectory](septic-trajectory.md)<br/>Septic polynomial trajectory
- [harmonic-path](harmonic-path.md)<br/>Harmonic (sinusoidal) path segment
- [cycloidal-path](cycloidal-path.md)<br/>Cycloidal path segment

## Multipoint Trajectories

Spline and B-spline interpolation through multiple waypoints.

- [cubic-spline](cubic-spline.md)<br/>Cubic spline interpolation (natural, clamped, periodic BCs)
- [smoothing-spline](smoothing-spline.md)<br/>Smoothing spline approximation with mu tradeoff
- [bspline-trajectory](bspline-trajectory.md)<br/>B-spline trajectory with compile-time degree

## Online Planners

Real-time trajectory filters for dynamic target tracking in control loops.

- [online-planner-2nd](online-planner-2nd.md)<br/>2nd-order planner (velocity + acceleration limits)
- [online-planner-3rd](online-planner-3rd.md)<br/>3rd-order planner (velocity + acceleration + jerk limits)
- `online_planner_diagnostics`<br/>Shared disposition report for both planners: which profile the last `update` built, against the one it was commanded. Documented on the [3rd-order page](online-planner-3rd.md#substitution-reporting) and the [2nd-order page](online-planner-2nd.md#substitution-reporting); header `ctrlpp/trajectory/online_planner_diagnostics.h`

## Operations

Multi-axis coordination and trajectory manipulation.

- [synchronize](synchronize.md)<br/>Multi-axis synchronization: every axis is checked at the slowest axis's duration before any axis is retimed

## See Also

- [Trajectory Generation Theory](../../background/trajectory-generation.md)<br/>Mathematical background for all trajectory algorithms
