# Trajectory API Reference

Trajectory generation primitives for point-to-point motion, multi-point interpolation, and multi-axis coordination. All trajectory types are domain-agnostic (scalar or vector-valued, no spatial awareness).

## Infrastructure

Core types, concepts, and building blocks used by all trajectory primitives.

- [trajectory](trajectory.md) -- Trajectory convenience header and concept
- [trajectory-segment](trajectory-segment.md) -- `trajectory_segment` concept definition
- [trajectory-types](trajectory-types.md) -- `trajectory_point`, `Vector` type aliases
- [path-segment](path-segment.md) -- `path_segment` concept for geometric paths
- [piecewise-path](piecewise-path.md) -- Piecewise path composition from path segments
- [piecewise-trajectory](piecewise-trajectory.md) -- Piecewise trajectory composition from trajectory segments
- [time-scaling](time-scaling.md) -- Time scaling functions for path-to-trajectory conversion

## Elementary Profiles

Single-segment point-to-point motion primitives with analytical velocity profiles.

- [trapezoidal-trajectory](trapezoidal-trajectory.md) -- Trapezoidal (bang-coast-bang) velocity profile
- [double-s-trajectory](double-s-trajectory.md) -- Double-S (jerk-limited) velocity profile
- [modified-trap-trajectory](modified-trap-trajectory.md) -- Modified trapezoidal velocity profile
- [modified-sin-trajectory](modified-sin-trajectory.md) -- Modified sinusoidal velocity profile

## Polynomial Paths and Trajectories

Polynomial motion primitives defined by boundary conditions on position and derivatives.

- [cubic-path](cubic-path.md) -- Cubic polynomial path (position + velocity BCs)
- [cubic-trajectory](cubic-trajectory.md) -- Cubic polynomial trajectory (time-parameterised)
- [quintic-path](quintic-path.md) -- Quintic polynomial path (position + velocity + acceleration BCs)
- [quintic-trajectory](quintic-trajectory.md) -- Quintic polynomial trajectory
- [septic-path](septic-path.md) -- Septic polynomial path (up to jerk BCs)
- [septic-trajectory](septic-trajectory.md) -- Septic polynomial trajectory
- [harmonic-path](harmonic-path.md) -- Harmonic (sinusoidal) path segment
- [cycloidal-path](cycloidal-path.md) -- Cycloidal path segment

## Multipoint Trajectories

Spline and B-spline interpolation through multiple waypoints.

- [cubic-spline](cubic-spline.md) -- Cubic spline interpolation (natural, clamped, periodic BCs)
- [smoothing-spline](smoothing-spline.md) -- Smoothing spline approximation with mu tradeoff
- [bspline-trajectory](bspline-trajectory.md) -- B-spline trajectory with compile-time degree

## Online Planners

Real-time trajectory filters for dynamic target tracking in control loops.

- [online-planner-2nd](online-planner-2nd.md) -- 2nd-order planner (velocity + acceleration limits)
- [online-planner-3rd](online-planner-3rd.md) -- 3rd-order planner (velocity + acceleration + jerk limits)

## Operations

Multi-axis coordination and trajectory manipulation.

- [synchronize](synchronize.md) -- Multi-axis synchronization (rescale to slowest axis)

## See Also

- [Trajectory Generation Theory](../../background/trajectory-generation.md) -- Mathematical background for all trajectory algorithms
