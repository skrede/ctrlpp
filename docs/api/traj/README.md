# Trajectory API Reference

API reference for trajectory generation primitives. All headers are under `ctrlpp/traj/`.

## Foundation

Core types, concepts, and composition primitives.

| Page | Header | Description |
|------|--------|-------------|
| [trajectory](trajectory.md) | `trajectory.h` | Adapter lifting normalized paths into trajectory segments |
| [trajectory-segment](trajectory-segment.md) | `trajectory_segment.h` | Concept constraining trajectory segment types |
| [trajectory-types](trajectory-types.md) | `trajectory_types.h` | Core output types (`trajectory_point`, `path_point`) |
| [path-segment](path-segment.md) | `path_segment.h` | Concept constraining normalized path segment types |
| [piecewise-path](piecewise-path.md) | `piecewise_path.h` | Variadic heterogeneous path composition |
| [piecewise-trajectory](piecewise-trajectory.md) | `piecewise_trajectory.h` | Variadic heterogeneous trajectory composition |

## Elementary -- Polynomial

Normalized motion laws based on polynomial basis functions.

| Page | Header | Description |
|------|--------|-------------|
| [cubic-path](cubic-path.md) | `cubic_path.h` | Cubic polynomial path (C1) |
| [cubic-trajectory](cubic-trajectory.md) | `cubic_trajectory.h` | Cubic polynomial with velocity BCs |
| [quintic-path](quintic-path.md) | `quintic_path.h` | Quintic polynomial path (C2) |
| [quintic-trajectory](quintic-trajectory.md) | `quintic_trajectory.h` | Quintic polynomial with velocity+acceleration BCs |
| [septic-path](septic-path.md) | `septic_path.h` | Septic polynomial path (C3) |
| [septic-trajectory](septic-trajectory.md) | `septic_trajectory.h` | Septic polynomial with velocity+acceleration+jerk BCs |

## Elementary -- Trigonometric

Normalized motion laws based on trigonometric basis functions.

| Page | Header | Description |
|------|--------|-------------|
| [harmonic-path](harmonic-path.md) | `harmonic_path.h` | Harmonic (sinusoidal velocity) motion law |
| [cycloidal-path](cycloidal-path.md) | `cycloidal_path.h` | Cycloidal motion law (zero endpoint acceleration) |

## Composite

Multi-phase velocity profiles solving for timing from kinematic limits.

| Page | Header | Description |
|------|--------|-------------|
| [trapezoidal-trajectory](trapezoidal-trajectory.md) | `trapezoidal_trajectory.h` | Trapezoidal velocity (LSPB) with degenerate handling |
| [double-s-trajectory](double-s-trajectory.md) | `double_s_trajectory.h` | Double-S (7-segment) jerk-limited profile |
| [modified-trap-trajectory](modified-trap-trajectory.md) | `modified_trap_trajectory.h` | Modified trapezoidal with cycloidal acceleration |
| [modified-sin-trajectory](modified-sin-trajectory.md) | `modified_sin_trajectory.h` | Modified sinusoidal with harmonic+cycloidal blend |

## Operations

Utilities for computing timing and synchronization.

| Page | Header | Description |
|------|--------|-------------|
| [time-scaling](time-scaling.md) | `time_scaling.h` | Kinematic time scaling from path derivatives and limits |

## See Also

- [Trajectory Generation Theory](../../reference/trajectory-generation.md)
