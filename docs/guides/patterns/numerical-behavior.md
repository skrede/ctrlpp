# Numerical Behavior

This page documents the library's contract around numerical edge cases so users know what to expect and how to build robust systems on top of it.

## NaN propagation

ctrlpp does **not** silently clamp, saturate, or replace NaN/Inf values that
arise from user-supplied inputs. If you feed a controller or estimator inputs
that produce intermediate overflow (for example, a PID with `Kp = 1e200` and
`error = 1e200`), the output will be NaN or Inf.

This is intentional. In a real-time control system, a NaN in the actuator
command is a loud signal that something upstream is broken. Silently clamping
the output to a "safe" value would mask the root cause and make debugging
harder.

**Your responsibility:** validate inputs at system boundaries (sensor readings,
setpoints, configuration parameters) before passing them to ctrlpp. The library
trusts that inputs are mathematically reasonable.

## What the library does guard against

ctrlpp hardens its internals against degenerate-but-valid inputs that would
otherwise produce silent corruption through intermediate overflow:

- **Spline near-zero spans:** B-spline and cubic spline evaluation uses
  relative epsilon thresholds (scaled by knot magnitude) rather than absolute
  machine epsilon, preventing overflow when adjacent knots are very close but
  not identical.

- **Thomas solver pivot collapse:** The tridiagonal solver detects near-zero
  pivots during forward elimination and returns zero rather than producing
  Inf from division.

- **DARE symplectic overflow:** The discrete algebraic Riccati equation solver
  checks that the symplectic matrix is finite before Schur decomposition, and
  that the extracted solution P is finite before eigenvalue validation. Returns
  `std::nullopt` for degenerate inputs.

- **L1 DC gain inversion:** The L1 adaptive controller validates that the
  predictor model's DC gain is invertible before computing the feedforward gain
  `K_r`. Falls back to identity if the inversion is ill-conditioned.

- **Particle filter weight degeneracy:** When all particles have negligible
  likelihood (complete weight collapse), the log-weight normalizer resets to
  uniform weights rather than producing NaN from `-inf - (-inf)`.

- **UKF covariance repair:** The unscented Kalman filter enforces positive
  semi-definiteness of the state covariance after the update step by clamping
  negative eigenvalues to zero. This prevents the covariance from going
  indefinite due to the negative Merwe sigma point weights.

- **N4SID degenerate data:** The subspace identification algorithm checks
  finiteness of extracted system matrices and zeros out non-finite results
  rather than propagating NaN through downstream computations.

- **RLS overflow guard:** The recursive least squares estimator skips the
  update step when the quadratic form `phi^T * P * phi` overflows, preventing
  NaN from `0 * Inf` in the gain computation.

## Summary

| Layer | Behavior |
|---|---|
| Inputs (your code) | Validate at system boundaries |
| Constructor/config | Library rejects or falls back on degenerate configs |
| Algorithm internals | Library uses numerically stable formulations |
| Outputs | NaN/Inf propagates faithfully, never silently clamped |
