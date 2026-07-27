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
  that the extracted solution P is finite before eigenvalue validation. It
  returns `ctrlpp::expected<dare_result<Scalar, NX>, dare_error>`, so a
  degenerate input comes back as a named rejection rather than a bare empty
  result: `dare_error::non_finite_input` when A, B, Q, R or the assembled
  symplectic Z contains NaN or Inf.

- **L1 DC gain inversion:** The L1 adaptive controller validates that the
  predictor model's DC gain is invertible before computing the feedforward gain
  `K_r`. Construction goes through the fallible factory
  `l1_controller::create`, which returns
  `ctrlpp::expected<l1_controller, l1_error>` and reports
  `l1_error::singular_predictor` when (I - A_m) is singular,
  `l1_error::singular_dc_gain` when the DC gain (I - A_m)^{-1} B is singular or
  non-finite, and `l1_error::non_finite_gain` when `K_r` is non-finite. It does
  not throw.

- **Particle filter weight degeneracy:** When all particles have negligible
  likelihood (complete weight collapse), the log-weight normalizer resets to
  uniform weights rather than producing NaN from `-inf - (-inf)`.

- **UKF covariance update:** The unscented Kalman filter uses the algebraically
  complete minimum mean-square-error reduction `P = P - K*S*K^T`, where
  `S = Pzz + R` and `K = Pxz*S^{-1}`. Because `K*S*K^T = K*Pxz^T`, this term is
  exactly the uncertainty the measurement removes and no extra `K*R*K^T` term is
  added. The sigma points are built from a permutation-correct covariance square
  root, so the reduction stays symmetric positive semidefinite without requiring
  an eigendecomposition.

- **MOESP degenerate data:** The subspace identification algorithm checks
  finiteness of extracted system matrices. When data is rank-deficient,
  returns a result with `condition_number = infinity` and default-initialized
  system matrices so the caller can detect and handle the failure.

- **RLS overflow guard:** The recursive least squares estimator skips the
  update step when the quadratic form `phi^T * P * phi` overflows, preventing
  NaN from `0 * Inf` in the gain computation. The `update()` method returns
  `false` when this occurs so the caller can detect skipped updates.

## Summary

| Layer | Behavior |
|---|---|
| Inputs (your code) | Validate at system boundaries |
| Constructor/config | Library rejects or falls back on degenerate configs |
| Algorithm internals | Library uses numerically stable formulations |
| Outputs | NaN/Inf propagates faithfully, never silently clamped |

The outputs row is being narrowed for the per-step surfaces currently under
conversion: where a surface gains a typed failure return, a degenerate step is
reported through that return instead of being left to propagate as NaN. See
[error-reporting.md](error-reporting.md) for the two reporting channels and for
which surfaces have been converted so far.
