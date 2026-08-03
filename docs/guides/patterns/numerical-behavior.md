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

### Where that contract has been narrowed

The paragraph above no longer describes the whole library, and the exceptions
are named individually rather than as a category. A surface that carries a
carried state, rather than merely producing an output, cannot honor
faithful propagation and remain usable: one poisoned sample destroys the memory
permanently, and every later step then produces NaN from valid data. Those
surfaces **reject** instead.

**Reject a non-finite operand before mutating anything** (a rejected step leaves
the carried estimate bitwise unchanged and reports the specific cause on the
failure channel):

| Surface | Return |
|---|---|
| `kalman_filter::update` | `expected<void, kalman_update_error>` |
| `ekf::update` | `expected<void, ekf_update_error>` |
| `ukf::update` | `expected<void, ukf_update_error>` |
| `mekf::update` | `expected<void, mekf_update_error>` |
| `manifold_ukf::update` | `expected<void, manifold_ukf_update_error>` |
| `luenberger_observer::update` | `expected<void, luenberger_update_error>` |
| `mhe::update` | `expected<void, ekf_update_error>` (forwarded from the embedded filter) |
| `nmhe::update` | `expected<void, ekf_update_error>` (forwarded from the embedded filter) |
| `complementary_filter::update`, all three overloads | `expected<void, cf_update_error>` |
| `pid::compute`, both overloads | `expected<vector_t, pid_step_error>` |
| `mrac_controller::evaluate` | `expected<input_type, mrac_step_error>` |
| `l1_controller::evaluate` | `expected<input_type, l1_step_error>` |

The last three carry a payload where the estimator updates carry none, because a
controller step produces the command it was asked for. A refusal there means the
caller has **no control output for this cycle**, which is a materially different
situation from an estimator declining a measurement and leaving its estimate
standing: an actuator will be driven by something regardless. The caller decides
what -- hold the last successful command, drive a configured safe value, or fail
over -- and the controller does not choose, because the right answer is a
property of the plant.

Two of those surfaces reject more than a non-finite operand. `pid::compute` also
rejects a step that is not a positive finite duration; it used to answer a
non-positive step with the stored output on the *success* path, which the caller
could not tell apart from a freshly computed command, so a stopped clock read as
a steady loop. And `pid::compute`'s four-argument overload rejects a non-finite
tracking signal before delegating, because that signal is back-assigned into the
integrator once the cycle succeeds.

**Rejected at construction, so the operand never enters the recursion at all.**
`kalman_filter`, `ekf`, `ukf` and `rls` are built through a fallible `create`
which is their only construction path, and which rejects a configuration
carrying a non-finite `Q`, `R`, `x0` or `P0` (`filter_error`), or -- for `rls` --
a forgetting factor outside `(0, 1]`, a non-finite initial covariance, or a
non-positive covariance bound (`rls_error`). `mhe`, `nmhe` and `recursive_arx`
embed one of those types and forward its rejection. This is the earliest point
at which such a mistake can be reported, and reporting it here is what keeps it
from surfacing as a non-finite estimate several steps later, in a place that
names nothing.

The boundary is **finiteness only**. An ill-conditioned but finite
configuration -- a covariance with entries many orders of magnitude apart, a
singular or zero `P0` -- is a legitimately posed problem and is accepted;
conditioning is a numerical-behavior question, which is what this page is about,
and not a domain violation. Symmetry and positive definiteness are not tested
either, because the types do not promise them of the configuration they are
handed.

**Still propagate faithfully,** with no rejection channel:

- **`predict` on every one of those seven types.** Its input is a command the
  caller already owns and the plant already took, so refusing it would leave the
  filter with no propagation for a step that happened. A prediction that poisons
  the carried estimate is not silent: the next `update` rejects, names the cause,
  and latches the type's `health()` query to `non_finite_estimate`.
- **The configuration and reset paths of all ten converted types.** `pid`'s
  `set_params`, `set_integral` and `freeze_integral`, the adaptive controllers'
  construction and `reset` -- none is fallible, so a non-finite gain, initial
  parameter matrix or reference model still enters through them. The first cycle
  that follows rejects and latches the type's `health()` query, which is what
  that query exists for.
- **`l1_controller`'s projection, in one specific case.** The elementwise clamp
  propagates a `NaN` unchanged, but it *replaces* an infinity with the configured
  bound whenever that bound is finite -- so a finite, in-range command can come
  out of an estimate that carries no information, and nothing downstream can
  detect it. That substitution is reported through `health()` as
  `projection_clamped_non_finite` rather than through the failure channel,
  because the cycle did produce the command the algorithm prescribes. See the
  L1 API page for the mechanism.
- **`rls::update`.** It returns a bare `bool` and skips a sample whose
  denominator is non-finite or near zero, so the caller learns that something was
  skipped but never why. Its *construction* is now validated; its per-step
  surface is not yet on the failure channel.
- **`ekf_config::numerical_eps`, and the `Q`/`R`/`P0` of `mekf`,
  `manifold_ukf` and `complementary_filter`.** Those three factories exist and
  validate their initial quaternion, but not their noise configuration; the
  finite-difference step is likewise unchecked. Both are the same defect class
  the four types above just closed, on surfaces this change did not open.
- Every other surface in the library, including `particle_filter::update` and the
  online planners. These are being moved onto the failure channel too, but they
  have not been moved yet, and until they are the propagation contract above is
  what governs them.

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

- **DARE arithmetic range, and what a common rescale carries:** The discrete
  algebraic Riccati equation solver equilibrates the weights before forming the
  symplectic matrix, then establishes two things about the answer before
  reporting success: that it retains more than half of the scalar type's
  significand, and that the closed loop it produces is inside the unit disk. The
  first is a statement about the ANSWER rather than about the residual, and the
  distinction is load-bearing: measured over 2.5 million poses against
  extended-precision solutions by a different algorithm, the answers that keep
  half the significand and the answers that do not are *contiguous* on the
  residual, so no threshold on that quantity separates them. What is compared
  instead is an estimate of the answer's own relative forward error, obtained by
  inverting the residual map's derivative -- the Stein operator built from the
  closed loop -- against the residual, and the margin it is compared against is
  `sqrt(eps)`, which for a radix-2 type is exactly the retention of half the
  fractional significand bits and is therefore derived rather than fitted. It
  returns
  `ctrlpp::expected<dare_result<Scalar, NX>, dare_error>`, so a refusal is
  named rather than represented by a bare empty result:
  `dare_error::non_finite_input` when A, B, Q or R contains NaN or Inf, and
  `dare_error::arithmetic_limit` when finite inputs cannot produce a verified
  stabilizing result at the scalar type's precision.

  Equilibration's claim is that multiplying both weights by a common positive
  factor poses the same arithmetic problem, and the verification has to honor
  that claim rather than re-introduce the range limit the equilibration removed.
  So the verification runs on the equilibrated problem and the claim about the
  caller's own scale is **carried** across the rescale: the equation is
  homogeneous of degree one in `(P, Q, R)` together and the gain is homogeneous
  of degree zero, so the accuracy estimate is homogeneous of degree zero -- the
  residual it inverts and the solution norm it divides by both scale by the same
  factor -- and the stabilizing spectrum and the definiteness transport unchanged
  as well. Only what homogeneity cannot supply is checked directly -- that the
  rescaled solution is finite, and that it is still positive semi-definite at the
  scale actually returned.

  The check is still repeated at the caller's scale wherever its operands
  resolve, with every magnitude formed by dividing out the operand's largest
  entry first, so it neither overflows nor underflows on finite operands. It
  cannot widen what is accepted; where it resolves and disagrees, **it
  declines**, and where it cannot be formed the carried claim stands alone.

  What is carried and what is refused. Every common-scaled pose whose answer the
  type can represent, and whose answer retains more than half the significand, is
  carried. The ceiling is a derived property -- `max_finite /
  max|P_equilibrated|`, the point at which the returned solution itself stops
  being representable -- and not a constant. At the bottom the binding quantity is
  the rescale's own precision: once the product is subnormal the answer loses
  significand, the gain it implies departs from the equilibrated gain, and the
  solver refuses once that departure exceeds the counted-operation margin. In the
  weight-RATIO direction the binding quantity is the accuracy estimate: a state
  weighting whose condition number reaches `1e8`, or an input weight ten decades
  under the state weight, produces an answer that has lost more than half the
  significand, and it is refused rather than returned.

- **CARE arithmetic range, at both ends:** The continuous algebraic Riccati
  solver verifies the solution it extracted -- residual bound and
  closed-loop spectrum -- before reporting success, on every solve. **That is
  true of every selectable method, not only of the default one**, and it is one
  rule reached from three paths rather than three rules that resemble each
  other. It used to be true of the default alone: on a controllable and
  detectable near-axis family the two Schur tags returned an unstable closed
  loop in 150 of 814 and 156 of 839 successes, while the published result
  contract promised a stabilizing matrix without qualifying by method. A Schur
  method that cannot certify what it extracted now declines with
  `care_error::unverified_solution`. Every
  magnitude entering that verification, and entering the iteration's own
  convergence test, is computed in a form that neither overflows nor underflows
  on finite operands: a magnitude is formed by dividing out the operand's
  largest entry first, so a matrix of finite entries whose sum of squares would
  leave the top of the range still yields a finite magnitude, and one whose sum
  of squares would leave the bottom still yields a nonzero one. **A magnitude
  that cannot be resolved declines rather than certifies** -- an infinite
  comparison on both sides is not evidence, and a residual scale of zero would
  turn the acceptance test into `0 <= 0`.

  The bottom of the range is where the default path actually failed, and the
  rest of this entry is about that path specifically. Its
  stable-subspace factorization compares a squared quantity against an absolute
  floor -- the scalar type's smallest normal value -- and the quantity in
  question is the square of the solution's own coupling. A solution smaller than
  the square root of the smallest normal value therefore had the information
  discarded, and the solver returned success carrying `P = 0`. It now rescales
  that factorization by an exact power of two derived from the operand, which
  the factorization is equivariant under, and the band is answered rather than
  lost. Where the rescale does not reach -- the Newton step forms its own
  inverse, whose entries are reciprocal squared magnitudes and reach zero
  independently -- the answer is refused, not returned.

  A solution below the representable band is therefore **declined**, and so is
  one the arithmetic can form but not verify. Note the practical cost: the
  default path performs no weight equilibration, so a common rescale of `Q`
  and `R` far above the dynamics' own scale is declined rather than carried. A
  rescale that does not increase the weights is carried at every magnitude
  swept.

  **That cost is method-specific, and one tag does not pay it.** Swept over
  eighteen decades of a common rescale in both directions, on a comfortably
  damped family and a structurally simple one, `balanced_schur_care_method`
  answers all 1,152 draws of each population and every answer agrees with the
  scale-invariant gain oracle; the default answers 632 and 641. The difference
  is the DGEBAL-style diagonal balance that variant applies before factorizing,
  which is exactly the equilibration the default lacks. A caller whose weights
  sit far from their dynamics' own scale, and who would rather pay for an answer
  than receive a decline, should select it -- see
  [lqr](../../api/control/lqr.md) for what that choice costs in arithmetic.

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
| Constructor/config | Library rejects or falls back on degenerate configs. `kalman_filter`, `ekf`, `ukf` and `rls` reject a non-finite noise or initial-condition field through a fallible `create`, which is their only construction path; finiteness only, so ill-conditioned-but-finite configurations are accepted |
| Algorithm internals | Library uses numerically stable formulations |
| Outputs | NaN/Inf propagates faithfully, never silently clamped, except on the surfaces named above |
| Carried estimator state | The seven `update` surfaces named above reject a non-finite operand before mutating anything |
| Carried controller state | `pid::compute` (both overloads), `mrac_controller::evaluate` and `l1_controller::evaluate` reject before mutating anything, and return the command on the success channel |

The outputs row is still being narrowed for the per-step surfaces under
conversion: where a surface gains a typed failure return, a degenerate step is
reported through that return instead of being left to propagate as NaN. See
[error-reporting.md](error-reporting.md) for the two reporting channels and for
which surfaces have been converted so far.
