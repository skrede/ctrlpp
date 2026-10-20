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
| `particle_filter::update` | `expected<void, pf_update_error>` (particles, weights and random generator unchanged on refusal) |
| `pid::set_params`, `pid::set_integral` | `expected<void, pid_config_error>` |
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

The boundary is **finiteness, symmetry and definiteness**. A `Q`, `R` or `P0`
whose asymmetry exceeds the rounding bound is refused as
`asymmetric_process_noise`, `asymmetric_measurement_noise` or
`asymmetric_initial_covariance`, and one whose asymmetry is within it is stored
as its symmetrized matrix (see
[the drift convention below](#rounding-drift-is-accepted-structure-is-refused)).
A matrix with an eigenvalue below `-floor`, with the floor
`max(N eps max|M_ij|, N denorm_min)` (the symmetric eigensolver's
backward-error order, kept from underflowing to zero), is refused as
`indefinite_process_noise`, `indefinite_measurement_noise` or
`indefinite_initial_covariance`: a matrix that is not symmetric positive
semidefinite is not a covariance, for any filter. The same checks run in `mekf`
and `manifold_ukf` after their initial quaternion, and in
`kalman_filter::set_noise` and `kalman_filter::reset_covariance`, which refuse
before assigning anything. A singular positive semidefinite matrix (a singular
or zero `P0`, a rank-one `Q`, a noiseless `R`) is accepted, and so is an
ill-conditioned one with entries many orders of magnitude apart: conditioning
is a numerical-behavior question, which is what this page is about, and not a
domain violation.

`particle_filter` is built the same way, through `create` or
`make_particle_filter`, and additionally refuses a singular `R`, because it
scores every particle with `R`'s inverse and log-determinant. A singular
positive semidefinite `Q` or `P0` is still accepted. The refusals are listed in
[particle-filter.md](../../api/estimation/particle-filter.md).

**Still propagate faithfully,** with no rejection channel:

- **`predict` on every one of those seven types.** Its input is a command the
  caller already owns and the plant already took, so refusing it would leave the
  filter with no propagation for a step that happened. A prediction that poisons
  the carried estimate is not silent: the next `update` rejects, names the cause,
  and latches the type's `health()` query to `non_finite_estimate`.
- **The configuration and reset paths of all ten converted types.** The adaptive
  controllers' construction and `reset` are not fallible, so a non-finite initial
  parameter matrix or reference model still enters through them. The first cycle
  that follows rejects and latches the type's `health()` query, which is what
  that query exists for. `pid` is the exception: `create`, `set_params` and
  `set_integral` refuse a non-finite or contradictory value with a
  `pid_config_error`, and `freeze_integral` and `reset` take no value that could
  be non-finite. `create` is its only construction path; the constructor is
  private.
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
- **`ekf_config::numerical_eps`, and the scalar options of the quaternion
  filters** (`mekf_config::dt` and `numerical_eps`, the geodesic-mean options of
  `manifold_ukf_config`, and the gains and step of `cf_config`). Their factories
  validate the initial quaternion and the noise and initial covariances, but not
  these scalars; they are the same defect class, on surfaces this change did
  not open.
- Every other surface in the library, including the online planners. These are
  being moved onto the failure channel too, but they have not been moved yet,
  and until they are the propagation contract above is what governs them.

## What the library does guard against

ctrlpp hardens its internals against degenerate-but-valid inputs that would
otherwise produce silent corruption through intermediate overflow:

- **Spline knot spans:** B-spline evaluation treats a knot span as empty only
  when it is exactly zero, the 0/0 := 0 convention of the B-spline recurrences,
  and evaluates every positive span as given. A span short enough that a
  derivative control point overflows is refused by `bspline_trajectory::create`
  as `spline_error::unrepresentable_spline` instead of being evaluated.

- **Thomas solver pivot collapse:** The tridiagonal solver checks every pivot
  once during forward elimination against a counted rounding bound scaled by
  the pivot's row, and refuses a pivot it cannot tell from zero (a NaN pivot
  included) as a typed error instead of returning a solution. The cyclic solver
  applies the same refusal to its Sherman-Morrison denominator, which vanishes
  exactly when the cyclic matrix is singular. `cubic_spline::create` reports
  either refusal as `spline_error::unrepresentable_spline`.

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
  `ctrlpp::expected<dare_result<Scalar, NX, NU>, dare_error>`, so a refusal is
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
  `care_error::unverified_solution`. Both legs of that verification are counted
  rather than tuned, and the residual leg is counted at the basis the extraction
  actually reads: every continuous path takes its subspace basis from a
  `2n`-by-`2n` factorization, the two Schur tags from the Hamiltonian and the
  default tag from the sign-function projector, so the bound counts a
  `2n`-by-`2n` Householder factorization and its solve. Counting it at `n` would
  under-count every path by roughly a factor of eight and refuse answers that are
  correct. The closed-loop leg walks the factor block by block and decides each
  2 by 2 block by its discriminant, because a block holding two real eigenvalues
  has a half-trace that is neither of them and testing the half-trace can certify
  a mode in the open right half-plane. Every
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
  one the arithmetic can form but not verify.

  **`ctrlpp::care` equilibrates both weightings before it builds the
  Hamiltonian**, for every method tag. It divides `Q` and `R` by their largest
  entry, solves at that scale, and multiplies the solution back. The equation is
  homogeneous of degree one in `(P, Q, R)` taken together, so the rescale changes
  nothing about the answer, and the gain `K = R^-1 B' P` does not move at all.
  What it changes is which poses have an answer: the rescale is a block-diagonal
  similarity of the Hamiltonian -- with `D = diag(I, sI)`, `D^-1 H D` is exactly
  the Hamiltonian of the divided pose -- so it cannot move the spectrum, and the
  divisor is chosen to land both off-diagonal blocks at unit order at once.

  Swept over eighteen decades of a common rescale in both directions, on a
  comfortably damped family and a structurally simple one: **all three tags now
  answer all 1,152 draws of each population.** Before the equilibration the
  default answered 652 and 655 and only `balanced_schur_care_method` answered
  every draw, because its DGEBAL-style balance was the only equilibration the
  library had. A caller whose weights sit far from their dynamics' own scale no
  longer needs to select a tag to get an answer -- see
  [lqr](../../api/control/lqr.md) for what the tags cost in arithmetic.

  **The acceptance check's cost is now TIMED, and the counted figure it replaces
  was not a cost.** The operation count put the continuous check at "at most
  8.3 percent" of the solve; timed with the check as the independent variable it
  exceeds that at every `NX` from 4 to 24, reaching 14.9 percent, and it runs from
  0.87x to 2.00x the counted figure depending on the tag and the size. Counts and
  timings are both published in `detail/care_methods.h`, each labeled, because a
  count is an honest statement about arithmetic and a poor predictor of wall time
  at sizes where the whole working set is cache-resident. Budget from the timed
  column. The selection between tags is unaffected: the default pays the larger
  RELATIVE overhead precisely because it has the cheaper solve to pay it out of,
  and it wins on absolute time at every swept size.

  Those counts are **pinned by an assertion**, not merely recorded here, so a
  future tightening of the acceptance rule that costs those answers breaks a
  test rather than leaving this recommendation wrong.

  Two criteria agree that the answers are right, and they are worth
  distinguishing. Each draw's gain is compared against the same pose at unit
  weight scale, which is a **scale-invariance consistency check**: the gain
  `K = R^-1 B' P` does not move under a common rescale, so the two must agree --
  but the reference comes from the same solver under the same tag, so it cannot
  detect an error common to both scales, and a weight-scale defect is exactly
  that error. The structurally simple family additionally admits a **closed
  form**: for `A = -I`, `B = I`, `Q = R = sI` the equation collapses to
  `p^2 + 2sp - s^2 = 0`, giving `P = s(sqrt(2) - 1) I` and a gain of
  `(sqrt(2) - 1) I` at every `s`. That criterion owes the solver nothing. Zero
  disagreements under either, for all three tags, with a worst relative error
  against the closed form of `4.02e-16`.

  **This does not mean the balanced variant answers everything.** It answers the
  poses a bad weight scale would otherwise lose. On the near-axis band -- where a
  plant eigenvalue is driven toward the imaginary axis and the stable and
  unstable invariant subspaces stop being separated -- the same variant declines
  every one of 112 poses, and no balance repairs that, because the obstacle is
  the pose rather than its scaling.

- **L1 DC gain inversion:** The L1 adaptive controller validates that the
  predictor model's DC gain is invertible before computing the feedforward gain
  `K_r`. Construction goes through the fallible factory
  `l1_controller::create`, which returns
  `ctrlpp::expected<l1_controller, l1_error>` and reports
  `l1_error::singular_predictor` when (I - A_m) is singular,
  `l1_error::singular_dc_gain` when the DC gain (I - A_m)^{-1} B is singular or
  non-finite, and `l1_error::non_finite_gain` when `K_r` is non-finite. It does
  not throw.

- **Particle filter weight depletion:** When no particle explains a finite
  measurement (every likelihood is zero), the normalizer resets to uniform
  weights rather than producing NaN from `-inf - (-inf)` or dividing by a zero
  total. The step succeeds and latches `pf_health::weights_reset`, so the
  depletion stays distinguishable from a refused measurement. A non-finite
  measurement never reaches this guard: `update` refuses it with
  `pf_update_error::non_finite_measurement` and leaves the particles, the weights
  and the random generator unchanged.

- **UKF covariance update:** The unscented Kalman filter uses the algebraically
  complete minimum mean-square-error reduction `P = P - K*S*K^T`, where
  `S = Pzz + R` and `K = Pxz*S^{-1}`. Because `K*S*K^T = K*Pxz^T`, this term is
  exactly the uncertainty the measurement removes and no extra `K*R*K^T` term is
  added. The reduction is taken from the covariance the sigma set represents,
  which is the matrix `S S^T` of the square root the points were built from.
  When the carried covariance was not positive definite and had to be repaired
  to the nearest symmetric positive definite matrix, that is the repaired
  matrix, so the posterior is the repaired covariance less the gain term and is
  positive semidefinite up to rounding. The repair is reported through
  `health()` as `covariance_repaired`. The manifold UKF follows the same rule in
  its tangent space and reports the repair through its own `health()`.

- **MOESP degenerate data:** The subspace identification algorithm refuses
  with `sysid_error::model_not_identified` when the condition number of the
  projection is not finite (its `NX`-th singular value is zero, as for an
  all-zero record) or when an extracted system matrix is non-finite. It
  returns no sentinel model, so a returned model always has a finite
  condition number.

- **RLS overflow guard:** The recursive least squares estimator skips the
  update step when the quadratic form `phi^T * P * phi` overflows, preventing
  NaN from `0 * Inf` in the gain computation. The `update()` method returns
  `false` when this occurs so the caller can detect skipped updates.

## Rounding drift is accepted, structure is refused

An input that must carry a structural property, such as the symmetry of a
covariance or the unit norm of a rotation quaternion, is judged by its defect:
how far it is from having the property.
The defect is compared with a counted rounding bound `c * N * eps * scale`,
never below the underflow term `N * denorm_min`.

- A defect within the bound is **drift**, the residue of the rounded arithmetic
  that formed the input. It is accepted and projected onto the property, and
  the projected input is what the library stores and uses. An input that
  already has the property is left bitwise untouched.
- A defect beyond the bound is **structure**: the input is a different object
  from the one the property describes. It is refused with a typed error that
  names the input.

`N` is a count in machine epsilons, `eps = 2 u`, derived for each application
from the convention that an `n`-term inner product carries at most about `n u`
of its scale, with `u` the unit roundoff (Higham, *Accuracy and Stability of
Numerical Algorithms*, 2nd ed., 2002, Sec. 3.1). The counted factor `c` is the
order-one prefactor of that count. It
is a defaulted parameter of the shared helpers, `symmetrize_within_drift` and
`psd_decomposition` in `ctrlpp/detail/covariance_ops.h` and
`checked_unit_quaternion` in `ctrlpp/estimation/detail/filter_configuration.h`,
default 1, and not a field of any configuration. The absolute term is the gradual-underflow part of the
standard model, `fl(x op y) = (x op y)(1 + delta) + eta` with
`|eta| <= denorm_min / 2` (Higham, Sec. 2.1). Without it the bound would
collapse to zero for a subnormal-scale input, which would then be judged with
no tolerance at all; at and above `scale = denorm_min / (c eps)` the bound is
the relative term exactly, so no normal-scale verdict depends on it.

### Covariances

A covariance is symmetric positive semidefinite by definition. Each `Q`, `R` and
`P0` is checked in this order:

1. A non-finite entry is refused (`non_finite_process_noise`,
   `non_finite_measurement_noise`, `non_finite_initial_covariance`).
2. The defect is `max|M_ij - M_ji|` and the scale is `max|M_ij|`. Each triangle
   of a matrix formed by products carries at most `N u` of that scale, so the
   two triangles of one exactly symmetric matrix differ by at most
   `2 N u = N eps` times it: that is the count behind `c = 1`. Within the bound,
   each differing pair is replaced in both positions by its midpoint
   `M_ij + (M_ji - M_ij) / 2`, which cannot overflow, and an equal pair is never
   written. Beyond it the matrix is refused (`asymmetric_process_noise`,
   `asymmetric_measurement_noise`, `asymmetric_initial_covariance`).
3. The definiteness verdict runs on the symmetrized matrix. An eigenvalue below
   `-floor`, with `floor = max(c N eps max|M_ij|, N denorm_min)`, is refused
   (`indefinite_process_noise`, `indefinite_measurement_noise`,
   `indefinite_initial_covariance`). The relative term is the order of the
   symmetric eigensolver's backward error (the Weyl perturbation bound). The
   absolute term exists because the relative one underflows to zero below
   `max|M_ij| = denorm_min / (c eps)`, where a positive semidefinite matrix
   whose rounded entries put an eigenvalue one subnormal unit below zero would
   otherwise be refused.

Across matrices the order is the order the configuration declares its fields;
the particle filter's singular-`R` check follows `R`'s definiteness.

This applies today to `Q`, `R` and `P0` at every filter factory
(`kalman_filter`, `ekf`, `ukf`, `mekf`, `manifold_ukf` and `particle_filter`,
whose sampling factors, inverse and log-determinant are formed from the
symmetrized matrices) and at `kalman_filter::set_noise` and
`kalman_filter::reset_covariance`. It does not yet apply to the weight matrices
of the regulators, the Riccati solvers and the predictive controllers, nor to
the moving-horizon estimators' own copies of `Q`, `R` and `P0`: `mhe` and
`nmhe` forward the embedded filter's refusal, but the copies they keep for
their own cost are stored as given.

### Unit quaternions

A rotation quaternion is unit by definition. The initial attitude `q0` of
`mekf`, `manifold_ukf` and `complementary_filter` is checked before any other
configuration field:

1. A zero or non-finite norm is refused (`degenerate_quaternion`): normalizing
   it would produce NaN.
2. The defect is `| ||q0|| - 1 |`, the scale is 1 and `N = 4`, so the bound is
   `max(4 c eps, 4 denorm_min)`, which is `4 c eps` for any `c` of order one.
   The count behind `c = 1`: a quaternion normalized once in working precision
   has a norm defect of at most `(N / 2 + 2) u` (the four-term sum of squares,
   the square root and the division), and the check's own evaluation of the
   norm adds at most `(N / 2 + 1) u`; their sum, `7 u`, is within
   `N eps = 8 u`. Beyond the bound `q0` is refused (`non_unit_quaternion`): a
   quaternion of norm 2 is not an inaccurate attitude but no attitude at all.
3. Within the bound the drift is projected away: the filter stores
   `q0.normalized()`. A `q0` whose computed norm is exactly one is divided by
   one and so stored bitwise as given.

The unit tests sweep quaternions from the producers an application uses, at
both precisions: normalized Gaussian 4-vectors, angle-axis conversions,
`so3::exp` of random rotation vectors, yaw-pitch-roll compositions of three
angle-axis rotations, and a double unit quaternion rounded to float. Every draw
is accepted at `c = 1`.

### Rotation matrices and rigid transforms

No ctrlpp input takes a rotation matrix or an element of SE(n) today. One that
does follows the same rule. The defects are the orthogonality defect
`||R^T R - I||` and the orientation defect `|det R - 1|`, each against a bound
of the same form `max(c N eps, N denorm_min)` with `N` counted from the
arithmetic that forms it (an entry of `R^T R` is an `n`-term inner product). Within
the bound the input is projected onto the group: `R` is replaced by its
orthogonal polar factor, the nearest orthogonal matrix in the Frobenius norm
(Higham, "Computing the Polar Decomposition with Applications", SIAM J. Sci.
Stat. Comput. 7(4), 1986), and the translation of a rigid transform is kept.
Anything else is refused with a typed error; a reflection, `det R = -1`, is
structure, not drift.

## Summary

| Layer | Behavior |
|---|---|
| Inputs (your code) | Validate at system boundaries |
| Constructor/config | Library rejects or falls back on degenerate configs. `kalman_filter`, `ekf`, `ukf` and `rls` reject a non-finite noise or initial-condition field through a fallible `create`, which is their only construction path; `kalman_filter`, `ekf`, `ukf`, `mekf`, `manifold_ukf` and `particle_filter` also refuse a structurally asymmetric or indefinite `Q`, `R` or `P0` and store a drift-asymmetric one symmetrized ([the drift convention](#rounding-drift-is-accepted-structure-is-refused)), and `kalman_filter::set_noise` and `reset_covariance` do the same before assigning. `mekf`, `manifold_ukf` and `complementary_filter` refuse an initial quaternion whose norm is not one within rounding and normalize one whose norm is. A singular positive semidefinite or an ill-conditioned configuration is accepted. `particle_filter` (through `create` or `make_particle_filter`, likewise its only path) additionally refuses a singular `R`, because it scores with `R`'s inverse and log-determinant; a singular positive semidefinite `Q` or `P0` is accepted and sampled without lifting its zero eigenvalues |
| Algorithm internals | Library uses numerically stable formulations |
| Outputs | NaN/Inf propagates faithfully, never silently clamped, except on the surfaces named above |
| Carried estimator state | The seven `update` surfaces named above reject a non-finite operand before mutating anything |
| Carried controller state | `pid::compute` (both overloads), `mrac_controller::evaluate` and `l1_controller::evaluate` reject before mutating anything, and return the command on the success channel |

The outputs row is still being narrowed for the per-step surfaces under
conversion: where a surface gains a typed failure return, a degenerate step is
reported through that return instead of being left to propagate as NaN. See
[error-reporting.md](error-reporting.md) for the two reporting channels and for
which surfaces have been converted so far.
