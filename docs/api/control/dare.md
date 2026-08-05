# dare

Discrete Algebraic Riccati Equation solver using symplectic Schur decomposition. Finds the stabilizing solution P to A'PA - P - A'PB(R + B'PB)^{-1}B'PA + Q = 0. This is the workhorse behind `lqr_gain` and the terminal cost computation in MPC.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::dare<Scalar, NX, NU>` | `#include <ctrlpp/control/dare.h>` |
| `ctrlpp::dare<Scalar, NX, NU>` | `#include <ctrlpp/dare.h>` (convenience) |
| `ctrlpp::dare_error`, `ctrlpp::dare_result` | `#include <ctrlpp/control/dare_types.h>` |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
| `NX` | `std::size_t` | State dimension |
| `NU` | `std::size_t` | Input dimension |

## Functions

### dare

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
auto dare(const Matrix<Scalar, NX, NX>& A,
          const Matrix<Scalar, NX, NU>& B,
          const Matrix<Scalar, NX, NX>& Q,
          const Matrix<Scalar, NU, NU>& R)
    -> ctrlpp::expected<dare_result<Scalar, NX, NU>, dare_error>;
```

Solves the standard DARE. Before forming the 2n x 2n symplectic matrix, it
applies a common positive divisor to Q and R. This does not change the optimal
gain: the equilibrated solution is multiplied by the divisor before it is
returned. The divisor is the largest-magnitude entry across both weights, so
one weight entry is near unity and independently common-scaled poses reduce to
the same canonical problem.

The divisor is chosen **before** anything is factorized, and the symplectic
operands are then factorized once, from the equilibrated weighting. That order
matters for more than cost. The input Gramian `G = B R^{-1} B'` is the quantity
equilibration exists to keep representable: on a weighting whose entries are
subnormal, `R^{-1}` overflows and `G` is infinite at the caller's scale while
being an ordinary number at the equilibrated one. Forming it at the caller's
scale first and correcting afterwards cannot recover those poses, because the
correction would be applied to an infinity.

The solver computes a real Schur decomposition, reorders eigenvalues inside the
unit disk to the top-left block, and extracts P = U21 * U11^{-1}. Before
reporting success it establishes two things about the answer: that it retains
more than half of the scalar type's significand, and that the closed-loop
spectrum it produces is inside the unit disk by more than its own backward-error
margin.

### What a success claims, and how it is decided

**A success claims accuracy of the ANSWER, not smallness of the residual.**
The quantity compared against a margin is an estimate of `P`'s own relative
forward error, and the margin is `sqrt(eps)`.

The estimate comes from the residual map's derivative rather than from the
residual itself. Writing `A_cl = A - BK` for the closed loop and
`Omega(X) = X - A_cl' X A_cl` for the Stein operator built from it, the residual
of a computed solution satisfies `R(P_hat) = -Omega(P_hat - P) + O(||P_hat -
P||^2)` to first order, so inverting `Omega` against the residual estimates the
error in the returned matrix directly:

```
estimated relative forward error  =  ||Omega^-1(R(P_hat))||_F / ||P||_F
```

That is solved on the symmetric subspace, whose dimension is `M = NX(NX+1)/2`,
from quantities the solve already holds. It allocates nothing. Assembling the
operator costs exactly `NX^4` operations -- `NX` diagonal columns at `NX^2` each
plus `NX(NX-1)/2` off-diagonal columns at `2 NX^2` each, since an off-diagonal
column carries a second outer product -- and the rank-revealing solve is
`(2/3) M^3`, with the coordinate round-trip `O(M^2)`. Against a seven-iteration
Schur solve that is 3.6% at `NX = 2`, 10.7% at `NX = 4`, 24.5% at `NX = 6` and
47.7% at `NX = 8`.

Storage for the operator is `M x M`, so it grows as the fourth power of the
state dimension. The decomposition factorizes into the operator's own array
rather than a copy of it, which removes the second live `M x M` array: the
estimator's stack frame is 864, 2,400, 6,432 and 15,056 bytes at
`NX = 2, 4, 6, 8`, and the peak of the whole `dare` call chain -- the quantity a
hard-real-time caller budgets -- is 5,352, 11,688, 23,144 and 42,760 bytes for an
`NX`-state, 3-input pose. The same chain with no accuracy estimate on it at all
still peaks at 17,112 bytes at `NX = 6` and 28,168 at `NX = 8`, so on a 4-16 KB
task stack the supported maximum is `NX = 4` whether the estimate is formed or
not; the estimator is not what decides the small-stack answer.

**`sqrt(eps)` is derived and not calibrated.** For a radix-2 type with
`eps = 2^-p`, `sqrt(eps) = 2^(-p/2)` is exactly the retention of `p/2` of the `p`
fractional significand bits, so "more than half the significand of the returned
answer is correct" *is* "relative forward error at most `sqrt(eps)`", in the
type's own radix. It carries to `float` and to any other radix-2 type without
being re-measured.

**Why the residual is not the quantity compared.** Measured over 2.5 million
poses against extended-precision solutions computed by a different algorithm, the
answers that keep more than half the significand and the answers that do not are
*contiguous* on the residual: the worst kept and the best lost are adjacent, and
the distribution is unimodal across ten decades with no gap anywhere in it. No
threshold on that quantity separates them. The margin this replaced -- a
counted-operation envelope on the residual -- refused none of the answers in that
population that had lost more than half their significand.

**The accepted set moved, and it moved inward.** Poses whose answer loses more
than half the significand are now refused with `arithmetic_limit` where they
previously returned a value. On the two recorded weight-ratio populations this
removes 531 of 13,114 accepted rows, and every one of the 531 is confirmed by an
extended-precision reference to have lost more than half the significand; no row
that the reference calls correct stopped being accepted, and no row that was
refused became accepted. A caller who was relying on a returned `P` at a weight
ratio around ten decades, or on a state weighting whose condition number is
`1e8` or worse, will now see `arithmetic_limit` instead of an answer that was
wrong in its eighth digit.

**What the estimate is not.** It is a first-order estimate, not a bound. Its
left tail under-predicts where cancellation in the residual happens to be
favorable, so a small population of answers that have genuinely lost more than
half the significand is still returned rather than refused. Closing that gap
would need a second solve in higher precision, which is not a postcondition.
This is a known and accepted exposure, not an omission.

**When the operator is too ill-conditioned to invert, that is an absence of
evidence.** `Omega` is exactly singular only when the closed loop carries an
eigenvalue pair with `lambda_i * lambda_j = 1`, which a spectrum strictly inside
the unit disk forbids. The solvability test actually performed is a numerical
rank test at a threshold relative to the largest pivot, so it reports rank
deficiency whenever `cond(Omega)` exceeds `1 / (M * eps)` -- roughly
`1 - |lambda|^2 < M * eps` -- and the stabilizing claim is established against a
much wider margin. There is therefore a band in which the spectrum check passes,
the answer may be perfectly good, and the estimate still cannot be formed. In
that band the estimate reports that it has no opinion rather than refuting the
answer, and the closed-loop spectrum check remains the sole owner of the
stabilizing claim. The refusal a caller sees is unchanged -- an estimate that
cannot be formed still declines, under the "a magnitude the verification needs
could not be formed" clause of `arithmetic_limit` rather than the
"forward error exceeded `sqrt(eps)`" clause. Measured over 1.5 million poses
outside the randomized target's conditioning filters, this fires on 115 of them
and on none of the 971,244 poses inside the filters.

`@cite laub1979`, and Higham, *Accuracy and Stability of Numerical Algorithms*,
2nd ed., Ch. 19, for the orthogonal-transformation error analysis the backward
error rests on.

### Where the claim about your own scale is made

Those checks run on the **equilibrated** problem, and the claim about the
problem you posed is then *carried* across the rescale rather than re-formed at
your scale. The equation is homogeneous of degree one in `(P, Q, R)` taken
together, so multiplying all three by a positive factor multiplies the residual
by exactly that factor; the gain `(R + B'PB)^{-1} B'PA` is homogeneous of degree
zero, so the closed loop is unchanged; and a positive factor cannot move an
eigenvalue across zero. The accuracy estimate is therefore homogeneous of degree
**zero**: the residual it inverts scales by the factor, the operator it inverts
does not, and the `||P||_F` it divides by scales by the factor as well, so the
estimated relative forward error is the same number at both scales. What
homogeneity does not supply is checked directly: the rescaled solution must be
finite, and it is re-tested for positive semi-definiteness at the scale actually
returned.

Re-forming the evidence at your scale instead of carrying it is what limited the
accepted range before. Every magnitude the verification forms is now computed by
dividing out the operand's largest entry first, so on a pose whose *answer* is an
ordinary normal number the *evidence* no longer leaves the range -- a comparison
between two infinities was refusing problems the solver could solve.

The direct check at your scale is still performed, and every magnitude it forms
is computed by dividing out the operand's largest entry first, so it neither
overflows nor underflows on finite operands. It cannot widen what is accepted,
but it can refuse: where the carried claim and the direct check disagree, **the
direct check declines**.

Where the direct check cannot be formed, what happens depends on **why**, and the
two cases are deliberately not merged:

- **The gain could not be formed at your scale**, because `R + B'PB` left the
  range there. That is a statement about the scale rather than about the answer,
  it is exactly the band the carried claim was written to cover, and the answers
  in it are bit-exact (see below). The carried claim stands and the answer is
  returned.
- **Anything else failed to resolve**, including the forward-error estimator
  declining. These come with no account of why, and the carried claim's premise
  is that both scales pose the same problem -- an unexplained dissolution of the
  direct check is the weakest place to assume it. **The solver declines.**

Merging the two would force one disposition on both: either discard a band of
exactly-correct answers, or accept an absence of evidence whose cause is unknown.

### The accepted range's ceiling

The ceiling is a derived property of the scalar type's range, not a constant.
The returned matrix is the equilibrated solution times the common divisor, so
the largest divisor whose answer is representable is

```
max_finite / max|P_equilibrated|
```

The quantity that limits the range is therefore the returned solution's own
largest entry. On the equal-weight scalar pose `A = 0.5`, `B = 1`, `Q = R = c`
that derivation puts the ceiling at `c = 1.586972e+308`, which is where a
bisected sweep measures the last carried magnitude. Above it the answer is not
representable and the refusal is `arithmetic_limit`.

The direct check's own reach stops slightly earlier -- 0.2748 decades earlier on
that pose -- because the gain it recomputes needs the sum `R + B'PB` formed at
your scale, and that sum leaves the range before the answer does. That is
exactly the band the carried claim covers, and the band is not a marginal one:
across it every returned answer reproduces the homogeneous truth `c * P_unit` to
**zero** relative error. Measured, a rule that declined it instead would drop the
ceiling from `1.586972e+308` to `8.428864e+307` and give up bit-exact answers for
nothing. `dare_hardening_test` pins this band, checking both that it is accepted
and that its lower edge is found by walking until the direct check resolves again
rather than by a hardcoded bound.

At the bottom the binding quantity is different again, and it is the rescale
itself: multiplying the equilibrated solution by the divisor keeps full relative
precision only while the product is normal, and once it is subnormal the answer
loses significand and the gain it implies departs from the equilibrated gain.
The solver refuses when that departure exceeds the counted-operation margin,
which on the same pose is measured at a common scale near `1e-316`.

On success `result->P` is the stabilizing solution and `result->K` is the
feedback gain `(R + B'PB)^{-1} B'PA` that the solve itself formed while verifying
that solution; `result->subspace_separation` and `result->reorder_complete` are
conditioning diagnostics.

`result->K` is the gain the acceptance decision was made on, and it is not
equivalent to re-forming the gain from `result->P`. The solve forms it at the
equilibrated scale and applies the same dimension-times-unit-roundoff rank test
to `R + B'PB` that the pencil build applies to `R`, then checks the result is
finite; a caller who re-forms it from `P` alone reproduces neither. The gain is
homogeneous of degree zero in `(P, Q, R)`, so the equilibrated gain **is** the
caller's gain -- but the caller-scale sum `R + B'PB` leaves the top of the range
on poses whose answer is an ordinary normal number, and a rank-revealing solve
handed such a sum returns a **zero** gain rather than refusing. Measured over the
scalar weight-ratio sweep, 3,546 faithfully realized accepted poses spanning
every decade the type supports at three weight ratios: re-forming the gain from
`P` reaches a relative error of `1.000` on two of them -- a returned gain of
zero, which is no feedback at all -- while `result->K` stays within
`1.547e-16`, below one unit of roundoff, on every one.

Refusals:

| Enumerator | Condition |
| --- | --- |
| `dare_error::non_stabilizable` | fewer than n eigenvalues of the symplectic spectrum lie inside the unit region |
| `dare_error::non_finite_input` | A, B, Q or R contains NaN/Inf |
| `dare_error::singular_a` | A is rank-deficient to a scale-relative reciprocal-pivot tolerance, so the `A^{-T}` the pencil build needs does not exist |
| `dare_error::singular_r` | R is rank-deficient to the same tolerance, so the `R^{-1}` the same pencil build needs for `G = B R^{-1} B'` does not exist. The test is applied to the **equilibrated** weighting, the operand the solve actually inverts, so it also covers a weighting that is nonzero on its own but vanishes against the state weighting -- `Q = 1e300 I` with `R = 1e-300` divides to a weighting that underflows to zero, and relative to the posed problem the input weighting is zero |
| `dare_error::singular_u11` | the top-left block of the reordered invariant-subspace basis is singular; P cannot be extracted. **Covers two different situations -- see below** |
| `dare_error::non_psd_solution` | the extracted P is not positive semi-definite. The test is an LDLT pivot-sign test against `N * eps * max\|P_ij\|`, the order of the factorization's own backward error, below which a pivot carries no sign information |
| `dare_error::schur_failed` | the real Schur factorization did not converge |
| `dare_error::arithmetic_limit` | finite inputs could not produce a verified stabilizing solution at the scalar type's precision. Specifically: the equilibrated answer's estimated relative forward error exceeded `sqrt(eps)`, so it retains less than half the significand; its closed-loop spectrum or gain solve did not verify; a magnitude the equilibrated verification needs could not be formed; the weights overflowed while being equilibrated; the solution overflowed while being rescaled; the rescaled solution failed the positive-semi-definiteness test at its returned scale; the direct check at the caller's scale resolved and refuted the carried claim; the direct check failed to resolve for any reason other than the gain leaving the range at the caller's scale; or the equilibrated and returned gains, both formed, disagreed by more than the counted-operation margin |

A second cause is no longer on that list either: a weighting that vanishes under
equilibration. That is now `singular_r`, which names the weight ratio the caller
chose rather than sending them to look at precision. Measured across sixteen
extreme ratio poses spanning `Q` from `1e150` to `1e308` against `R` from
`1e-150` to `1e-320`, fifteen moved from `arithmetic_limit` to `singular_r`; no
pose moved between accepted and declined.

The one cause that is **no longer** on that list is the gain the check at the
caller's scale could not form because `R + B'PB` left the range there. That used
to be a refusal; it is now an absence of evidence with a known cause, and the
carried claim is what decides. Every *other* magnitude the direct check fails to
form is still a refusal.

Positive semi-definiteness is tested on both matrices that exist, and
deliberately so. The extraction primitive tests the P the solve produced; when
common-weight equilibration is active that P is the *equilibrated* one, and the
value handed back to the caller is `P * weight_scale`, which the first test never
saw. The rescale is therefore re-tested at the returned scale. Scaling by a
positive factor preserves definiteness mathematically, so this second test only
ever fires on the rounding the rescale itself introduces -- which is why it
reports `arithmetic_limit` rather than `non_psd_solution`, and why the floor has
to carry the factorization's backward error rather than assume there is none.

An indefinite state weighting is not a supported input. `Q = -I` on an otherwise
well-posed pair leaves fewer than n eigenvalues of the symplectic spectrum inside
the unit disk, so it is refused as `non_stabilizable`.

### What `singular_u11` covers, and what `non_stabilizable` misses

These two are worth reading together, because the placement count cannot see every pair that has no stabilizing solution.

An uncontrollable mode at `|lambda| > 1` contributes **both** `lambda` and its reciprocal to the symplectic spectrum, so n eigenvalues do lie inside the unit disk and `non_stabilizable` does not fire; the invariant subspace fails to project instead, and the refusal arrives as `singular_u11`. An uncontrollable mode at `|lambda| = 1` contributes two eigenvalues **on** the circle, neither inside, so that case does reach `non_stabilizable`. Both are refusals -- **no gain is ever returned for a pair with no stabilizing solution** -- but only the second names the structural cause.

The implication that makes `singular_u11` informative is standard and exact: a stabilizable and detectable pair has a nonsingular `U11` (Laub 1979 Sec. III), so **in exact arithmetic** a singular `U11` implies the pair is not both stabilizable and detectable.

**That qualifier is load-bearing, which is why the enumerator is not named after the structural cause.** The test performed is a numerical rank test with a threshold relative to the largest pivot, so a genuinely well-posed pair whose invariant subspace is severely ill-conditioned reaches the same branch. Measured on the one-parameter family `A = diag(2, 1/2)`, `B = [delta; 1]`, `Q = I`, `R = 1`, which is controllable and detectable for **every** `delta != 0`:

| `delta` | `rank[B AB]` | result |
| --- | --- | --- |
| `1e-1` … `1e-3` | 2 | value, `\|\|P\|\|` rising from `8.9e2` to `8.9e6` |
| `1e-4` … `1e-7` | 2 | **`arithmetic_limit`** -- a well-posed pair refused for accuracy |
| `1e-8` … `1e-14` | 2 | **`singular_u11`** -- the same pair refused for rank |
| `1e-16`, `0` | 1 | `singular_u11` -- genuinely uncontrollable |

There are two transitions and they are three and a half decades apart. `\|\|P\|\|` grows as `K / delta^2` with a measured `K = 8.87`, so the estimated relative forward error `\|\|P\|\| * eps` reaches `sqrt(eps)` -- half the significand, the one rule the accuracy gate applies -- at `delta = sqrt(K * sqrt(eps)) = 3.6e-4`, and the pair is declined for accuracy from there down. The smallest pivot of the computed `U11` falls as `C * delta^2` with a measured `C = 0.170`, so it reaches the rank test's `n * eps` threshold at `delta = sqrt(n * eps / C) = 5.1e-8`, and the refusal changes name there. Naming this enumerator after the structural cause would state something false about every row in the middle two bands.

So `singular_u11` covers a pair with no stabilizing solution **and** a numerical failure to separate the invariant subspace, and **it does not distinguish them**. Telling them apart needs a stabilizability test the solver does not perform.

`care_error::singular_u11` has the identical shape for the identical reason, with the unit disk replaced by the open left half-plane and the reciprocal by the reflection `-lambda`.

#### How far the enumerator can be relied on

**Branch on `has_value()`. Treat the enumerator as a diagnostic for a human reader rather than as a control signal.**

The tables above say which condition produces which enumerator, and they hold wherever the condition itself is decidable in the scalar type. They stop holding below the precision's resolution boundary -- the point at which the quantity deciding the problem falls under `sqrt(epsilon)` relative to the operands carrying it, which for the family above is where the two transitions in the previous section sit. Below the boundary the solve is refused, and *which* refusal it carries is decided by instruction selection and by the linear-algebra implementation rather than by the input.

Both transitions were located to one unit in the last place per configuration, and the configurations disagree. Over four compilers (g++ 16.1.1, clang 22.1.8, 21.1.8 and 18.1.8), four optimization levels, fused multiply-add off and on, and two patch releases of Eigen:

| Boundary | Located between | Band width |
| --- | --- | --- |
| accepted to `arithmetic_limit` | `3.279e-04` and `3.738e-04` | 0.057 decades, 14% of the value |
| `arithmetic_limit` to `singular_u11` | `4.740e-08` and `5.741e-08` | 0.084 decades, 21% of the value |

Outside that band every configuration agrees, including at both ends of the family. Inside it, the answer belongs to the build. Measured on one bit-identical input inside the band -- the pair `A = diag(2, 1/2)`, `B = [0; 1]`, `Q = I`, `R = 1` with the cross weight `N = [0.1; 0.2]`, whose unstable mode is exactly uncontrollable so that no stabilizing solution exists:

| Configuration | with Eigen 3.4.0 | with Eigen 3.4.1 |
| --- | --- | --- |
| every compiler and level, no fused multiply-add (16 of 32) | `non_psd_solution` | `singular_u11` |
| every compiler at `-O0`, fused multiply-add enabled (4) | `non_psd_solution` | `singular_u11` |
| g++ at `-O1`, fused multiply-add enabled (1) | `non_psd_solution` | `singular_u11` |
| g++ at `-O2`, fused multiply-add enabled (1) | `singular_u11` | `non_psd_solution` |
| g++ at `-O3`, fused multiply-add enabled (1) | `singular_u11` | `singular_u11` |
| clang at `-O1` and above, fused multiply-add enabled (9) | `singular_u11` | `arithmetic_limit` |

Three enumerators over sixty-four configurations on one input -- `singular_u11` on 33, `non_psd_solution` on 22, `arithmetic_limit` on 9 -- and the majority verdict flips with a patch release of a dependency.

**Every configuration in that table is x86-64, and a green column is not evidence of stability on arm64.** The continuous solver's twin of this measurement produced a third enumerator on arm64 that neither x86-64 compiler produced. The axis that moves the discrete answer is fused multiply-add contraction, which is enabled by default on that architecture, so the contracting rows are the closest available proxy -- a prediction, not a measurement, and not a substitute for running it.

Three enumerators report a property of the operands rather than of an iteration and stay meaningful at any distance from the boundary: `non_finite_input`, `singular_a` and `singular_r`. The first is an exact predicate on the entries; the other two are rank tests relative to the operand's own largest pivot, so they are decided by the input wherever that operand is not itself near rank deficiency, and on an exactly rank-deficient operand the threshold is exactly zero and the verdict is exact. A finite, well-formed input with a nonsingular `A` and a nonsingular `R` never carries any of the three.

### dare (with cross-weight)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
auto dare(const Matrix<Scalar, NX, NX>& A,
          const Matrix<Scalar, NX, NU>& B,
          const Matrix<Scalar, NX, NX>& Q,
          const Matrix<Scalar, NU, NU>& R,
          const Matrix<Scalar, NX, NU>& N)
    -> ctrlpp::expected<dare_result<Scalar, NX, NU>, dare_error>;
```

DARE with state-input cross-weight N. Transforms to standard form via Q' = Q - NR^{-1}N', A' = A - BR^{-1}N' and delegates to the standard solver.

## Usage Example

```cpp
#include <ctrlpp/control/dare.h>

#include <Eigen/Dense>

#include <iostream>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;

    Eigen::Matrix2d A;
    A << 1.0, 0.1, 0.0, 1.0;

    Eigen::Matrix<double, 2, 1> B;
    B << 0.005, 0.1;

    Eigen::Matrix2d Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 1> R;
    R << 1.0;

    auto P_opt = ctrlpp::dare<double, NX, NU>(A, B, Q, R);
    if (!P_opt.has_value()) {
        // P_opt.error() names why the problem was refused.
        std::cerr << "the Riccati solve refused the problem\n";
        return 1;
    }

    std::cout << "P =\n" << P_opt->P << "\n";
}
```

## See Also

- [lqr](lqr.md)<br/> uses DARE internally to compute optimal gains
- [place](place.md)<br/> pole placement as an alternative design method
- [mpc](../mpc/mpc.md)<br/> MPC uses DARE for terminal cost computation
