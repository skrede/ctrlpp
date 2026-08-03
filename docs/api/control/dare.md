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
    -> ctrlpp::expected<dare_result<Scalar, NX>, dare_error>;
```

Solves the standard DARE. Before forming the 2n x 2n symplectic matrix, it
applies a common positive divisor to Q and R. This does not change the optimal
gain: the equilibrated solution is multiplied by the divisor before it is
returned. The divisor is the largest-magnitude entry across both weights, so
one weight entry is near unity and independently common-scaled poses reduce to
the same canonical problem.

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

That is solved on the symmetric subspace, whose dimension is `NX(NX+1)/2`, from
quantities the solve already holds. It allocates nothing. Its cost is
`NX(NX+1)/2` rank-two updates to assemble the operator plus a rank-revealing
solve of about `NX^6 / 12` operations, which against a seven-iteration Schur
solve is 0.4% at `NX = 2`, 3.6% at `NX = 4` and 28.6% at `NX = 8`.

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

The direct check at your scale is still performed **wherever its operands
resolve**, and every magnitude it forms is computed by dividing out the
operand's largest entry first, so it neither overflows nor underflows on finite
operands. It cannot widen what is accepted, but it can refuse: where the carried
claim and the direct check disagree, **the direct check declines**. Where the
direct check cannot be formed at all it has produced no evidence, and the
carried claim stands alone.

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
exactly the band the carried claim covers.

At the bottom the binding quantity is different again, and it is the rescale
itself: multiplying the equilibrated solution by the divisor keeps full relative
precision only while the product is normal, and once it is subnormal the answer
loses significand and the gain it implies departs from the equilibrated gain.
The solver refuses when that departure exceeds the counted-operation margin,
which on the same pose is measured at a common scale near `1e-316`.

On success `result->P` is the stabilizing solution;
`result->subspace_separation` and `result->reorder_complete` are conditioning
diagnostics.

Refusals:

| Enumerator | Condition |
| --- | --- |
| `dare_error::non_stabilizable` | fewer than n eigenvalues of the symplectic spectrum lie inside the unit region |
| `dare_error::non_finite_input` | A, B, Q or R contains NaN/Inf |
| `dare_error::singular_a` | A is rank-deficient to a scale-relative reciprocal-pivot tolerance, so the `A^{-T}` the pencil build needs does not exist |
| `dare_error::singular_r` | R is rank-deficient to the same tolerance, so the `R^{-1}` the same pencil build needs for `G = B R^{-1} B'` does not exist |
| `dare_error::singular_u11` | the top-left block of the reordered invariant-subspace basis is singular; P cannot be extracted. **Covers two different situations -- see below** |
| `dare_error::non_psd_solution` | the extracted P is not positive semi-definite. The test is an LDLT pivot-sign test against `N * eps * max\|P_ij\|`, the order of the factorization's own backward error, below which a pivot carries no sign information |
| `dare_error::schur_failed` | the real Schur factorization did not converge |
| `dare_error::arithmetic_limit` | finite inputs could not produce a verified stabilizing solution at the scalar type's precision. Specifically: the equilibrated answer's estimated relative forward error exceeded `sqrt(eps)`, so it retains less than half the significand; its closed-loop spectrum or gain solve did not verify; a magnitude the equilibrated verification needs could not be formed; the weights overflowed while being equilibrated; the solution overflowed while being rescaled; the rescaled solution failed the positive-semi-definiteness test at its returned scale; the direct check at the caller's scale resolved and refuted the carried claim; or the equilibrated and returned gains, both formed, disagreed by more than the counted-operation margin |

The one cause that is **no longer** on that list is a magnitude the check at the
caller's scale could not form. That used to be a refusal; it is now an absence of
evidence, and the carried claim is what decides.

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
| `1e-1` … `1e-7` | 2 | value, `\|\|P\|\|` rising from `8.9e2` to `8.1e14` |
| `1e-8` … `1e-14` | 2 | **`singular_u11`** -- a well-posed pair refused |
| `1e-16`, `0` | 1 | `singular_u11` -- genuinely uncontrollable |

The transition sits where `\|\|P\|\|` approaches `1/eps` and the smallest pivot of the computed `U11` falls under the rank test's threshold. Naming this enumerator after the structural cause would state something false about every row in the middle band.

So `singular_u11` covers a pair with no stabilizing solution **and** a numerical failure to separate the invariant subspace, and **it does not distinguish them**. Telling them apart needs a stabilizability test the solver does not perform.

`care_error::singular_u11` has the identical shape for the identical reason, with the unit disk replaced by the open left half-plane and the reciprocal by the reflection `-lambda`.

### dare (with cross-weight)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
auto dare(const Matrix<Scalar, NX, NX>& A,
          const Matrix<Scalar, NX, NU>& B,
          const Matrix<Scalar, NX, NX>& Q,
          const Matrix<Scalar, NU, NU>& R,
          const Matrix<Scalar, NX, NU>& N)
    -> ctrlpp::expected<dare_result<Scalar, NX>, dare_error>;
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
