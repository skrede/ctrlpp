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
reporting success, it verifies the Riccati residual against a dimension- and
precision-derived forward-error margin and verifies that the resulting
closed-loop spectrum is inside the unit disk by more than its backward-error
margin. The checks run on both the equilibrated problem and the returned
solution in the caller's original scale.

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
| `dare_error::arithmetic_limit` | finite inputs could not produce a residual-verified stabilizing solution at the scalar type's precision, including overflow while equilibrating or unscaling, and including a returned P that fails the positive-semi-definiteness test at its returned scale after unscaling |

Positive semi-definiteness is tested twice, and deliberately so. The extraction
primitive tests the P the solve produced; when common-weight equilibration is
active that P is the *scaled* one, and the value handed back to the caller is
`P * weight_scale`, which the first test never saw. The postcondition therefore
re-tests at the returned scale. Scaling by a positive factor preserves
definiteness mathematically, so this second test only ever fires on the rounding
the rescale itself introduces -- which is why it reports `arithmetic_limit`
rather than `non_psd_solution`, and why the floor has to carry the
factorization's backward error rather than assume there is none.

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
