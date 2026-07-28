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

Solves the standard DARE. Forms the 2n x 2n symplectic matrix, computes its real Schur decomposition, reorders eigenvalues inside the unit disk to the top-left block, and extracts P = U21 * U11^{-1}. On success `result->P` is the stabilizing solution; `result->subspace_separation` and `result->reorder_complete` are conditioning diagnostics.

Refusals:

| Enumerator | Condition |
| --- | --- |
| `dare_error::non_stabilisable` | fewer than n eigenvalues of the symplectic spectrum lie inside the unit region |
| `dare_error::non_finite_input` | A, B, Q, R **or the assembled symplectic Z** contains NaN/Inf |
| `dare_error::singular_a` | A is rank-deficient to a scale-relative reciprocal-pivot tolerance, so the `A^{-T}` the pencil build needs does not exist |
| `dare_error::singular_u11` | the top-left block of the reordered invariant-subspace basis is singular; P cannot be extracted |
| `dare_error::non_psd_solution` | the extracted P is not positive semi-definite within an epsilon-scaled tolerance |
| `dare_error::schur_failed` | the real Schur factorization did not converge |

Two of these are worth reading together, because the count test cannot see every unstabilizable pair. An uncontrollable mode at `|lambda| > 1` contributes **both** `lambda` and its reciprocal to the symplectic spectrum, so n eigenvalues do lie inside the unit disk and `non_stabilisable` does not fire; the subspace they span fails to project instead, and the refusal arrives as `singular_u11`. An uncontrollable mode at `|lambda| = 1` contributes two eigenvalues **on** the circle, neither inside, so that case does reach `non_stabilisable`. Both are refusals -- no gain is returned for an unstabilizable pair either way -- but only the second names the structural cause.

A singular `R` is `non_finite_input` rather than an enumerator of its own: the pencil build forms `R^{-1}` and the assembled Z goes infinite, which is the condition that enumerator's own definition covers.

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
        // P_opt.error() names which of the six conditions refused the problem.
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
