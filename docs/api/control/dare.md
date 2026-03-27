# dare

Discrete Algebraic Riccati Equation solver using symplectic Schur decomposition. Finds the stabilising solution P to A'PA - P - A'PB(R + B'PB)^{-1}B'PA + Q = 0. This is the workhorse behind `lqr_gain` and the terminal cost computation in MPC.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::dare<Scalar, NX, NU>` | `#include <ctrlpp/control/dare.h>` |
| `ctrlpp::dare<Scalar, NX, NU>` | `#include <ctrlpp/dare.h>` (convenience) |

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
    -> std::optional<Matrix<Scalar, NX, NX>>;
```

Solves the standard DARE. Forms the 2n x 2n symplectic matrix, computes its complex Schur decomposition, reorders stable eigenvalues to the top-left block, and extracts P = U21 * U11^{-1}. Returns `std::nullopt` if A is singular, the system is not stabilisable, or the solution is not positive semi-definite.

### dare (with cross-weight)

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Matrix<Scalar, NX, NX>>
dare(const Matrix<Scalar, NX, NX>& A,
     const Matrix<Scalar, NX, NU>& B,
     const Matrix<Scalar, NX, NX>& Q,
     const Matrix<Scalar, NU, NU>& R,
     const Matrix<Scalar, NX, NU>& N);
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
    if (!P_opt) {
        std::cerr << "DARE failed to converge\n";
        return 1;
    }

    std::cout << "P =\n" << *P_opt << "\n";
}
```

## See Also

- [lqr](lqr.md) -- uses DARE internally to compute optimal gains
- [place](place.md) -- pole placement as an alternative design method
- [mpc](../mpc/mpc.md) -- MPC uses DARE for terminal cost computation
