# place

Pole placement via Ackermann's formula for single-input systems. Computes a state-feedback gain K such that the eigenvalues of (A - BK) match the desired closed-loop poles. Also provides a dual `place_observer` function for computing observer gains via the duality (A', C') -> L'.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::place<Scalar, NX, NU>` | `#include <ctrlpp/control/place.h>` |
| `ctrlpp::place<Scalar, NX, NU>` | `#include <ctrlpp/place.h>` (convenience) |

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (e.g. `double`, `float`) |
| `NX` | `std::size_t` | State dimension |
| `NU` | `std::size_t` | Input dimension (must be 1 for `place`) |

## Functions

### place

```cpp
template <typename Scalar, std::size_t NX, std::size_t NU>
std::optional<Eigen::Matrix<Scalar, int(NU), int(NX)>>
place(const Matrix<Scalar, NX, NX>& A,
      const Matrix<Scalar, NX, NU>& B,
      const std::array<std::complex<Scalar>, NX>& desired_poles);
```

Computes K such that eig(A - BK) = desired_poles using Ackermann's formula. Complex poles must appear in conjugate pairs. Returns `std::nullopt` if NU > 1 (multi-input not supported), if the conjugate pair requirement is violated, or if the system is uncontrollable.

### place_observer

```cpp
template <typename Scalar, std::size_t NX, std::size_t NY>
std::optional<Eigen::Matrix<Scalar, int(NX), int(NY)>>
place_observer(const Matrix<Scalar, NX, NX>& A,
               const Matrix<Scalar, NY, NX>& C,
               const std::array<std::complex<Scalar>, NX>& desired_poles);
```

Computes observer gain L via the duality L = place(A', C', poles)'. Requires NY = 1 (single-output). The observer update becomes x_hat += L * (z - C * x_hat).

## Usage Example

```cpp
#include <ctrlpp/control/place.h>

#include <Eigen/Dense>

#include <array>
#include <complex>
#include <iostream>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;

    Eigen::Matrix2d A;
    A << 0.0, 1.0, -2.0, -3.0;

    Eigen::Matrix<double, 2, 1> B;
    B << 0.0, 1.0;

    // Place closed-loop poles at -1 +/- j
    std::array<std::complex<double>, NX> poles = {
        std::complex<double>{-1.0, 1.0},
        std::complex<double>{-1.0, -1.0}
    };

    auto K_opt = ctrlpp::place<double, NX, NU>(A, B, poles);
    if (!K_opt) {
        std::cerr << "Pole placement failed\n";
        return 1;
    }

    std::cout << "K = " << *K_opt << "\n";

    // Verify closed-loop eigenvalues
    Eigen::Matrix2d Acl = A - B * (*K_opt);
    Eigen::EigenSolver<Eigen::Matrix2d> es(Acl);
    std::cout << "Closed-loop eigenvalues:\n" << es.eigenvalues() << "\n";
}
```

## See Also

- [lqr](lqr.md) -- optimal gain design as an alternative to pole placement
- [dare](dare.md) -- Riccati equation solver used by LQR
- [reference/optimal-control-theory](../reference/optimal-control-theory.md) -- mathematical background
