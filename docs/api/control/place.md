# place

Pole placement via Ackermann's formula for single-input systems. Computes a state-feedback gain K such that the eigenvalues of (A - BK) match the desired closed-loop poles. Also provides a dual `place_observer` function for computing observer gains via the duality (A', C') -> L'.

## Header and Alias

| Form | Header |
|------|--------|
| `ctrlpp::place<Scalar, NX, NU>` | `#include <ctrlpp/control/place.h>` |
| `ctrlpp::place<Scalar, NX, NU>` | `#include <ctrlpp/place.h>` (convenience) |
| `ctrlpp::place_error` | `#include <ctrlpp/control/place_types.h>` |

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
auto place(const Matrix<Scalar, NX, NX>& A,
           const Matrix<Scalar, NX, NU>& B,
           const std::array<std::complex<Scalar>, NX>& desired_poles)
    -> ctrlpp::expected<Eigen::Matrix<Scalar, int(NU), int(NX)>, place_error>;
```

Computes K such that eig(A - BK) = desired_poles using Ackermann's formula. Complex poles must appear in conjugate pairs.

Refusals, checked in order:

| Condition | Enumerator |
| --- | --- |
| `NU > 1` | `place_error::multi_input_not_supported` |
| the pole set is not closed under conjugation | `place_error::poles_not_conjugate_symmetric` |
| the controllability matrix is rank-deficient | `place_error::uncontrollable_pair` |

### place_observer

```cpp
template <typename Scalar, std::size_t NX, std::size_t NY>
auto place_observer(const Matrix<Scalar, NX, NX>& A,
                    const Matrix<Scalar, NY, NX>& C,
                    const std::array<std::complex<Scalar>, NX>& desired_poles)
    -> ctrlpp::expected<Eigen::Matrix<Scalar, int(NX), int(NY)>, place_error>;
```

Computes observer gain L via the duality L = place(A', C', poles)'. Requires NY = 1 (single-output). The observer update becomes x_hat += L * (z - C * x_hat).

`NY > 1` is `place_error::multi_output_not_supported`; everything else is the dual placement's own refusal forwarded, so `place_error::uncontrollable_pair` on the transposed pair reports that `(A, C)` is unobservable.

### place_error

```cpp
enum class place_error
{
    multi_input_not_supported,
    multi_output_not_supported,
    poles_not_conjugate_symmetric,
    uncontrollable_pair,
};
```

Declared in `<ctrlpp/control/place_types.h>`, which `place.h` includes.

The first two are **structural**: no choice of data or poles makes them succeed, because the single-channel Ackermann formula does not cover the shape at all. Multi-input assignment is not even a unique problem -- it has a free subspace, which is what a robust-assignment method exists to choose within. The last two are **conditions on the data**, which a different pole set or a different pair can satisfy. They are separate enumerators because telling a caller to change the poles when the input dimension is the obstacle sends them after something that cannot help.

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
    if (!K_opt.has_value()) {
        // K_opt.error() names which of the four conditions refused the design.
        std::cerr << "Pole placement refused the design\n";
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

- [lqr](lqr.md)<br/> optimal gain design as an alternative to pole placement
- [dare](dare.md)<br/> Riccati equation solver used by LQR
- [luenberger](../estimation/luenberger.md)<br/> observer using place_observer for gain design
