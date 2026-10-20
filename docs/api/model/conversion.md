# conversion

Conversion between transfer function and state-space representations. `tf2ss` produces controllable canonical form; `ss2tf` forms both polynomials from the eigenvalues of a diagonally balanced matrix. Both operate on SISO systems and both are fallible: an input they cannot convert is refused with a `conversion_error`, never returned as a result with infinities or a wrong polynomial.

## Header and Alias

| Form | Header |
|------|--------|
| `tf2ss(tf)` | `#include <ctrlpp/model/conversion.h>` |
| `ss2tf(sys)` | `#include <ctrlpp/model/conversion.h>` |
| `ctrlpp::conversion_error` | `#include <ctrlpp/model/conversion_types.h>` |
| (convenience) | `#include <ctrlpp/conversion.h>` |

## Functions

### tf2ss

```cpp
template <typename Scalar, std::size_t NumDeg, std::size_t DenDeg>
    requires (NumDeg <= DenDeg)
auto tf2ss(const transfer_function<Scalar, NumDeg, DenDeg>& tf)
    -> ctrlpp::expected<continuous_state_space<Scalar, DenDeg, 1, 1>, conversion_error>;
```

Converts a transfer function to controllable canonical form state-space. Requires a proper transfer function (`NumDeg <= DenDeg`). The state dimension equals `DenDeg`. Coefficients are normalized to a monic denominator.

Refusals, checked in order:

| Condition | Enumerator |
| --- | --- |
| a numerator or denominator coefficient is NaN or infinite | `conversion_error::non_finite_input` |
| `tf.denominator[0] == 0`: the denominator's true degree is below `DenDeg`, so no companion form of that dimension exists | `conversion_error::zero_leading_coefficient` |
| the coefficients are finite but a normalized one overflows (a tiny leading coefficient) | `conversion_error::unrepresentable_result` |

### ss2tf

```cpp
template <typename Scalar, std::size_t NX>
auto ss2tf(const continuous_state_space<Scalar, NX, 1, 1>& sys)
    -> ctrlpp::expected<transfer_function<Scalar, NX, NX>, conversion_error>;
```

Converts a SISO continuous state-space to H(s) = C (sI - A)^-1 B + D, coefficients highest-degree-first.

The denominator is the characteristic polynomial of A. A copy of A is balanced by a power-of-two diagonal similarity (LAPACK DGEBAL phase 2, exact, so the eigenvalues are A's own), its eigenvalues are computed by the QR algorithm, and the monic product of the real linear factors and the real quadratic factors of conjugate pairs gives the coefficients. The numerator goes through the same routine by the rank-one identity

det(sI - A + alpha B C) = det(sI - A) (1 + alpha C (sI - A)^-1 B),

so the numerator is D poly(A) + (poly(A - alpha B C) - poly(A)) / alpha for any nonzero alpha. `ss2tf` takes alpha as the power of two nearest ||A||_F / (||B|| ||C||), so that alpha B C is on the scale of A and the division by alpha is exact. The leading numerator coefficient is D exactly, and B = 0 or C = 0 returns D poly(A).

Refusals, checked in order:

| Condition | Enumerator |
| --- | --- |
| an entry of A, B, C or D is NaN or infinite | `conversion_error::non_finite_input` |
| the eigensolver does not converge | `conversion_error::eigenvalues_not_converged` |
| the input is finite but a coefficient is not | `conversion_error::unrepresentable_result` |

**Accuracy.** The coefficients carry the eigenvalues' rounding, so a round trip through `tf2ss` is not bit-exact even on small integer cases; it is bounded. The QR algorithm's eigenvalues are exact for the balanced matrix plus a perturbation of norm about n u ||A_bal||_F (u the unit roundoff), so each eigenvalue moves by at most r_i = kappa_i n u ||A_bal||_F to first order, with kappa_i its condition number. Writing e_k for the k-th coefficient of the polynomial with roots |lambda_i|, coefficient k moves by at most e_k(|lambda| + r) - e_k(|lambda|), plus the factor product's own rounding, gamma_m e_k(|lambda| + r) with gamma_m = m u / (1 - m u) and m counting two roundings per real factor and five per conjugate pair. The numerator inherits the two polynomials' bounds divided by alpha. On eight real poles from -0.1 to -1000 the recovered coefficients agree with the originals to about 4e-13 relative in the numerator and 2e-13 in the denominator.

### conversion_error

```cpp
enum class conversion_error
{
    non_finite_input,
    zero_leading_coefficient,
    eigenvalues_not_converged,
    unrepresentable_result,
};
```

Declared in `<ctrlpp/model/conversion_types.h>`, which `conversion.h` includes.

## Usage Example

```cpp
#include <ctrlpp/model/conversion.h>
#include <ctrlpp/model/transfer_function.h>
#include <ctrlpp/model/state_space.h>

#include <iostream>

int main()
{
    // H(s) = (2s + 3) / (s^2 + 5s + 6)
    ctrlpp::transfer_function<double, 1, 2> tf{
        .numerator = {2.0, 3.0},
        .denominator = {1.0, 5.0, 6.0}};

    auto sys = ctrlpp::tf2ss(tf);
    if(!sys)
    {
        std::cerr << "tf2ss refused the transfer function\n";
        return 1;
    }

    std::cout << "tf2ss result:\n"
              << "  A =\n" << sys->A << "\n"
              << "  B = " << sys->B.transpose() << "\n"
              << "  C = " << sys->C << "\n"
              << "  D = " << sys->D << "\n\n";

    // Round-trip back
    auto tf2 = ctrlpp::ss2tf(*sys);
    if(!tf2)
    {
        std::cerr << "ss2tf refused the system\n";
        return 1;
    }
    std::cout << "ss2tf round-trip:\n  num = [";
    for(auto c : tf2->numerator)
        std::cout << " " << c;
    std::cout << " ]\n  den = [";
    for(auto c : tf2->denominator)
        std::cout << " " << c;
    std::cout << " ]\n";
}
```

## See Also

- [state-space](state-space.md)<br/> state-space representation
- [transfer-function](transfer-function.md)<br/> transfer function representation
- [discretize](discretize.md)<br/> continuous-to-discrete conversion
