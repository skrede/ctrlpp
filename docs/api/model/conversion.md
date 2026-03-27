# conversion

Conversion between transfer function and state-space representations. `tf2ss` produces controllable canonical form; `ss2tf` uses the Leverrier-Faddeev algorithm. Both operate on SISO systems.

## Header and Alias

| Form | Header |
|------|--------|
| `tf2ss(tf)` | `#include <ctrlpp/model/conversion.h>` |
| `ss2tf(sys)` | `#include <ctrlpp/model/conversion.h>` |
| (convenience) | `#include <ctrlpp/conversion.h>` |

## Functions

### tf2ss

```cpp
template <typename Scalar, std::size_t NumDeg, std::size_t DenDeg>
    requires (NumDeg <= DenDeg)
constexpr continuous_state_space<Scalar, DenDeg, 1, 1>
tf2ss(const transfer_function<Scalar, NumDeg, DenDeg>& tf);
```

Converts a transfer function to controllable canonical form state-space. Requires a proper transfer function (`NumDeg <= DenDeg`). The state dimension equals `DenDeg`. Coefficients are automatically normalised to a monic denominator.

### ss2tf

```cpp
template <typename Scalar, std::size_t NX>
transfer_function<Scalar, NX, NX>
ss2tf(const continuous_state_space<Scalar, NX, 1, 1>& sys);
```

Converts a SISO continuous state-space to a transfer function using the Leverrier-Faddeev algorithm. Returns coefficients in highest-degree-first order.

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

    std::cout << "tf2ss result:\n"
              << "  A =\n" << sys.A << "\n"
              << "  B = " << sys.B.transpose() << "\n"
              << "  C = " << sys.C << "\n"
              << "  D = " << sys.D << "\n\n";

    // Round-trip back
    auto tf2 = ctrlpp::ss2tf(sys);
    std::cout << "ss2tf round-trip:\n  num = [";
    for(auto c : tf2.numerator)
        std::cout << " " << c;
    std::cout << " ]\n  den = [";
    for(auto c : tf2.denominator)
        std::cout << " " << c;
    std::cout << " ]\n";
}
```

## See Also

- [state-space](state-space.md) -- state-space representation
- [transfer-function](transfer-function.md) -- transfer function representation
- [discretise](discretise.md) -- continuous-to-discrete conversion
