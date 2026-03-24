# transfer_function

SISO transfer function representation as numerator and denominator polynomial coefficient arrays. Coefficients are stored in highest-degree-first order (MATLAB convention).

## Header and Alias

| Form | Header |
|------|--------|
| `transfer_function<Scalar, NumDeg, DenDeg>` | `#include <ctrlpp/model/transfer_function.h>` |

```cpp
template <typename Scalar, std::size_t NumDeg, std::size_t DenDeg>
struct transfer_function {
    std::array<Scalar, NumDeg + 1> numerator;
    std::array<Scalar, DenDeg + 1> denominator;
};
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type |
| `NumDeg` | `>= 0` | Numerator polynomial degree |
| `DenDeg` | `>= 1` | Denominator polynomial degree |

## Fields

| Field | Type | Description |
|-------|------|-------------|
| `numerator` | `std::array<Scalar, NumDeg + 1>` | Numerator coefficients, highest degree first |
| `denominator` | `std::array<Scalar, DenDeg + 1>` | Denominator coefficients, highest degree first |

## Usage Example

```cpp
#include "ctrlpp/model/transfer_function.h"
#include "ctrlpp/model/conversion.h"

#include <iostream>

int main()
{
    // Second-order system: H(s) = 1 / (s^2 + 2*s + 1)
    ctrlpp::transfer_function<double, 0, 2> tf{
        .numerator = {1.0},
        .denominator = {1.0, 2.0, 1.0}};

    // Convert to state-space
    auto sys = ctrlpp::tf2ss(tf);

    std::cout << "State-space from TF:\n"
              << "  A =\n" << sys.A << "\n"
              << "  B = " << sys.B.transpose() << "\n"
              << "  C = " << sys.C << "\n"
              << "  D = " << sys.D << "\n";

    // Round-trip: convert back to TF
    auto tf_back = ctrlpp::ss2tf(sys);

    std::cout << "\nRound-trip TF:\n  num = [";
    for(auto c : tf_back.numerator)
        std::cout << " " << c;
    std::cout << " ]\n  den = [";
    for(auto c : tf_back.denominator)
        std::cout << " " << c;
    std::cout << " ]\n";
}
```

## See Also

- [state-space](state-space.md) -- state-space representation
- [conversion](conversion.md) -- TF to SS and SS to TF conversion
