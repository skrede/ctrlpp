# fir

Finite Impulse Response filter with compile-time tap count. Implements a direct-form FIR convolution using a circular delay line and satisfies the `discrete_filter` concept for composability.

## Header and Alias

| Form | Header |
|------|--------|
| `fir<Scalar, N>` | `#include <ctrlpp/dsp/fir.h>` |
| (convenience) | `#include <ctrlpp/dsp.h>` |

```cpp
template <typename Scalar, std::size_t N>
class fir;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type (`double`, `float`) |
| `N` | `>= 1` | Number of filter taps (coefficients) |

## Constructors

```cpp
explicit constexpr fir(std::array<Scalar, N> taps);
```

Constructs the filter from an array of N tap coefficients. A CTAD deduction guide is provided:

```cpp
fir(std::array<Scalar, N>) -> fir<Scalar, N>;
```

## Methods

### process

```cpp
auto process(Scalar x) -> Scalar;
```

Shifts the new sample into the delay line and computes the convolution sum. Returns the filtered output.

### reset

```cpp
void reset();
```

Clears the delay line to zero.

### taps

```cpp
[[nodiscard]] auto taps() const -> std::array<Scalar, N> const&;
```

Returns the tap coefficient array.

## Usage Example

```cpp
#include "ctrlpp/dsp/fir.h"

#include <array>
#include <cmath>
#include <iostream>
#include <numbers>

int main()
{
    // 5-tap moving average filter
    constexpr std::size_t N = 5;
    std::array<double, N> taps{};
    for(auto& t : taps)
        t = 1.0 / static_cast<double>(N);

    ctrlpp::fir filter(taps);

    // Filter a noisy step signal
    constexpr double sample_hz = 100.0;
    constexpr double dt = 1.0 / sample_hz;

    for(int k = 0; k < 100; ++k)
    {
        double t = static_cast<double>(k) * dt;
        double step = (k >= 20) ? 1.0 : 0.0;
        double noise = 0.1 * std::sin(2.0 * std::numbers::pi * 40.0 * t);
        double noisy = step + noise;

        double filtered = filter.process(noisy);

        std::cout << t << "," << noisy << "," << filtered << "\n";
    }
}
```

## See Also

- [biquad](biquad.md) -- second-order IIR filter
- [discrete-filter](discrete-filter.md) -- discrete filter concept
- [reference/dsp-theory](../reference/dsp-theory.md) -- DSP theory and background
