# biquad

Second-order IIR (biquad) filter implemented with transposed direct form II. Provides factory functions for common filter types and satisfies the `discrete_filter` concept for composability. Individual biquad sections can be cascaded via `cascaded_biquad` for higher-order filters, and convenience functions `make_butterworth` and `make_chebyshev1` build complete cascaded designs from cutoff frequency and sample rate.

## Header and Alias

| Form | Header |
|------|--------|
| `biquad<Scalar>` | `#include <ctrlpp/dsp/biquad.h>` |
| (convenience) | `#include <ctrlpp/dsp.h>` |

```cpp
template <typename Scalar>
class biquad;
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | floating-point | Numeric type (`double`, `float`) |

## Supporting Types

### biquad_coeffs

```cpp
template <typename Scalar>
struct biquad_coeffs {
    Scalar b0{}, b1{}, b2{};
    Scalar a1{}, a2{};
};
```

Normalized biquad coefficients. The `a0` coefficient is implicitly 1 (already divided out in the factory functions).

## Constructors

```cpp
constexpr biquad() = default;
explicit constexpr biquad(biquad_coeffs<Scalar> c);
```

Construct from explicit coefficients, or use one of the factory functions below.

## Factory Functions

### low_pass

```cpp
static auto low_pass(Scalar cutoff_hz, Scalar sample_hz) -> biquad;
```

Creates a second-order Butterworth low-pass filter at the given cutoff frequency.

### notch

```cpp
static auto notch(Scalar freq_hz, Scalar sample_hz, Scalar q) -> biquad;
```

Creates a notch (band-reject) filter centred at `freq_hz` with quality factor `q`.

### dirty_derivative

```cpp
static auto dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz) -> biquad;
```

Creates a filtered derivative (high-pass with roll-off) via bilinear transform.

## Methods

### process

```cpp
auto process(Scalar x) -> Scalar;
```

Filters a single sample through the biquad section and returns the output.

### reset

```cpp
void reset();
void reset(Scalar value);
```

Resets internal state to zero, or initialises the filter state such that a constant input of `value` would produce the corresponding steady-state output.

### coefficients

```cpp
[[nodiscard]] auto coefficients() const -> biquad_coeffs<Scalar> const&;
```

Returns the current filter coefficients.

## cascaded_biquad

```cpp
template <typename Scalar, std::size_t N>
    requires (N >= 1)
class cascaded_biquad;
```

Chains N biquad sections in series. Each `process()` call passes the sample through all sections sequentially.

### Methods

- `process(Scalar x) -> Scalar` -- filters one sample through all sections
- `reset()` -- resets all sections
- `reset(Scalar value)` -- steady-state initialisation propagated through the cascade
- `section(std::size_t i) -> biquad<Scalar>&` -- access individual sections

## Convenience Design Functions

### make_butterworth

```cpp
template <std::size_t Order, typename Scalar>
    requires (Order % 2 == 0 && Order >= 2)
auto make_butterworth(Scalar cutoff_hz, Scalar sample_hz)
    -> cascaded_biquad<Scalar, Order / 2>;
```

Designs an `Order`-th order Butterworth low-pass filter as a cascade of `Order/2` biquad sections.

### make_chebyshev1

```cpp
template <std::size_t Order, typename Scalar>
    requires (Order % 2 == 0 && Order >= 2)
auto make_chebyshev1(Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db)
    -> cascaded_biquad<Scalar, Order / 2>;
```

Designs an `Order`-th order Chebyshev Type I low-pass filter with the specified passband ripple.

## Usage Example

```cpp
#include "ctrlpp/dsp/biquad.h"

#include <cmath>
#include <iostream>
#include <numbers>

int main()
{
    constexpr double sample_hz = 100.0;
    constexpr double cutoff_hz = 10.0;

    // Create a second-order low-pass filter
    auto lp = ctrlpp::biquad<double>::low_pass(cutoff_hz, sample_hz);

    // Create a 4th-order Butterworth low-pass
    auto butter4 = ctrlpp::make_butterworth<4>(cutoff_hz, sample_hz);

    // Filter a noisy sine wave
    constexpr double signal_hz = 5.0;
    constexpr double noise_hz = 45.0;
    constexpr double dt = 1.0 / sample_hz;

    for(int k = 0; k < 200; ++k)
    {
        double t = static_cast<double>(k) * dt;
        double signal = std::sin(2.0 * std::numbers::pi * signal_hz * t);
        double noise = 0.3 * std::sin(2.0 * std::numbers::pi * noise_hz * t);
        double noisy = signal + noise;

        double filtered_2nd = lp.process(noisy);
        double filtered_4th = butter4.process(noisy);

        std::cout << t << "," << noisy << ","
                  << filtered_2nd << "," << filtered_4th << "\n";
    }
}
```

## See Also

- [fir](fir.md) -- finite impulse response filter
- [discrete-filter](discrete-filter.md) -- discrete filter concept
- [reference/dsp-theory](../reference/dsp-theory.md) -- DSP theory and background
