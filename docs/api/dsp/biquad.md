# biquad

Second-order IIR (biquad) filter implemented with transposed direct form II. Provides factory functions for common filter types and satisfies the `discrete_filter` concept for composability. Individual biquad sections can be cascaded via `cascaded_biquad` for higher-order filters, and convenience functions `make_butterworth` and `make_chebyshev1` build complete cascaded designs from cutoff frequency and sample rate.

Every design factory validates its parameters and returns `ctrlpp::expected<Filter, dsp_error>`: a filter on success, or the specific `dsp_error` enumerator describing the rejected design. See [Design Validation: dsp_error](#design-validation-dsp_error).

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

## Design Validation: dsp_error

Header: `#include <ctrlpp/dsp/dsp_types.h>` (pulled in by `biquad.h`)

```cpp
enum class dsp_error {
    non_positive_sample_rate,
    cutoff_exceeds_nyquist,
    non_positive_q,
    non_positive_ripple,
    non_finite_input,
};
```

Every design factory returns `ctrlpp::expected<Filter, dsp_error>` and rejects invalid parameters instead of computing garbage coefficients:

| Enumerator | Rejected input |
|------------|----------------|
| `non_finite_input` | Any design parameter (frequency, sample rate, quality factor, ripple) is NaN or infinite, or the design chain produced a non-finite coefficient |
| `non_positive_sample_rate` | `sample_hz <= 0` |
| `cutoff_exceeds_nyquist` | The design frequency lies outside the open interval `(0, sample_hz / 2)`. A discrete-time filter can only realize a response strictly below half the sample rate; at or above it the design frequency aliases (Nyquist criterion) |
| `non_positive_q` | Notch quality factor `q <= 0` |
| `non_positive_ripple` | Chebyshev Type I passband ripple `ripple_db <= 0`. The ripple factor is `eps = sqrt(10^(ripple_db / 10) - 1)`, whose radicand is non-positive for every `ripple_db <= 0`, and at exactly zero the following `asinh(1 / eps)` takes an infinite argument. An equiripple passband is defined by a strictly positive ripple, so this is an exact domain bound and carries no tolerance |

The checks run in the order listed, so a design with several defects reports the first matching enumerator. The Nyquist bound is strict on both sides: `cutoff == sample_hz / 2` is rejected. The zero-ripple boundary is rejected too, along with negative zero: a Chebyshev Type I design with no ripple is not a degenerate Butterworth, it is outside the design's domain. Use `make_butterworth` for a maximally flat passband.

`make_chebyshev1` additionally sweeps the coefficients it is about to emit and reports `non_finite_input` rather than returning them, so a successful design never hands back a filter whose difference equation immediately contaminates its state.

## Factory Functions

### low_pass

```cpp
static auto low_pass(Scalar cutoff_hz, Scalar sample_hz)
    -> ctrlpp::expected<biquad, dsp_error>;
```

Creates a second-order Butterworth (maximally flat, quality factor `Q = 1/sqrt(2)`) low-pass filter at the given cutoff frequency. The response is monotone across the passband with no peaking and reaches -3.01 dB at the cutoff. `cutoff_hz` must lie in the open interval `(0, sample_hz / 2)`.

### notch

```cpp
static auto notch(Scalar freq_hz, Scalar sample_hz, Scalar q)
    -> ctrlpp::expected<biquad, dsp_error>;
```

Creates a notch (band-reject) filter centered at `freq_hz` with quality factor `q`. `freq_hz` must lie in the open interval `(0, sample_hz / 2)` and `q` must be strictly positive.

### dirty_derivative

```cpp
static auto dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz)
    -> ctrlpp::expected<biquad, dsp_error>;
```

Creates a band-limited differentiator with analog prototype `H(s) = wc*s / (s + wc)`, discretized via the bilinear transform. It acts as a true differentiator (`|H(f)| -> 2*pi*f`) up to the bandwidth `wc = 2*pi*bandwidth_hz`, above which it rolls off. `bandwidth_hz` must lie in the open interval `(0, sample_hz / 2)`.

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

Resets internal state to zero, or initializes the filter state such that a constant input of `value` would produce the corresponding steady-state output.

### coefficients

```cpp
auto coefficients() const -> biquad_coeffs<Scalar> const&;
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

- `process(Scalar x) -> Scalar`, filters one sample through all sections
- `reset()`, resets all sections
- `reset(Scalar value)`, steady-state initialization propagated through the cascade
- `section(std::size_t i) -> biquad<Scalar>&`, access individual sections

## Convenience Design Functions

### make_butterworth

```cpp
template <std::size_t Order, typename Scalar>
    requires (Order % 2 == 0 && Order >= 2)
auto make_butterworth(Scalar cutoff_hz, Scalar sample_hz)
    -> ctrlpp::expected<cascaded_biquad<Scalar, Order / 2>, dsp_error>;
```

Designs an `Order`-th order Butterworth low-pass filter as a cascade of `Order/2` biquad sections. `cutoff_hz` must lie in the open interval `(0, sample_hz / 2)`.

### make_chebyshev1

```cpp
template <std::size_t Order, typename Scalar>
    requires (Order % 2 == 0 && Order >= 2)
auto make_chebyshev1(Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db)
    -> ctrlpp::expected<cascaded_biquad<Scalar, Order / 2>, dsp_error>;
```

Designs an `Order`-th order Chebyshev Type I low-pass filter with the specified passband ripple. `cutoff_hz` must lie in the open interval `(0, sample_hz / 2)`, and `ripple_db` must be finite and strictly positive. A design that succeeds carries only finite coefficients.

## Usage Example

```cpp
// gnuplot: plot "< ./biquad_filter" using 1:2 with lines title "noisy", "" using 1:3 with lines title "filtered"
#include <ctrlpp/dsp/biquad.h>

#include <cmath>
#include <iostream>
#include <numbers>

int main()
{
    constexpr double sample_hz = 100.0;
    constexpr double cutoff_hz = 10.0;

    // Create a second-order low-pass filter. The factory returns
    // ctrlpp::expected<biquad, dsp_error>; .value() unwraps in this
    // exceptions-on snippet (embedded-clean code checks has_value()
    // and dereferences instead).
    auto lp = ctrlpp::biquad<double>::low_pass(cutoff_hz, sample_hz).value();

    // Create a 4th-order Butterworth low-pass
    auto butter4 = ctrlpp::make_butterworth<4>(cutoff_hz, sample_hz).value();

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

- [fir](fir.md)<br/> finite impulse response filter
- [discrete-filter](discrete-filter.md)<br/> discrete filter concept
- [background/dsp](../../background/dsp.md)<br/> DSP theory and background
