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
    invalid_resolution_factor,
    cutoff_below_coefficient_resolution,
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
| `invalid_resolution_factor` | The `resolution_factor` of `low_pass`, `notch`, `dirty_derivative`, `make_butterworth` or `make_chebyshev1` is NaN, infinite, or not strictly positive, so the refusal bound it scales has no meaning |
| `cutoff_below_coefficient_resolution` | A section's design gap `D = 1 + a1 + a2` (`1 + a1` for the first-order dirty derivative) is at or below `resolution_factor * eps`, so the stored coefficients cannot carry one correct digit of it. See [Coefficient-resolution floor](#coefficient-resolution-floor) |

The low-pass family (`low_pass`, `make_butterworth`, `make_chebyshev1`) checks, in order: the ripple (Chebyshev only), non-finite parameters, the sample rate, the Nyquist interval, the resolution factor, a cutoff ratio that rounds onto zero or Nyquist in the pre-warp `tan(pi * cutoff_hz / sample_hz)` (`cutoff_exceeds_nyquist`), and then, per section, a non-finite coefficient followed by the coefficient-resolution floor; a design with several defects reports the first matching enumerator. The Nyquist bound is strict on both sides: `cutoff == sample_hz / 2` is rejected. The zero-ripple boundary is rejected too, along with negative zero: a Chebyshev Type I design with no ripple is not a degenerate Butterworth, it is outside the design's domain. Use `make_butterworth` for a maximally flat passband.

`notch` checks a non-finite `q` first, then the frequency, sample rate and Nyquist checks above, then `q <= 0`, then the same resolution factor, pre-warp, non-finite and floor checks. `dirty_derivative` checks in the low-pass order. A `q` so small that `1 / q` overflows, or a `sample_hz` so large that the derivative gain `2 * sample_hz` overflows, is refused with `non_finite_input`.

Every design factory sweeps the coefficients it is about to emit and reports `non_finite_input` rather than returning them, so a successful design never hands back a filter whose difference equation immediately contaminates its state.

## Coefficient-resolution floor

Every design factory refuses a design whose gap its stored coefficients cannot carry: the low-pass family below, the notch and the dirty derivative in their own subsections.

A low-pass section's poles approach `z = 1` as the cutoff ratio `fc / fs` falls, and its behavior lives in the gap `D = 1 + a1 + a2` (the DC gain is `(b0 + b1 + b2) / D`), which shrinks as `(2 pi fc / fs)^2`. Stored in a scalar type with machine epsilon `eps`, `a1` near -2 and `a2` near 1 are each rounded to half an ulp, `eps / 2` and `eps / 4`, so even the best coefficients the type can hold carry an error of up to `3/4 eps` in `D` (Oppenheim and Schafer, 3rd ed., Sec. 6.8, coefficient quantization with poles near `z = 1`). The relative error of the design's gap is then

```
rho = resolution_factor * eps / D
```

and every low-pass-family factory refuses a design with `rho >= 1`, where `D = 4 beta0 / (1 + beta1 + beta0)` is computed from the design (never from the stored coefficients) and, for a cascade, the smallest section gap decides. `resolution_factor` defaults to the counted `3/4`. A caller who needs the gap resolved to a relative accuracy `tau` passes `3/4 / tau`; for example `0.75f / 1e-3f` refuses a float low-pass at 1 Hz / 1000 Hz, which the default accepts.

The sections are formed at unit scale from `Omega = tan(pi * fc / fs)` with `a1` and `a2` stored as offsets from -2 and 1, so the gap is never formed by cancellation and the representation, not the formula, is what limits a design: across the float sweep below, every float coefficient is within its counted rounding bound of the double design rounded to float.

| Design | `float` floor `fc / fs` | `double` floor `fc / fs` |
|--------|-------------------------|--------------------------|
| `low_pass` | `4.76e-5` | `2.05e-9` |
| `make_butterworth<4>` | `4.76e-5` | `2.05e-9` |
| `make_chebyshev1<4>`, 1 dB ripple | `9.00e-5` (second section, pole radius squared 0.279) | `3.89e-9` |

`rho` and the measured error of a float `low_pass` at `fs = 1000`, evaluated in double from the stored coefficients:

| `fc / fs` | `rho` | DC gain error | `|H(fc)|` error |
|-----------|-------|---------------|-----------------|
| `1e-1` | `3.3e-7` | `2.2e-7` | `8.7e-8` |
| `1e-2` | `2.4e-5` | `3.0e-6` | `-2.3e-7` |
| `1e-3` | `2.3e-3` | `-9.0e-4` | `1.3e-7` |
| `5e-4` | `9.1e-3` | `-4.7e-3` | `-7.2e-6` |
| `2e-4` | `5.7e-2` | `-2.0e-2` | `-9.3e-5` |
| `1e-4` | `0.23` | `-5.4e-2` | `-8.1e-4` |
| `7e-5` | `0.46` | `-0.19` | `-1.3e-2` |
| `5e-5` | `0.91` | `-0.17` | `-1.1e-2` |
| `4.8e-5` | `0.98` | `-0.24` | `-2.3e-2` |
| `3e-5` | `2.5` | refused | refused |

Just above the floor the design is accepted but not accurate: the storage error alone can move the DC gain by up to `rho / (1 - rho)`, which grows without bound as `rho` approaches 1, and the measured DC error is 17 percent at `5e-5` and 24 percent at `4.8e-5`. The floor refuses a design whose gap has no correct digit; it does not promise accuracy near it. Pass a larger `resolution_factor` for that. `double` reaches the floor only at ratios about four decades below `float`'s.

### Notch floor

The notch shares the low-pass section's denominator (`g1 = 1 / q`, `g0 = 1`), and its zeros sit at `exp(+-j w0)`, also near `z = 1`. Its numerator is stored as offsets too: `b1 = a1` and `b0 = b2 = 1 - beta1 / (1 + beta1 + beta0)`, so the numerator gap `b0 + b1 + b2` equals the design gap `D` and the DC gain is 1. The DC gain is the ratio of the two gaps, and storage moves the numerator gap by up to `eps` (`b1` `eps / 2`, `b0` and `b2` `eps / 4` each) on top of the denominator's `3/4 eps`, so `resolution_factor` defaults to the counted `7/4`. The floor depends on `q` only through `beta1 = Omega / q`, which is small next to 1 at the floor: the float floor is `7.27e-5` for `q = 0.7` and for `q = 10` alike.

`rho` and the measured error of a float `notch` with `q = 0.7` at `fs = 1000`; the response at `fs / 4`, far from the notch, is compared with the double design:

| `f / fs` | `rho` | DC gain error | `|H(fs / 4)|` error |
|----------|-------|---------------|---------------------|
| `1e-1` | `7.8e-7` | `0` | `-2.0e-8` |
| `1e-2` | `5.5e-5` | `-1.6e-5` | `5.7e-10` |
| `1e-3` | `5.3e-3` | `0` | `1.6e-11` |
| `3e-4` | `5.9e-2` | `1.7e-2` | `-7.5e-12` |
| `2e-4` | `0.13` | `-3.7e-2` | `9.9e-12` |
| `1e-4` | `0.53` | `0` | `5.9e-12` |
| `8e-5` | `0.83` | `-0.20` | `3.3e-13` |
| `7.5e-5` | `0.94` | `0.33` | `-4.1e-12` |
| `7e-5` | `1.08` | refused | refused |

The DC gain error is quantized: near the floor the stored numerator and denominator gaps differ by `2 b0 - 1 - a2`, a multiple of `eps / 2`, so the DC gain is either exactly 1 or off by at least `eps / 2` over the stored gap, 33 percent at `7.5e-5`. The response away from the notch is unaffected.

### Dirty-derivative floor

The dirty derivative is first order: its pole sits at `-a1` near `z = 1`, and its behavior lives in the gap `D1 = 1 + a1 = 2 Omega / (1 + Omega)`, which equals `w0 = 2 pi bandwidth_hz / sample_hz` to first order. `a1` is stored as the offset `-1 + D1`, so the gap is never formed by cancellation, and half an ulp of `a1` near -1 is `eps / 4`, its only storage error in the gap. `resolution_factor` defaults to the counted `1/4`, and the factory refuses `rho = resolution_factor * eps / D1 >= 1`.

`rho` and the measured error of a float `dirty_derivative` at `fs = 1000`: the stored gap `1 + a1` against the design's `D1`, and the response at the bandwidth against the double design:

| `bandwidth / fs` | `rho` | gap error | `|H(bandwidth)|` error |
|------------------|-------|-----------|------------------------|
| `1e-3` | `4.8e-6` | `4.1e-6` | `-2.1e-6` |
| `1e-5` | `4.7e-4` | `-1.1e-4` | `5.2e-5` |
| `1e-6` | `4.7e-3` | `-3.9e-3` | `2.0e-3` |
| `1e-7` | `4.7e-2` | `4.4e-2` | `-2.2e-2` |
| `1e-8` | `0.47` | `-5.1e-2` | `2.6e-2` |
| `7e-9` | `0.68` | `0.36` | `-0.16` |
| `5e-9` | `0.95` | `0.90` | `-0.34` |
| `4.8e-9` | `0.99` | `0.98` | `-0.36` |
| `4.5e-9` | `1.05` | refused | refused |

Just above the floor the stored gap can be off by up to `rho`, nearly 100 percent, and the response at the bandwidth by a third. As for the low-pass, pass a larger `resolution_factor` for accuracy near the floor. The floors with their default factors:

| Design | Default `resolution_factor` | `float` floor | `double` floor |
|--------|-----------------------------|---------------|----------------|
| `notch`, `q = 0.7` | `7/4` | `7.27e-5` | `3.14e-9` |
| `dirty_derivative` | `1/4` | `4.74e-9` | `8.83e-18` |

### Vector forwarders

The vector forwarders take the same trailing parameter and forward the scalar factory's verdict unchanged:

```cpp
static auto vector_biquad<Scalar, N>::low_pass(Scalar cutoff_hz, Scalar sample_hz,
    Scalar resolution_factor = /* 3/4 */) -> ctrlpp::expected<vector_biquad, dsp_error>;

static auto vector_biquad<Scalar, N>::notch(Scalar freq_hz, Scalar sample_hz, Scalar q,
    Scalar resolution_factor = /* 7/4 */) -> ctrlpp::expected<vector_biquad, dsp_error>;

static auto vector_biquad<Scalar, N>::dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz,
    Scalar resolution_factor = /* 1/4 */) -> ctrlpp::expected<vector_biquad, dsp_error>;

template <std::size_t Order, std::size_t N, typename Scalar>
auto make_vector_butterworth(Scalar cutoff_hz, Scalar sample_hz,
    Scalar resolution_factor = /* 3/4 */)
    -> ctrlpp::expected<vector_cascaded_biquad<Scalar, N, Order / 2>, dsp_error>;

template <std::size_t Order, std::size_t N, typename Scalar>
auto make_vector_chebyshev1(Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db,
    Scalar resolution_factor = /* 3/4 */)
    -> ctrlpp::expected<vector_cascaded_biquad<Scalar, N, Order / 2>, dsp_error>;
```

## Factory Functions

### low_pass

```cpp
static auto low_pass(Scalar cutoff_hz, Scalar sample_hz,
    Scalar resolution_factor = /* 3/4 */)
    -> ctrlpp::expected<biquad, dsp_error>;
```

Creates a second-order Butterworth (maximally flat, quality factor `Q = 1/sqrt(2)`) low-pass filter at the given cutoff frequency. The response is monotone across the passband with no peaking and reaches -3.01 dB at the cutoff. `cutoff_hz` must lie in the open interval `(0, sample_hz / 2)` and above the [coefficient-resolution floor](#coefficient-resolution-floor) of `Scalar`.

### notch

```cpp
static auto notch(Scalar freq_hz, Scalar sample_hz, Scalar q,
    Scalar resolution_factor = /* 7/4 */)
    -> ctrlpp::expected<biquad, dsp_error>;
```

Creates a notch (band-reject) filter centered at `freq_hz` with quality factor `q`. `freq_hz` must lie in the open interval `(0, sample_hz / 2)` and above the [notch floor](#notch-floor) of `Scalar`, and `q` must be strictly positive.

### dirty_derivative

```cpp
static auto dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz,
    Scalar resolution_factor = /* 1/4 */)
    -> ctrlpp::expected<biquad, dsp_error>;
```

Creates a band-limited differentiator with analog prototype `H(s) = wc*s / (s + wc)`, discretized via the bilinear transform. It acts as a true differentiator (`|H(f)| -> 2*pi*f`) up to the bandwidth `wc = 2*pi*bandwidth_hz`, above which it rolls off. `bandwidth_hz` must lie in the open interval `(0, sample_hz / 2)` and above the [first-order floor](#dirty-derivative-floor) of `Scalar`.

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
auto make_butterworth(Scalar cutoff_hz, Scalar sample_hz,
    Scalar resolution_factor = /* 3/4 */)
    -> ctrlpp::expected<cascaded_biquad<Scalar, Order / 2>, dsp_error>;
```

Designs an `Order`-th order Butterworth low-pass filter as a cascade of `Order/2` biquad sections. `cutoff_hz` must lie in the open interval `(0, sample_hz / 2)` and every section above the [coefficient-resolution floor](#coefficient-resolution-floor) of `Scalar`.

### make_chebyshev1

```cpp
template <std::size_t Order, typename Scalar>
    requires (Order % 2 == 0 && Order >= 2)
auto make_chebyshev1(Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db,
    Scalar resolution_factor = /* 3/4 */)
    -> ctrlpp::expected<cascaded_biquad<Scalar, Order / 2>, dsp_error>;
```

Designs an `Order`-th order Chebyshev Type I low-pass filter with the specified passband ripple. `cutoff_hz` must lie in the open interval `(0, sample_hz / 2)` with every section above the [coefficient-resolution floor](#coefficient-resolution-floor) of `Scalar`, and `ripple_db` must be finite and strictly positive. The section with the smaller pole radius sets the floor, so it sits above the Butterworth floor of the same order. A design that succeeds carries only finite coefficients.

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
