# discrete_filter

C++20 concepts defining the minimal interface for digital filters: `discrete_filter` for a scalar filter and `vector_discrete_filter` for a filter over a vector sample. Generic code constrained on them accepts any filter type with that interface.

## Header

| Form | Header |
|------|--------|
| `discrete_filter`, `vector_discrete_filter` (concepts) | `#include <ctrlpp/dsp/discrete_filter.h>` |
| (convenience) | `#include <ctrlpp/dsp.h>` |

## Concept Definition

```cpp
template <typename F>
concept discrete_filter = requires(F f, typename F::scalar_type x) {
    { f.process(x) } -> std::convertible_to<typename F::scalar_type>;
};

template <typename F, typename Vec>
concept vector_discrete_filter = requires(F f, Vec v) {
    { f.process(v) } -> std::convertible_to<Vec>;
};
```

## What Satisfies It

A type satisfies `discrete_filter` when it provides:

1. A nested type alias `scalar_type` naming the numeric type (e.g., `double`, `float`).
2. A member function `process(scalar_type) -> scalar_type` (or returning something convertible to `scalar_type`).

`vector_discrete_filter<F, Vec>` asks only that `process(Vec)` return something convertible to `Vec`; the vector type is a parameter of the concept, not a member of `F`.

Neither concept requires `reset()`; this is a convention followed by the library's built-in filters but not enforced at the concept level.

## Types That Satisfy discrete_filter

- `biquad<Scalar>`, a second-order IIR filter
- `fir<Scalar, N>`, a finite impulse response filter
- `cascaded_biquad<Scalar, N>`, a cascade of biquad sections

`vector_biquad<Scalar, N>` and `vector_cascaded_biquad<Scalar, N, Sections>` satisfy `vector_discrete_filter` with `Vec = Vector<Scalar, N>`; see [biquad](biquad.md#vector_biquad).

## Example Model Implementation

```cpp
#include <ctrlpp/dsp/discrete_filter.h>

#include <span>
#include <array>
#include <iostream>

// A simple exponential moving average satisfying the concept
struct ema
{
    using scalar_type = double;

    double alpha{0.1};
    double state{0.0};

    double process(double x)
    {
        state = alpha * x + (1.0 - alpha) * state;
        return state;
    }
};

static_assert(ctrlpp::discrete_filter<ema>);

// Generic function that works with any discrete_filter
template <ctrlpp::discrete_filter F>
void filter_signal(F& f, std::span<double> data)
{
    for(auto& v : data)
        v = f.process(v);
}

int main()
{
    ema filter{.alpha = 0.2};

    std::array<double, 8> data{1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0};
    filter_signal(filter, data);

    for(auto v : data)
        std::cout << v << " ";
    std::cout << "\n";
}
```

## See Also

- [biquad](biquad.md)<br/> second-order IIR filter
- [fir](fir.md)<br/> finite impulse response filter
- [background/dsp](../../background/dsp.md)<br/> DSP theory and background
