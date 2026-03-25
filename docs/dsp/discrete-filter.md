# discrete_filter

C++23 concept defining the minimal interface for composable digital filters. Any type satisfying this concept can be used interchangeably in filter chains, generic signal processing pipelines, and policy-based compositions.

## Header

| Form | Header |
|------|--------|
| `discrete_filter` (concept) | `#include <ctrlpp/dsp/discrete_filter.h>` |
| (convenience) | `#include <ctrlpp/dsp.h>` |

## Concept Definition

```cpp
template <typename F>
concept discrete_filter = requires(F f, typename F::scalar_type x) {
    { f.process(x) } -> std::convertible_to<typename F::scalar_type>;
};
```

## What Satisfies It

A type satisfies `discrete_filter` when it provides:

1. A nested type alias `scalar_type` naming the numeric type (e.g., `double`, `float`).
2. A member function `process(scalar_type) -> scalar_type` (or returning something convertible to `scalar_type`).

The concept does not require `reset()` -- that is a convention followed by the library's built-in filters but not enforced at the concept level.

## Types That Satisfy discrete_filter

- `biquad<Scalar>` -- second-order IIR filter
- `fir<Scalar, N>` -- finite impulse response filter
- `cascaded_biquad<Scalar, N>` -- cascaded biquad sections

## Example Model Implementation

```cpp
#include <ctrlpp/dsp/discrete_filter.h>

#include <iostream>

// A simple exponential moving average satisfying the concept
struct ema
{
    using scalar_type = double;

    double alpha{0.1};
    double state{0.0};

    auto process(double x) -> double
    {
        state = alpha * x + (1.0 - alpha) * state;
        return state;
    }
};

static_assert(ctrlpp::discrete_filter<ema>);

// Generic function that works with any discrete_filter
template <ctrlpp::discrete_filter F>
void filter_signal(F& f, double* data, std::size_t n)
{
    for(std::size_t i = 0; i < n; ++i)
        data[i] = f.process(data[i]);
}

int main()
{
    ema filter{.alpha = 0.2};

    double data[] = {1.0, 2.0, 3.0, 2.0, 1.0, 0.0, 1.0, 2.0};
    filter_signal(filter, data, 8);

    for(auto v : data)
        std::cout << v << " ";
    std::cout << "\n";
}
```

## See Also

- [biquad](biquad.md) -- second-order IIR filter
- [fir](fir.md) -- finite impulse response filter
