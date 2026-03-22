#ifndef HPP_GUARD_CTRLPP_DSP_FIR_H
#define HPP_GUARD_CTRLPP_DSP_FIR_H

#include "ctrlpp/dsp/discrete_filter.h"

#include <array>
#include <cstddef>

namespace ctrlpp {

template<typename Scalar, std::size_t N>
class fir {
public:
    using scalar_type = Scalar;

    explicit constexpr fir(std::array<Scalar, N> taps) : taps_{taps} {}

    auto process(Scalar x) -> Scalar {
        for (std::size_t i = N - 1; i > 0; --i) {
            delay_[i] = delay_[i - 1];
        }
        delay_[0] = x;

        Scalar y{};
        for (std::size_t i = 0; i < N; ++i) {
            y += taps_[i] * delay_[i];
        }
        return y;
    }

    void reset() { delay_.fill(Scalar{0}); }

    [[nodiscard]] auto taps() const -> std::array<Scalar, N> const& {
        return taps_;
    }

private:
    std::array<Scalar, N> taps_{};
    std::array<Scalar, N> delay_{};
};

template<typename Scalar, std::size_t N>
fir(std::array<Scalar, N>) -> fir<Scalar, N>;

namespace detail {

static_assert(discrete_filter<fir<double, 3>>);
static_assert(discrete_filter<fir<float, 4>>);

}

}

#endif
