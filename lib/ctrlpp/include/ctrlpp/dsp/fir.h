#ifndef HPP_GUARD_CTRLPP_DSP_FIR_H
#define HPP_GUARD_CTRLPP_DSP_FIR_H

// @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Sec. 6.5.1, Eq. (6.46) (direct-form FIR convolution)

#include "ctrlpp/dsp/discrete_filter.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cstddef>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t N>
class fir
{
    static_assert(N >= 1, "fir requires at least one tap");

public:
    using scalar_type = Scalar;

    explicit constexpr fir(std::array<Scalar, N> taps)
        : m_taps{taps}
        , m_delay{}
    {
    }

    Scalar process(Scalar x)
    {
        for(std::size_t i = N - 1; i > 0; --i)
        {
            m_delay[i] = m_delay[i - 1];
        }
        m_delay[0] = x;

        Scalar y{};
        for(std::size_t i = 0; i < N; ++i)
        {
            y += m_taps[i] * m_delay[i];
        }
        return y;
    }

    void reset() { m_delay.fill(Scalar{0}); }

    std::array<Scalar, N> const& taps() const { return m_taps; }

private:
    std::array<Scalar, N> m_taps;
    std::array<Scalar, N> m_delay;
};

template <typename Scalar, std::size_t N>
fir(std::array<Scalar, N>) -> fir<Scalar, N>;

namespace detail
{

static_assert(discrete_filter<fir<double, 3>>);
static_assert(discrete_filter<fir<float, 4>>);

}

}

#endif
