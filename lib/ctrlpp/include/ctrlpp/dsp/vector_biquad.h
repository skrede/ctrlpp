#ifndef HPP_GUARD_CTRLPP_DSP_VECTOR_BIQUAD_H
#define HPP_GUARD_CTRLPP_DSP_VECTOR_BIQUAD_H

/// @brief Vector-valued biquad filters wrapping N scalar biquads behind Eigen vector interface.
///
/// @cite oppenheim1997 -- Oppenheim & Willsky, "Signals and Systems", 1997

#include "ctrlpp/types.h"

#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/dsp/discrete_filter.h"

#include <array>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

namespace detail
{

template <typename T, std::size_t N, std::size_t... Is>
constexpr auto filled_array_impl(T const& value, std::index_sequence<Is...>) -> std::array<T, N>
{
    return {{(static_cast<void>(Is), value)...}};
}

template <typename T, std::size_t N>
constexpr auto filled_array(T const& value) -> std::array<T, N>
{
    return filled_array_impl<T, N>(value, std::make_index_sequence<N>{});
}

}

template <typename Scalar, std::size_t N>
    requires(N >= 1)
class vector_biquad
{
public:
    using scalar_type = Scalar;

    constexpr vector_biquad() = default;
    explicit constexpr vector_biquad(std::array<biquad<Scalar>, N> channels) : channels_{channels} {}

    auto process(Vector<Scalar, N> x) -> Vector<Scalar, N>
    {
        Vector<Scalar, N> y;
        for(std::size_t i = 0; i < N; ++i)
        {
            y[static_cast<int>(i)] = channels_[i].process(x[static_cast<int>(i)]);
        }
        return y;
    }

    void reset()
    {
        for(auto& ch : channels_)
            ch.reset();
    }

    void reset(Vector<Scalar, N> value)
    {
        for(std::size_t i = 0; i < N; ++i)
        {
            channels_[i].reset(value[static_cast<int>(i)]);
        }
    }

    static auto low_pass(Scalar cutoff_hz, Scalar sample_hz) -> vector_biquad
    {
        auto const proto = biquad<Scalar>::low_pass(cutoff_hz, sample_hz);
        return vector_biquad{detail::filled_array<biquad<Scalar>, N>(proto)};
    }

    static auto notch(Scalar freq_hz, Scalar sample_hz, Scalar q) -> vector_biquad
    {
        auto const proto = biquad<Scalar>::notch(freq_hz, sample_hz, q);
        return vector_biquad{detail::filled_array<biquad<Scalar>, N>(proto)};
    }

    static auto dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz) -> vector_biquad
    {
        auto const proto = biquad<Scalar>::dirty_derivative(bandwidth_hz, sample_hz);
        return vector_biquad{detail::filled_array<biquad<Scalar>, N>(proto)};
    }

private:
    std::array<biquad<Scalar>, N> channels_{};
};

template <typename Scalar, std::size_t N, std::size_t Sections>
    requires(N >= 1 && Sections >= 1)
class vector_cascaded_biquad
{
public:
    using scalar_type = Scalar;

    explicit constexpr vector_cascaded_biquad(std::array<cascaded_biquad<Scalar, Sections>, N> channels)
        : channels_{channels}
    {
    }

    auto process(Vector<Scalar, N> x) -> Vector<Scalar, N>
    {
        Vector<Scalar, N> y;
        for(std::size_t i = 0; i < N; ++i)
        {
            y[static_cast<int>(i)] = channels_[i].process(x[static_cast<int>(i)]);
        }
        return y;
    }

    void reset()
    {
        for(auto& ch : channels_)
            ch.reset();
    }

    void reset(Vector<Scalar, N> value)
    {
        for(std::size_t i = 0; i < N; ++i)
        {
            channels_[i].reset(value[static_cast<int>(i)]);
        }
    }

private:
    std::array<cascaded_biquad<Scalar, Sections>, N> channels_;
};

template <std::size_t Order, std::size_t N, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
auto make_vector_butterworth(Scalar cutoff_hz, Scalar sample_hz) -> vector_cascaded_biquad<Scalar, N, Order / 2>
{
    auto const proto = make_butterworth<Order>(cutoff_hz, sample_hz);
    return vector_cascaded_biquad<Scalar, N, Order / 2>{
        detail::filled_array<cascaded_biquad<Scalar, Order / 2>, N>(proto)};
}

template <std::size_t Order, std::size_t N, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
auto make_vector_chebyshev1(Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db)
    -> vector_cascaded_biquad<Scalar, N, Order / 2>
{
    auto const proto = make_chebyshev1<Order>(cutoff_hz, sample_hz, ripple_db);
    return vector_cascaded_biquad<Scalar, N, Order / 2>{
        detail::filled_array<cascaded_biquad<Scalar, Order / 2>, N>(proto)};
}

namespace detail
{

static_assert(vector_discrete_filter<vector_biquad<double, 1>, Vector<double, 1>>);
static_assert(vector_discrete_filter<vector_biquad<double, 3>, Vector<double, 3>>);
static_assert(vector_discrete_filter<vector_cascaded_biquad<double, 2, 2>, Vector<double, 2>>);

}

}

#endif
