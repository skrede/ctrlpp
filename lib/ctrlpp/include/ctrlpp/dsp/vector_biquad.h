#ifndef HPP_GUARD_CTRLPP_DSP_VECTOR_BIQUAD_H
#define HPP_GUARD_CTRLPP_DSP_VECTOR_BIQUAD_H

// @cite oppenheim2010dsp -- Oppenheim & Schafer, "Discrete-Time Signal Processing", 3rd ed., 2010, Sec. 6.4, Eq. (6.44) (transposed direct form II)

#include "ctrlpp/dsp/biquad.h"
#include "ctrlpp/dsp/dsp_types.h"
#include "ctrlpp/dsp/discrete_filter.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cstddef>
#include <utility>

namespace ctrlpp
{

namespace detail
{

template <typename T, std::size_t N, std::size_t... Is>
constexpr std::array<T, N> filled_array_impl(T const& value, std::index_sequence<Is...>)
{
    return {{(static_cast<void>(Is), value)...}};
}

template <typename T, std::size_t N>
constexpr std::array<T, N> filled_array(T const& value)
{
    return filled_array_impl<T, N>(value, std::make_index_sequence<N>{});
}

}

template <ctrlpp_floating_scalar Scalar, std::size_t N>
    requires(N >= 1)
class vector_biquad
{
public:
    using scalar_type = Scalar;

    constexpr vector_biquad()
        : m_channels{}
    {
    }

    explicit constexpr vector_biquad(std::array<biquad<Scalar>, N> channels)
        : m_channels{channels}
    {
    }

    Vector<Scalar, N> process(Vector<Scalar, N> x)
    {
        Vector<Scalar, N> y;
        for(std::size_t i = 0; i < N; ++i)
        {
            y[static_cast<Eigen::Index>(i)] = m_channels[i].process(x[static_cast<Eigen::Index>(i)]);
        }
        return y;
    }

    void reset()
    {
        for(auto& ch : m_channels)
            ch.reset();
    }

    void reset(Vector<Scalar, N> value)
    {
        for(std::size_t i = 0; i < N; ++i)
        {
            m_channels[i].reset(value[static_cast<Eigen::Index>(i)]);
        }
    }

    static expected<vector_biquad, dsp_error> low_pass(Scalar cutoff_hz, Scalar sample_hz, Scalar resolution_factor = detail::second_order_resolution_count<Scalar>())
    {
        auto const proto = biquad<Scalar>::low_pass(cutoff_hz, sample_hz, resolution_factor);
        if(!proto.has_value())
            return unexpected(proto.error());
        return vector_biquad{detail::filled_array<biquad<Scalar>, N>(*proto)};
    }

    static expected<vector_biquad, dsp_error> notch(Scalar freq_hz, Scalar sample_hz, Scalar q, Scalar resolution_factor = detail::notch_resolution_count<Scalar>())
    {
        auto const proto = biquad<Scalar>::notch(freq_hz, sample_hz, q, resolution_factor);
        if(!proto.has_value())
            return unexpected(proto.error());
        return vector_biquad{detail::filled_array<biquad<Scalar>, N>(*proto)};
    }

    static expected<vector_biquad, dsp_error> dirty_derivative(Scalar bandwidth_hz, Scalar sample_hz, Scalar resolution_factor = detail::first_order_resolution_count<Scalar>())
    {
        auto const proto = biquad<Scalar>::dirty_derivative(bandwidth_hz, sample_hz, resolution_factor);
        if(!proto.has_value())
            return unexpected(proto.error());
        return vector_biquad{detail::filled_array<biquad<Scalar>, N>(*proto)};
    }

private:
    std::array<biquad<Scalar>, N> m_channels;
};

template <ctrlpp_floating_scalar Scalar, std::size_t N, std::size_t Sections>
    requires(N >= 1 && Sections >= 1)
class vector_cascaded_biquad
{
public:
    using scalar_type = Scalar;

    explicit constexpr vector_cascaded_biquad(std::array<cascaded_biquad<Scalar, Sections>, N> channels)
        : m_channels{channels}
    {
    }

    Vector<Scalar, N> process(Vector<Scalar, N> x)
    {
        Vector<Scalar, N> y;
        for(std::size_t i = 0; i < N; ++i)
        {
            y[static_cast<Eigen::Index>(i)] = m_channels[i].process(x[static_cast<Eigen::Index>(i)]);
        }
        return y;
    }

    void reset()
    {
        for(auto& ch : m_channels)
            ch.reset();
    }

    void reset(Vector<Scalar, N> value)
    {
        for(std::size_t i = 0; i < N; ++i)
        {
            m_channels[i].reset(value[static_cast<Eigen::Index>(i)]);
        }
    }

private:
    std::array<cascaded_biquad<Scalar, Sections>, N> m_channels;
};

template <std::size_t Order, std::size_t N, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
expected<vector_cascaded_biquad<Scalar, N, Order / 2>, dsp_error> make_vector_butterworth(
        Scalar cutoff_hz, Scalar sample_hz, Scalar resolution_factor = detail::second_order_resolution_count<Scalar>())
{
    auto const proto = make_butterworth<Order>(cutoff_hz, sample_hz, resolution_factor);
    if(!proto.has_value())
        return unexpected(proto.error());
    return vector_cascaded_biquad<Scalar, N, Order / 2>{detail::filled_array<cascaded_biquad<Scalar, Order / 2>, N>(*proto)};
}

template <std::size_t Order, std::size_t N, typename Scalar>
    requires(Order % 2 == 0 && Order >= 2)
expected<vector_cascaded_biquad<Scalar, N, Order / 2>, dsp_error> make_vector_chebyshev1(
        Scalar cutoff_hz, Scalar sample_hz, Scalar ripple_db, Scalar resolution_factor = detail::second_order_resolution_count<Scalar>())
{
    auto const proto = make_chebyshev1<Order>(cutoff_hz, sample_hz, ripple_db, resolution_factor);
    if(!proto.has_value())
        return unexpected(proto.error());
    return vector_cascaded_biquad<Scalar, N, Order / 2>{detail::filled_array<cascaded_biquad<Scalar, Order / 2>, N>(*proto)};
}

namespace detail
{

static_assert(vector_discrete_filter<vector_biquad<double, 1>, Vector<double, 1>>);
static_assert(vector_discrete_filter<vector_biquad<double, 3>, Vector<double, 3>>);
static_assert(vector_discrete_filter<vector_cascaded_biquad<double, 2, 2>, Vector<double, 2>>);

}

}

#endif
