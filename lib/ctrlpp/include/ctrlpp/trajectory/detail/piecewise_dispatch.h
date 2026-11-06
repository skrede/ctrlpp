#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PIECEWISE_DISPATCH_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_PIECEWISE_DISPATCH_H

// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
// Machines and Robots", 2008, Ch. 3, Sec. 3.6, p.117-119

#include <array>
#include <tuple>
#include <cstddef>
#include <utility>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Derived, typename Scalar, typename PointType, typename... Segments>
class piecewise_base
{
public:
    explicit piecewise_base(Segments... segs)
        : m_segments{std::move(segs)...}
        , m_breakpoints{}
        , m_total_duration{}
    {
        compute_breakpoints(std::index_sequence_for<Segments...>{});
    }

    PointType evaluate(Scalar t) const
    {
        auto const tc = std::clamp(t, Scalar{0}, m_total_duration);
        return dispatch(tc, std::index_sequence_for<Segments...>{});
    }

    Scalar duration() const { return m_total_duration; }

    static constexpr std::size_t segment_count = sizeof...(Segments);

protected:
    std::tuple<Segments...> m_segments;
    std::array<Scalar, sizeof...(Segments)> m_breakpoints;
    Scalar m_total_duration;

private:
    template <std::size_t... Is>
    void compute_breakpoints(std::index_sequence<Is...>)
    {
        Scalar cumulative{0};
        ((m_breakpoints[Is] = (cumulative += std::get<Is>(m_segments).duration())), ...);
        m_total_duration = cumulative;
    }

    template <std::size_t I>
    PointType evaluate_segment(Scalar t) const
    {
        Scalar local_t = (I == 0) ? t : t - m_breakpoints[I - 1];
        return std::get<I>(m_segments).evaluate(local_t);
    }

    PointType evaluate_last() const
    {
        constexpr auto last = sizeof...(Segments) - 1;
        return std::get<last>(m_segments).evaluate(std::get<last>(m_segments).duration());
    }

    template <std::size_t... Is>
    PointType dispatch(Scalar t, std::index_sequence<Is...>) const
    {
        PointType result{};
        bool found = false;
        auto try_seg = [&]<std::size_t I>(std::integral_constant<std::size_t, I>)
        {
            if(found)
                return;
            if constexpr(I < sizeof...(Segments) - 1)
            {
                if(t < m_breakpoints[I])
                {
                    result = evaluate_segment<I>(t);
                    found = true;
                }
            }
            else
            {
                result = evaluate_segment<I>(t);
                found = true;
            }
        };
        (try_seg(std::integral_constant<std::size_t, Is>{}), ...);
        if(!found)
            result = evaluate_last();
        return result;
    }
};

}

#endif
