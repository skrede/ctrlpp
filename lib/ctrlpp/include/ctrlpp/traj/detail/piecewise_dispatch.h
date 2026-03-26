#ifndef HPP_GUARD_CTRLPP_TRAJ_DETAIL_PIECEWISE_DISPATCH_H
#define HPP_GUARD_CTRLPP_TRAJ_DETAIL_PIECEWISE_DISPATCH_H

/// @brief Shared dispatch logic for piecewise path and trajectory compositions.
///
/// Provides a CRTP base template containing breakpoint computation, time dispatch,
/// and segment evaluation. Both piecewise_path and piecewise_trajectory inherit
/// from this base, differing only in their point type and concept constraints.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Ch. 4 (composite trajectories)

#include <algorithm>
#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

namespace ctrlpp::detail
{

/// @brief CRTP base for piecewise compositions with shared dispatch logic.
///
/// Derived must provide a type alias `point_type` for the return type of evaluate().
template <typename Derived, typename Scalar, typename PointType, typename... Segments>
class piecewise_base
{
public:
    explicit piecewise_base(Segments... segs)
        : segments_{std::move(segs)...}
    {
        compute_breakpoints(std::index_sequence_for<Segments...>{});
    }

    auto evaluate(Scalar t) const -> PointType
    {
        auto const tc = std::clamp(t, Scalar{0}, total_dur_);
        return dispatch(tc, std::index_sequence_for<Segments...>{});
    }

    auto duration() const -> Scalar { return total_dur_; }

    static constexpr std::size_t segment_count = sizeof...(Segments);

protected:
    std::tuple<Segments...> segments_;
    std::array<Scalar, sizeof...(Segments)> breakpoints_{};
    Scalar total_dur_{};

private:
    template <std::size_t... Is>
    void compute_breakpoints(std::index_sequence<Is...>)
    {
        Scalar cumulative{0};
        ((breakpoints_[Is] = (cumulative += std::get<Is>(segments_).duration())), ...);
        total_dur_ = cumulative;
    }

    template <std::size_t I>
    auto evaluate_segment(Scalar t) const -> PointType
    {
        Scalar local_t = (I == 0) ? t : t - breakpoints_[I - 1];
        return std::get<I>(segments_).evaluate(local_t);
    }

    auto evaluate_last() const -> PointType
    {
        constexpr auto last = sizeof...(Segments) - 1;
        return std::get<last>(segments_).evaluate(std::get<last>(segments_).duration());
    }

    template <std::size_t... Is>
    auto dispatch(Scalar t, std::index_sequence<Is...>) const -> PointType
    {
        PointType result{};
        bool found = false;
        auto try_seg = [&]<std::size_t I>(std::integral_constant<std::size_t, I>)
        {
            if (found) return;
            if constexpr (I < sizeof...(Segments) - 1)
            {
                if (t < breakpoints_[I])
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
        if (!found)
        {
            result = evaluate_last();
        }
        return result;
    }
};

} // namespace ctrlpp::detail

#endif
