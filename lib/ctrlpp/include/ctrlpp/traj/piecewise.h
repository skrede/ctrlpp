#ifndef HPP_GUARD_CTRLPP_TRAJ_PIECEWISE_H
#define HPP_GUARD_CTRLPP_TRAJ_PIECEWISE_H

/// @brief Variadic piecewise trajectory composition.
///
/// Composes heterogeneous trajectory segments into a single trajectory that
/// evaluates seamlessly across segment boundaries with automatic time offset.
/// The piecewise compositor itself satisfies trajectory_segment, enabling
/// recursive composition.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic
/// Machines and Robots", 2009, Ch. 4 (composite trajectories)

#include "ctrlpp/traj/cubic_segment.h"
#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/trajectory_types.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <tuple>
#include <utility>

namespace ctrlpp
{

/// @brief Composes heterogeneous trajectory segments into a single trajectory.
///
/// Segments are evaluated in sequence with automatic time offset. At breakpoints,
/// the next segment is evaluated at local t=0 to ensure continuity.
template <typename Scalar, std::size_t ND, typename... Segments>
    requires(trajectory_segment<Segments, Scalar, ND> && ...)
class piecewise
{
public:
    explicit piecewise(Segments... segs)
        : segments_{std::move(segs)...}
    {
        compute_breakpoints(std::index_sequence_for<Segments...>{});
    }

    auto evaluate(Scalar t) const -> trajectory_point<Scalar, ND>
    {
        auto const tc = std::clamp(t, Scalar{0}, total_dur_);
        return dispatch(tc, std::index_sequence_for<Segments...>{});
    }

    auto duration() const -> Scalar { return total_dur_; }

    static constexpr std::size_t segment_count = sizeof...(Segments);

private:
    std::tuple<Segments...> segments_;
    std::array<Scalar, sizeof...(Segments)> breakpoints_{};
    Scalar total_dur_{};

    template <std::size_t... Is>
    void compute_breakpoints(std::index_sequence<Is...>)
    {
        Scalar cumulative{0};
        ((breakpoints_[Is] = (cumulative += std::get<Is>(segments_).duration())), ...);
        total_dur_ = cumulative;
    }

    template <std::size_t I>
    auto evaluate_segment(Scalar t) const -> trajectory_point<Scalar, ND>
    {
        Scalar local_t = (I == 0) ? t : t - breakpoints_[I - 1];
        return std::get<I>(segments_).evaluate(local_t);
    }

    auto evaluate_last() const -> trajectory_point<Scalar, ND>
    {
        constexpr auto last = sizeof...(Segments) - 1;
        return std::get<last>(segments_).evaluate(std::get<last>(segments_).duration());
    }

    template <std::size_t... Is>
    auto dispatch(Scalar t, std::index_sequence<Is...>) const -> trajectory_point<Scalar, ND>
    {
        trajectory_point<Scalar, ND> result{};
        bool found = false;
        auto try_seg = [&]<std::size_t I>(std::integral_constant<std::size_t, I>)
        {
            if (found) return;
            if constexpr (I < sizeof...(Segments) - 1)
            {
                // Non-last segments: at breakpoint, hand off to next segment
                if (t < breakpoints_[I])
                {
                    result = evaluate_segment<I>(t);
                    found = true;
                }
            }
            else
            {
                // Last segment: evaluate up to and including its endpoint
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

/// Verify piecewise itself satisfies trajectory_segment (recursive composition).
static_assert(trajectory_segment<
    piecewise<double, 1, cubic_segment<double, 1>, cubic_segment<double, 1>>, double, 1>);

} // namespace ctrlpp

#endif
