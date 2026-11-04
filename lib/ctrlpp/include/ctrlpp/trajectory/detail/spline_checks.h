#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_SPLINE_CHECKS_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_SPLINE_CHECKS_H

#include "ctrlpp/trajectory/trajectory_types.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

// Three multiply-adds in the Horner form, two roundings each, and the six chained operations behind each of the quadratic
// and cubic coefficients, each worth up to one unit in the last place at the scale of the waypoint positions.
inline constexpr std::int32_t spline_position_rounding_ops = 3 * 2 + 6 + 6;

// The resolution an evaluated position is known at, and so the reach a term needs to move one.
template <typename Scalar>
inline Scalar spline_position_resolution(Scalar scale)
{
    return static_cast<Scalar>(spline_position_rounding_ops) * std::numeric_limits<Scalar>::epsilon() * scale;
}

template <typename Scalar>
inline Scalar largest_magnitude(std::vector<Scalar> const& values)
{
    Scalar largest{0};
    for(auto const value : values)
    {
        largest = std::max(largest, std::abs(value));
    }
    return largest;
}

template <typename Scalar>
inline ctrlpp::expected<void, spline_error> check_knot_spans(std::vector<Scalar> const& times)
{
    for(std::size_t i = 0; i + 1 < times.size(); ++i)
    {
        auto const span = times[i + 1] - times[i];
        if(!std::isfinite(span))
        {
            return ctrlpp::unexpected(spline_error::unrepresentable_spline);
        }
        if(!(span > Scalar{0}))
        {
            return ctrlpp::unexpected(spline_error::non_increasing_times);
        }
    }
    return {};
}

// v0 and vn are the clamped end velocities, refused with the knots and waypoints when not finite; a spline without them passes zeros.
template <typename Scalar>
inline ctrlpp::expected<void, spline_error> check_spline_knots(std::vector<Scalar> const& times, std::vector<Scalar> const& positions, Scalar v0, Scalar vn)
{
    if(times.size() < 2)
    {
        return ctrlpp::unexpected(spline_error::too_few_points);
    }
    if(positions.size() != times.size())
    {
        return ctrlpp::unexpected(spline_error::size_mismatch);
    }
    auto const finite = [](Scalar value) { return std::isfinite(value); };
    if(!std::all_of(times.begin(), times.end(), finite) || !std::all_of(positions.begin(), positions.end(), finite) || !std::isfinite(v0) || !std::isfinite(vn))
    {
        return ctrlpp::unexpected(spline_error::non_finite_input);
    }
    return check_knot_spans(times);
}

}

#endif
