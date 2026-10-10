#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_BSPLINE_CHECKS_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_BSPLINE_CHECKS_H

#include "ctrlpp/trajectory/trajectory_types.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <vector>
#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
auto all_finite(std::vector<Scalar> const& values) -> bool
{
    return std::all_of(values.begin(), values.end(), [](Scalar value) { return std::isfinite(value); });
}

template <typename Scalar>
auto check_interpolation_input(std::vector<Scalar> const& times, std::vector<Scalar> const& positions)
    -> ctrlpp::expected<void, spline_error>
{
    if (!all_finite(times) || !all_finite(positions)) {
        return ctrlpp::unexpected(spline_error::non_finite_input);
    }
    for (std::size_t i = 0; i + 1 < times.size(); ++i) {
        if (!(times[i + 1] > times[i])) {
            return ctrlpp::unexpected(spline_error::non_increasing_times);
        }
        if (!std::isfinite(times[i + 1] - times[i])) {
            return ctrlpp::unexpected(spline_error::unrepresentable_spline);
        }
    }
    return {};
}

template <int Degree, typename Scalar>
auto check_knots(std::size_t control_point_count, std::vector<Scalar> const& knot_vector)
    -> ctrlpp::expected<void, spline_error>
{
    if (knot_vector.empty()) {
        return {};
    }
    if (knot_vector.size() != control_point_count + static_cast<std::size_t>(Degree) + 1) {
        return ctrlpp::unexpected(spline_error::bad_knot_count);
    }
    if (!std::is_sorted(knot_vector.begin(), knot_vector.end())) {
        return ctrlpp::unexpected(spline_error::non_monotonic_knots);
    }
    if (!std::isfinite(knot_vector.back() - knot_vector.front())) {
        return ctrlpp::unexpected(spline_error::unrepresentable_spline);
    }
    if (!(knot_vector[static_cast<std::size_t>(Degree)] < knot_vector[control_point_count])) {
        return ctrlpp::unexpected(spline_error::non_increasing_times);
    }
    return {};
}

}

#endif
