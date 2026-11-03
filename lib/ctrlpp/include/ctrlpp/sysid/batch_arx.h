#ifndef HPP_GUARD_CTRLPP_SYSID_BATCH_ARX_H
#define HPP_GUARD_CTRLPP_SYSID_BATCH_ARX_H

// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 4

#include "ctrlpp/sysid/fit_metrics.h"
#include "ctrlpp/sysid/sysid_types.h"
#include "ctrlpp/sysid/sysid_result.h"

#include "ctrlpp/sysid/detail/batch_arx_detail.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <cstddef>
#include <utility>
#include <optional>
#include <algorithm>

namespace ctrlpp
{

// The checks are exact preconditions of the index arithmetic, not identifiability conditions: a fit with fewer
// regressor rows than NA + NB, or a rank-deficient one, is returned and its diagnostics report the deficiency.
template <std::size_t NA, std::size_t NB, typename Derived1, typename Derived2>
expected<arx_result<typename Derived1::Scalar, std::max(NA, NB), 1, 1>, sysid_error> batch_arx(
        const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U, std::optional<typename Derived1::Scalar> rank_threshold = std::nullopt)
{
    using Scalar = typename Derived1::Scalar;
    static constexpr std::size_t NX = std::max(NA, NB);

    if(Y.cols() != U.cols())
        return unexpected(sysid_error::record_length_mismatch);
    if(Y.rows() != 1 || U.rows() != 1)
        return unexpected(sysid_error::record_not_single_row);
    if(Y.cols() <= static_cast<Eigen::Index>(NX))
        return unexpected(sysid_error::too_few_samples);
    if(!Y.allFinite() || !U.allFinite())
        return unexpected(sysid_error::non_finite_sample);

    auto const regression = detail::assemble_arx_regression<NA, NB>(Y, U);
    auto const fit = detail::solve_arx_least_squares(regression, rank_threshold);

    // Built in place inside the expected, so the result is not moved on return.
    expected<arx_result<Scalar, NX, 1, 1>, sysid_error> result{std::in_place};
    auto& identified = *result;
    identified.system = detail::realize_arx_observer_form<NA, NB>(fit.theta);
    auto const y_predicted = detail::simulate_arx_open_loop(identified.system, U);
    Eigen::VectorX<Scalar> const y_actual = Y.row(0).transpose();
    identified.metrics = detail::compute_fit_metrics_unchecked(y_actual, y_predicted);
    identified.diagnostics = fit.diagnostics;
    return result;
}

}

#endif
