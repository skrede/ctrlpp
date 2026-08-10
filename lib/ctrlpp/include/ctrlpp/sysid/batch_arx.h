#ifndef HPP_GUARD_CTRLPP_SYSID_BATCH_ARX_H
#define HPP_GUARD_CTRLPP_SYSID_BATCH_ARX_H

/// @brief Batch ARX model identification via QR decomposition.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 4

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/sysid/fit_metrics.h"
#include "ctrlpp/sysid/sysid_types.h"
#include "ctrlpp/sysid/sysid_result.h"

#include "ctrlpp/sysid/detail/batch_arx_detail.h"

#include <Eigen/Dense>

#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

/// @brief Identify an ARX(NA, NB) model from a single-input single-output record.
///
/// Rejections, checked in order, all of them exact preconditions of the routine's
/// own index arithmetic and data layout rather than tolerances:
///  * Y and U hold different sample counts -> sysid_error::record_length_mismatch.
///    Sample k of one is paired with sample k of the other, so unequal lengths
///    have no consistent pairing.
///  * either record has a row count other than one -> sysid_error::record_not_single_row.
///    The body reads only row zero, so a multi-row record would be identified
///    from a fraction of its data.
///  * the sample count is at or below max(NA, NB) -> sysid_error::too_few_samples.
///    One regressor row is formed per sample beyond max(NA, NB), so a strictly
///    greater sample count is exactly the condition for a nonempty regressor
///    matrix. At the order the matrix is empty; below it the row count is the
///    result of an unsigned subtraction that wraps.
///  * a sample in either record is NaN or infinite -> sysid_error::non_finite_sample.
///    Samples enter the regressor and the least-squares solve directly.
///
/// These are validity conditions on the arithmetic, not identifiability
/// conditions on the fit. In particular the routine does NOT require the
/// regressor row count to reach the parameter count NA + NB, so a record that
/// passes every check can still yield an under-determined or rank-deficient
/// fit; the returned `fit_metrics` are the place to judge that.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 4
template <std::size_t NA, std::size_t NB, typename Derived1, typename Derived2>
auto batch_arx(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U)
    -> expected<arx_result<typename Derived1::Scalar, std::max(NA, NB), 1, 1>, sysid_error>
{
    using Scalar = typename Derived1::Scalar;
    static constexpr std::size_t NP = NA + NB;
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
    Eigen::Matrix<Scalar, static_cast<int>(NP), 1> theta = regression.Phi.colPivHouseholderQr().solve(regression.target);

    auto const sys = detail::realize_arx_observer_form<NA, NB>(theta);
    auto const y_predicted = detail::simulate_arx_open_loop(sys, U);
    Eigen::VectorX<Scalar> const y_actual = Y.row(0).transpose();

    auto metrics = detail::compute_fit_metrics_unchecked(y_actual, y_predicted);

    return arx_result<Scalar, NX, 1, 1>{.system = sys, .metrics = metrics};
}

}

#endif
