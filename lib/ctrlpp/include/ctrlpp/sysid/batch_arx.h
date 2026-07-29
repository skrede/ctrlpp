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

    // Observer canonical realization dimension = max(deg A, deg B) (Ljung 1999, Ch. 4).
    // When NB > NA the extra b-coefficients b_{NA+1..NB} require additional states;
    // truncating to NA states would silently drop them from the transfer function.
    static constexpr std::size_t NX = std::max(NA, NB);

    if(Y.cols() != U.cols())
        return unexpected(sysid_error::record_length_mismatch);
    if(Y.rows() != 1 || U.rows() != 1)
        return unexpected(sysid_error::record_not_single_row);
    if(Y.cols() <= static_cast<Eigen::Index>(NX))
        return unexpected(sysid_error::too_few_samples);
    if(!Y.allFinite() || !U.allFinite())
        return unexpected(sysid_error::non_finite_sample);

    auto N = static_cast<std::size_t>(Y.cols());
    std::size_t k = NX;
    auto n_eff = static_cast<Eigen::Index>(N - k);

    // Build regressor matrix Phi (n_eff x NP) and target vector
    Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NP)> Phi(n_eff, static_cast<int>(NP));
    Eigen::VectorX<Scalar> Y_target(n_eff);

    for(Eigen::Index i = 0; i < n_eff; ++i)
    {
        auto row = static_cast<Eigen::Index>(k) + i;

        // y-regressors: [Y(0, row-1), Y(0, row-2), ..., Y(0, row-NA)]
        for(std::size_t j = 0; j < NA; ++j)
        {
            Phi(i, static_cast<Eigen::Index>(j)) = Y(0, row - static_cast<Eigen::Index>(j + 1));
        }
        // u-regressors: [U(0, row-1), U(0, row-2), ..., U(0, row-NB)]
        for(std::size_t j = 0; j < NB; ++j)
            Phi(i, static_cast<Eigen::Index>(NA + j)) = U(0, row - static_cast<Eigen::Index>(j + 1));

        Y_target(i) = Y(0, row);
    }

    // Solve via QR decomposition
    Eigen::Matrix<Scalar, static_cast<int>(NP), 1> theta = Phi.colPivHouseholderQr().solve(Y_target);

    // Observer canonical form for ARX (NX = max(NA, NB) states)
    Matrix<Scalar, NX, NX> A = Matrix<Scalar, NX, NX>::Zero();
    Matrix<Scalar, NX, 1> B = Matrix<Scalar, NX, 1>::Zero();
    Matrix<Scalar, 1, NX> C = Matrix<Scalar, 1, NX>::Zero();
    Matrix<Scalar, 1, 1> D = Matrix<Scalar, 1, 1>::Zero();

    // A: first column = a-coefficients (rows NA..NX-1 stay zero), superdiagonal = 1
    for(std::size_t i = 0; i < NA; ++i)
        A(static_cast<int>(i), 0) = theta(static_cast<int>(i));

    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(static_cast<int>(i), static_cast<int>(i + 1)) = Scalar{1};

    // B: all b-coefficients (rows NB..NX-1 stay zero)
    for(std::size_t i = 0; i < NB; ++i)
        B(static_cast<int>(i), 0) = theta(static_cast<int>(NA + i));

    C(0, 0) = Scalar{1};

    discrete_state_space<Scalar, NX, 1, 1> sys{.A = A, .B = B, .C = C, .D = D};

    // Simulate identified model to compute fit metrics
    Eigen::Matrix<Scalar, static_cast<int>(NX), 1> x = Eigen::Matrix<Scalar, static_cast<int>(NX), 1>::Zero();
    Eigen::VectorX<Scalar> y_predicted(static_cast<Eigen::Index>(N));

    for(std::size_t t = 0; t < N; ++t)
    {
        Eigen::Matrix<Scalar, 1, 1> u_vec;
        u_vec << U(0, static_cast<Eigen::Index>(t));
        auto y_hat = (C * x + D * u_vec).eval();
        y_predicted(static_cast<Eigen::Index>(t)) = y_hat(0, 0);
        x = (A * x + B * u_vec).eval();
    }

    // Flatten Y to column vector for fit metric computation
    Eigen::VectorX<Scalar> y_actual(static_cast<Eigen::Index>(N));
    for(std::size_t t = 0; t < N; ++t)
        y_actual(static_cast<Eigen::Index>(t)) = Y(0, static_cast<Eigen::Index>(t));

    auto metrics = detail::compute_fit_metrics_unchecked(y_actual, y_predicted);

    return arx_result<Scalar, NX, 1, 1>{.system = sys, .metrics = metrics};
}

}

#endif
