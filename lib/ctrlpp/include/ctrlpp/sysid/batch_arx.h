#ifndef HPP_GUARD_CTRLPP_SYSID_BATCH_ARX_H
#define HPP_GUARD_CTRLPP_SYSID_BATCH_ARX_H

/// @brief Batch ARX model identification via QR decomposition.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 4

#include "ctrlpp/types.h"
#include "ctrlpp/model/state_space.h"

#include "ctrlpp/sysid/fit_metrics.h"
#include "ctrlpp/sysid/sysid_result.h"

#include <Eigen/Dense>

#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

template <std::size_t NA, std::size_t NB, typename Derived1, typename Derived2>
arx_result<typename Derived1::Scalar, NA, 1, 1> batch_arx(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U)
{
    using Scalar = typename Derived1::Scalar;
    static constexpr std::size_t NP = NA + NB;

    auto N = static_cast<std::size_t>(Y.cols());
    std::size_t k = std::max(NA, NB);
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

    // Observer canonical form for ARX
    Matrix<Scalar, NA, NA> A = Matrix<Scalar, NA, NA>::Zero();
    Matrix<Scalar, NA, 1> B = Matrix<Scalar, NA, 1>::Zero();
    Matrix<Scalar, 1, NA> C = Matrix<Scalar, 1, NA>::Zero();
    Matrix<Scalar, 1, 1> D = Matrix<Scalar, 1, 1>::Zero();

    // A: first column = a-coefficients, superdiagonal = 1
    for(std::size_t i = 0; i < NA; ++i)
        A(static_cast<int>(i), 0) = theta(static_cast<int>(i));

    for(std::size_t i = 0; i + 1 < NA; ++i)
        A(static_cast<int>(i), static_cast<int>(i + 1)) = Scalar{1};

    // B: b-coefficients, zero-padded if NB < NA
    for(std::size_t i = 0; i < NB && i < NA; ++i)
        B(static_cast<int>(i), 0) = theta(static_cast<int>(NA + i));

    C(0, 0) = Scalar{1};

    discrete_state_space<Scalar, NA, 1, 1> sys{.A = A, .B = B, .C = C, .D = D};

    // Simulate identified model to compute fit metrics
    Eigen::Matrix<Scalar, static_cast<int>(NA), 1> x = Eigen::Matrix<Scalar, static_cast<int>(NA), 1>::Zero();
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

    auto metrics = compute_fit_metrics(y_actual, y_predicted);

    return {.system = sys, .metrics = metrics};
}

} // namespace ctrlpp

#endif
