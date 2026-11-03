#ifndef HPP_GUARD_CTRLPP_SYSID_DETAIL_BATCH_ARX_DETAIL_H
#define HPP_GUARD_CTRLPP_SYSID_DETAIL_BATCH_ARX_DETAIL_H

// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 4

#include "ctrlpp/types.h"

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/sysid/sysid_result.h"

#include <Eigen/Dense>

#include <cstddef>
#include <optional>
#include <algorithm>

namespace ctrlpp
{

namespace detail
{

template <typename Scalar, std::size_t NP>
struct arx_regression
{
    arx_regression()
        : Phi{}
        , target{}
    {
    }

    Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NP)> Phi;
    Eigen::VectorX<Scalar> target;
};

template <std::size_t NA, std::size_t NB, typename Derived1, typename Derived2>
arx_regression<typename Derived1::Scalar, NA + NB> assemble_arx_regression(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U)
{
    using Scalar = typename Derived1::Scalar;
    static constexpr std::size_t NP = NA + NB;
    static constexpr auto k = static_cast<Eigen::Index>(std::max(NA, NB));

    auto const n_eff = Y.cols() - k;
    arx_regression<Scalar, NP> regression;
    regression.Phi.resize(n_eff, Eigen::NoChange);
    regression.target.resize(n_eff);

    for(Eigen::Index i = 0; i < n_eff; ++i)
    {
        auto const row = k + i;
        for(std::size_t j = 0; j < NA; ++j)
            regression.Phi(i, static_cast<Eigen::Index>(j)) = Y(0, row - static_cast<Eigen::Index>(j + 1));
        for(std::size_t j = 0; j < NB; ++j)
            regression.Phi(i, static_cast<Eigen::Index>(NA + j)) = U(0, row - static_cast<Eigen::Index>(j + 1));
        regression.target(i) = Y(0, row);
    }
    return regression;
}

template <typename Scalar, std::size_t NP>
struct arx_least_squares_fit
{
    arx_least_squares_fit()
        : theta{}
        , diagnostics{}
    {
    }

    Eigen::Matrix<Scalar, static_cast<int>(NP), 1> theta;
    arx_diagnostics<Scalar> diagnostics;
};

// The threshold reaches rank() only: Eigen fixes the number of pivots the solve uses from epsilon during the
// factorization, so the threshold changes what is reported and never what is returned. Left disengaged, Eigen's
// default applies (Eigen/src/QR/ColPivHouseholderQR.h, threshold()).
template <typename Scalar, std::size_t NP>
arx_least_squares_fit<Scalar, NP> solve_arx_least_squares(const arx_regression<Scalar, NP>& regression, std::optional<Scalar> rank_threshold)
{
    Eigen::ColPivHouseholderQR<Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NP)>> qr(regression.Phi);
    if(rank_threshold)
        qr.setThreshold(*rank_threshold);

    arx_least_squares_fit<Scalar, NP> fit{};
    fit.theta = qr.solve(regression.target);
    fit.diagnostics.numerical_rank = static_cast<std::size_t>(qr.rank());
    fit.diagnostics.parameter_count = NP;
    fit.diagnostics.effective_samples = static_cast<std::size_t>(regression.Phi.rows());
    fit.diagnostics.residual_norm = (regression.Phi * fit.theta - regression.target).norm();
    return fit;
}

// Observer canonical realization of dimension max(deg A, deg B) (Ljung 1999, Ch. 4): when NB > NA the extra
// b-coefficients need states of their own, and truncating to NA states would drop them from the transfer function.
template <std::size_t NA, std::size_t NB, typename Scalar>
discrete_state_space<Scalar, std::max(NA, NB), 1, 1> realize_arx_observer_form(const Eigen::Matrix<Scalar, static_cast<int>(NA + NB), 1>& theta)
{
    static constexpr std::size_t NX = std::max(NA, NB);

    Matrix<Scalar, NX, NX> A = Matrix<Scalar, NX, NX>::Zero();
    Matrix<Scalar, NX, 1> B = Matrix<Scalar, NX, 1>::Zero();
    Matrix<Scalar, 1, NX> C = Matrix<Scalar, 1, NX>::Zero();
    Matrix<Scalar, 1, 1> D = Matrix<Scalar, 1, 1>::Zero();

    for(std::size_t i = 0; i < NA; ++i)
        A(static_cast<Eigen::Index>(i), 0) = theta(static_cast<Eigen::Index>(i));

    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i + 1)) = Scalar{1};

    for(std::size_t i = 0; i < NB; ++i)
        B(static_cast<Eigen::Index>(i), 0) = theta(static_cast<Eigen::Index>(NA + i));

    C(0, 0) = Scalar{1};

    return {A, B, C, D};
}

template <typename Scalar, std::size_t NX, typename Derived>
Eigen::VectorX<Scalar> simulate_arx_open_loop(const discrete_state_space<Scalar, NX, 1, 1>& sys, const Eigen::MatrixBase<Derived>& U)
{
    Eigen::Matrix<Scalar, static_cast<int>(NX), 1> x = Eigen::Matrix<Scalar, static_cast<int>(NX), 1>::Zero();
    Eigen::VectorX<Scalar> y_predicted(U.cols());

    for(Eigen::Index t = 0; t < U.cols(); ++t)
    {
        Eigen::Matrix<Scalar, 1, 1> u_vec;
        u_vec << U(0, t);
        auto y_hat = (sys.C * x + sys.D * u_vec).eval();
        y_predicted(t) = y_hat(0, 0);
        x = (sys.A * x + sys.B * u_vec).eval();
    }
    return y_predicted;
}

}

}

#endif
