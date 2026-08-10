#ifndef HPP_GUARD_CTRLPP_SYSID_DETAIL_BATCH_ARX_DETAIL_H
#define HPP_GUARD_CTRLPP_SYSID_DETAIL_BATCH_ARX_DETAIL_H

/// @brief Stages of batch ARX identification: regressor assembly, observer
/// canonical realization and open-loop simulation.
///
/// @cite ljung1999 -- Ljung, "System Identification: Theory for the User", 1999, Ch. 4

#include "ctrlpp/types.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

namespace detail
{

template <typename Scalar, std::size_t NP>
struct arx_regression
{
    Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NP)> Phi{};
    Eigen::VectorX<Scalar> target{};
};

template <std::size_t NA, std::size_t NB, typename Derived1, typename Derived2>
auto assemble_arx_regression(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U)
    -> arx_regression<typename Derived1::Scalar, NA + NB>
{
    using Scalar = typename Derived1::Scalar;
    static constexpr std::size_t NP = NA + NB;
    static constexpr auto k = static_cast<Eigen::Index>(std::max(NA, NB));

    auto const n_eff = Y.cols() - k;
    arx_regression<Scalar, NP> regression{
        Eigen::Matrix<Scalar, Eigen::Dynamic, static_cast<int>(NP)>(n_eff, static_cast<int>(NP)),
        Eigen::VectorX<Scalar>(n_eff)};

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

/// Observer canonical realization dimension = max(deg A, deg B) (Ljung 1999, Ch. 4).
/// When NB > NA the extra b-coefficients b_{NA+1..NB} require additional states;
/// truncating to NA states would silently drop them from the transfer function.
template <std::size_t NA, std::size_t NB, typename Scalar>
auto realize_arx_observer_form(const Eigen::Matrix<Scalar, static_cast<int>(NA + NB), 1>& theta)
    -> discrete_state_space<Scalar, std::max(NA, NB), 1, 1>
{
    static constexpr std::size_t NX = std::max(NA, NB);

    Matrix<Scalar, NX, NX> A = Matrix<Scalar, NX, NX>::Zero();
    Matrix<Scalar, NX, 1> B = Matrix<Scalar, NX, 1>::Zero();
    Matrix<Scalar, 1, NX> C = Matrix<Scalar, 1, NX>::Zero();
    Matrix<Scalar, 1, 1> D = Matrix<Scalar, 1, 1>::Zero();

    for(std::size_t i = 0; i < NA; ++i)
        A(static_cast<int>(i), 0) = theta(static_cast<int>(i));

    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(static_cast<int>(i), static_cast<int>(i + 1)) = Scalar{1};

    for(std::size_t i = 0; i < NB; ++i)
        B(static_cast<int>(i), 0) = theta(static_cast<int>(NA + i));

    C(0, 0) = Scalar{1};

    return {.A = A, .B = B, .C = C, .D = D};
}

/// The state starts at zero, so the first max(deg A, deg B) predictions are
/// produced from a state the identification never estimated.
template <typename Scalar, std::size_t NX, typename Derived>
auto simulate_arx_open_loop(const discrete_state_space<Scalar, NX, 1, 1>& sys, const Eigen::MatrixBase<Derived>& U)
    -> Eigen::VectorX<Scalar>
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
