#ifndef HPP_GUARD_CTRLPP_SYSID_DETAIL_MOESP_STEPS_H
#define HPP_GUARD_CTRLPP_SYSID_DETAIL_MOESP_STEPS_H

// @cite verhaegen2007 -- Verhaegen & Verdult, "Filtering and System Identification: A Least Squares Approach", 2007, Sec. 9.6 (PO-MOESP)

#include "ctrlpp/types.h"

#include "ctrlpp/model/state_space.h"

#include "ctrlpp/sysid/fit_metrics.h"

#include <Eigen/Dense>

#include <cstddef>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Derived>
Eigen::MatrixX<typename Derived::Scalar> block_hankel(const Eigen::MatrixBase<Derived>& signal, Eigen::Index start, Eigen::Index block_rows, Eigen::Index cols)
{
    auto ny = signal.rows();
    Eigen::MatrixX<typename Derived::Scalar> H(ny * block_rows, cols);

    for(Eigen::Index bi = 0; bi < block_rows; ++bi)
        for(Eigen::Index c = 0; c < cols; ++c)
            H.block(bi * ny, c, ny, 1) = signal.col(start + bi + c);
    return H;
}

template <typename Scalar>
struct moesp_lq_result
{
    Eigen::MatrixX<Scalar> O_i;
    Eigen::MatrixX<Scalar> Y_future;
    Eigen::MatrixX<Scalar> U_future;
    Eigen::Index i;
    Eigen::Index j;
};

// @cite verhaegen2007 -- stacked Hankel matrix [U_f; U_p; Y_p; Y_f]
template <typename Scalar, typename Derived1, typename Derived2>
Eigen::MatrixX<Scalar> assemble_moesp_hankel(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U, Eigen::Index i, Eigen::Index j)
{
    auto ny = Y.rows();
    auto nu = U.rows();

    auto Y_past = block_hankel(Y, 0, i, j);
    auto Y_future = block_hankel(Y, i, i, j);
    auto U_past = block_hankel(U, 0, i, j);
    auto U_future = block_hankel(U, i, i, j);

    auto r1 = nu * i;
    auto r2 = (nu + ny) * i;
    auto r3 = ny * i;

    Eigen::MatrixX<Scalar> H(r1 + r2 + r3, j);
    H.topRows(r1) = U_future;
    H.middleRows(r1, nu * i) = U_past;
    H.middleRows(r1 + nu * i, ny * i) = Y_past;
    H.bottomRows(r3) = Y_future;
    return H;
}

// @cite verhaegen2007 -- PO-MOESP projection L32 of the LQ factorization
template <typename Scalar>
Eigen::MatrixX<Scalar> extract_moesp_projection(const Eigen::MatrixX<Scalar>& H, Eigen::Index r1, Eigen::Index r2, Eigen::Index r3, Eigen::Index j)
{
    auto total_rows = r1 + r2 + r3;
    auto r_size = std::min(j, total_rows);
    Eigen::HouseholderQR<Eigen::MatrixX<Scalar>> qr(H.transpose());
    Eigen::MatrixX<Scalar> R_h = qr.matrixQR().topRows(r_size).template triangularView<Eigen::Upper>();
    Eigen::MatrixX<Scalar> L = R_h.transpose();

    auto c1 = std::min(r1, L.cols());
    auto c2 = std::min(r2, L.cols() - c1);
    auto br3 = std::min(r3, L.rows() - r1 - r2);

    if(c2 <= 0 || br3 <= 0)
        return Eigen::MatrixX<Scalar>::Zero(r3, 1);

    return L.block(r1 + r2, c1, br3, c2);
}

// @cite verhaegen2007 -- PO-MOESP projection via LQ decomposition
template <typename Derived1, typename Derived2>
moesp_lq_result<typename Derived1::Scalar> moesp_projection(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U, Eigen::Index i)
{
    using Scalar = typename Derived1::Scalar;

    auto N = Y.cols();
    auto ny = Y.rows();
    auto nu = U.rows();
    Eigen::Index j = N - 2 * i + 1;

    auto r1 = nu * i;
    auto r2 = (nu + ny) * i;
    auto r3 = ny * i;

    auto H = assemble_moesp_hankel<Scalar>(Y, U, i, j);
    auto O_i = extract_moesp_projection<Scalar>(H, r1, r2, r3, j);

    auto Y_future = block_hankel(Y, i, i, j);
    auto U_future = block_hankel(U, i, i, j);

    return {O_i, Y_future, U_future, i, j};
}

// @cite verhaegen2007 -- A from the shift relation Gamma_shifted = Gamma_trunc * A
template <typename Scalar, std::size_t NX>
Matrix<Scalar, NX, NX> extract_system_A(const Eigen::MatrixX<Scalar>& Gamma, Eigen::Index ny)
{
    static constexpr auto nx = static_cast<Eigen::Index>(NX);
    auto gamma_rows = Gamma.rows();
    Eigen::MatrixX<Scalar> Gamma_trunc = Gamma.topRows(gamma_rows - ny);
    Eigen::MatrixX<Scalar> Gamma_shifted = Gamma.bottomRows(gamma_rows - ny);

    Eigen::MatrixX<Scalar> A_dyn = Gamma_trunc.colPivHouseholderQr().solve(Gamma_shifted);

    Matrix<Scalar, NX, NX> A;
    for(Eigen::Index r = 0; r < nx; ++r)
        for(Eigen::Index c = 0; c < nx; ++c)
            A(r, c) = A_dyn(r, c);
    return A;
}

template <typename Scalar, std::size_t NX>
struct moesp_input_response
{
    Matrix<Scalar, NX, 1> B;
    Matrix<Scalar, 1, 1> D;
    Vector<Scalar, NX> x0;
};

// @cite verhaegen2007 -- the initial state is estimated jointly with B and D:
// row t of the regression is [C X_t, u_t, C A^t], with X_{t+1} = A X_t + u_t I.
template <typename Scalar, std::size_t NX, typename Derived1, typename Derived2>
moesp_input_response<Scalar, NX> recover_initial_state_and_BD(
        const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, 1, NX>& C, const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U)
{
    static constexpr auto nx = static_cast<Eigen::Index>(NX);
    auto const samples = Y.cols();
    Eigen::MatrixX<Scalar> Phi(samples, 2 * nx + 1);
    Matrix<Scalar, NX, NX> X = Matrix<Scalar, NX, NX>::Zero();
    Matrix<Scalar, 1, NX> free_response = C;
    for(Eigen::Index t = 0; t < samples; ++t)
    {
        Phi.block(t, 0, 1, nx) = C * X;
        Phi(t, nx) = U(0, t);
        Phi.block(t, nx + 1, 1, nx) = free_response;
        X = (A * X).eval();
        X.diagonal().array() += U(0, t);
        free_response = (free_response * A).eval();
    }
    Eigen::VectorX<Scalar> const theta = Phi.colPivHouseholderQr().solve(Y.row(0).transpose());
    return {theta.head(nx), theta.segment(nx, 1), theta.tail(nx)};
}

template <typename Scalar, std::size_t NX, typename Derived1, typename Derived2>
fit_metrics<Scalar> compute_moesp_metrics(
        const discrete_state_space<Scalar, NX, 1, 1>& sys, const Vector<Scalar, NX>& x0, const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U)
{
    auto const samples = Y.cols();
    Vector<Scalar, NX> x = x0;
    Eigen::VectorX<Scalar> y_predicted(samples);
    for(Eigen::Index t = 0; t < samples; ++t)
    {
        y_predicted(t) = (sys.C * x)(0, 0) + sys.D(0, 0) * U(0, t);
        x = (sys.A * x + sys.B * U(0, t)).eval();
    }
    return compute_fit_metrics_unchecked(Y.row(0).transpose(), y_predicted);
}

}

#endif
