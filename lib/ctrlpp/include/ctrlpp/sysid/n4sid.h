#ifndef HPP_GUARD_CTRLPP_SYSID_N4SID_H
#define HPP_GUARD_CTRLPP_SYSID_N4SID_H

/// @brief N4SID subspace system identification via BDCSVD.
///
/// @cite vanoverschee1994 -- Van Overschee & De Moor, "N4SID: Subspace Algorithms for the Identification of Combined Deterministic-Stochastic Systems", 1994

#include "ctrlpp/types.h"
#include "ctrlpp/model/state_space.h"

#include "ctrlpp/sysid/fit_metrics.h"
#include "ctrlpp/sysid/sysid_result.h"

#include <Eigen/Dense>

#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

namespace detail
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
struct n4sid_lq_result
{
    Eigen::MatrixX<Scalar> O_i;
    Eigen::MatrixX<Scalar> Y_future;
    Eigen::MatrixX<Scalar> U_future;
    Eigen::Index i;
    Eigen::Index j;
};

/// @cite vanoverschee1994 -- Oblique projection via LQ decomposition
template <typename Derived1, typename Derived2>
n4sid_lq_result<typename Derived1::Scalar> n4sid_oblique_projection(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U, Eigen::Index i)
{
    using Scalar = typename Derived1::Scalar;

    auto N = Y.cols();
    auto ny = Y.rows();
    auto nu = U.rows();
    Eigen::Index j = N - 2 * i + 1;

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

    Eigen::HouseholderQR<Eigen::MatrixX<Scalar>> qr(H.transpose());
    auto total_rows = r1 + r2 + r3;
    auto r_size = std::min(j, total_rows);
    Eigen::MatrixX<Scalar> R_h = qr.matrixQR().topRows(r_size).template triangularView<Eigen::Upper>();
    Eigen::MatrixX<Scalar> L = R_h.transpose();

    auto l_rows = L.rows();
    auto l_cols = L.cols();
    auto c1 = std::min(r1, l_cols);
    auto c2 = std::min(r2, l_cols - c1);
    auto br3 = std::min(r3, l_rows - r1 - r2);

    if(c2 <= 0 || br3 <= 0)
        return {Eigen::MatrixX<Scalar>::Zero(r3, 1), Y_future, U_future, i, j};

    Eigen::MatrixX<Scalar> L32 = L.block(r1 + r2, c1, br3, c2);

    return {L32, Y_future, U_future, i, j};
}

/// Extract system matrix A from observability matrix via shift relation.
/// @cite vanoverschee1994 -- Gamma_shifted = Gamma_trunc * A
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

/// Recover B and D via least-squares on state basis propagation.
template <typename Scalar, std::size_t NX, typename Derived1, typename Derived2>
std::pair<Matrix<Scalar, NX, 1>, Matrix<Scalar, 1, 1>> recover_BD(const Matrix<Scalar, NX, NX>& A,
                                                                    const Matrix<Scalar, 1, NX>& C,
                                                                    const Eigen::MatrixBase<Derived1>& Y,
                                                                    const Eigen::MatrixBase<Derived2>& U)
{
    static constexpr auto nx = static_cast<Eigen::Index>(NX);
    auto N_data = Y.cols();

    Eigen::MatrixX<Scalar> Phi(N_data, nx + 1);
    Phi.setZero();
    Eigen::VectorX<Scalar> y_target(N_data);

    Eigen::MatrixX<Scalar> x_basis = Eigen::MatrixX<Scalar>::Zero(nx, nx);

    for(Eigen::Index t = 0; t < N_data; ++t)
    {
        Scalar u_t = U(0, t);
        y_target(t) = Y(0, t);

        for(Eigen::Index j = 0; j < nx; ++j)
            Phi(t, j) = (C * x_basis.col(j))(0, 0);
        Phi(t, nx) = u_t;

        Eigen::MatrixX<Scalar> x_basis_new = A * x_basis;
        for(Eigen::Index j = 0; j < nx; ++j)
            x_basis_new(j, j) += u_t;
        x_basis = x_basis_new;
    }

    Eigen::VectorX<Scalar> bd = Phi.colPivHouseholderQr().solve(y_target);

    Matrix<Scalar, NX, 1> B;
    for(Eigen::Index r = 0; r < nx; ++r)
        B(r, 0) = bd(r);

    Matrix<Scalar, 1, 1> D;
    D(0, 0) = bd(nx);

    return {B, D};
}

/// Simulate identified model and compute fit metrics.
template <typename Scalar, std::size_t NX, typename Derived1, typename Derived2>
auto compute_n4sid_metrics(const Matrix<Scalar, NX, NX>& A,
                           const Matrix<Scalar, NX, 1>& B,
                           const Matrix<Scalar, 1, NX>& C,
                           const Matrix<Scalar, 1, 1>& D,
                           const Eigen::MatrixBase<Derived1>& Y,
                           const Eigen::MatrixBase<Derived2>& U) -> fit_metrics<Scalar>
{
    auto N_data = Y.cols();
    Vector<Scalar, NX> x = Vector<Scalar, NX>::Zero();
    Eigen::VectorX<Scalar> y_predicted(N_data);
    Eigen::VectorX<Scalar> y_actual(N_data);

    for(Eigen::Index t = 0; t < N_data; ++t)
    {
        Matrix<Scalar, 1, 1> u_vec;
        u_vec(0, 0) = U(0, t);
        auto y_hat = (C * x + D * u_vec).eval();
        y_predicted(t) = y_hat(0, 0);
        x = (A * x + B * u_vec).eval();
        y_actual(t) = Y(0, t);
    }

    return compute_fit_metrics(y_actual, y_predicted);
}

}

template <typename Derived1, typename Derived2>
Eigen::VectorX<typename Derived1::Scalar> n4sid_singular_values(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U, std::size_t block_rows = 0)
{
    using Scalar = typename Derived1::Scalar;

    auto N = Y.cols();
    auto i = static_cast<Eigen::Index>(block_rows);
    if(i == 0)
        i = std::min(static_cast<Eigen::Index>(N / 4), Eigen::Index{30});

    Eigen::Index j = N - 2 * i + 1;
    if(j <= 0)
        return Eigen::VectorX<Scalar>{};

    auto [O_i, Y_f, U_f, bi, bj] = detail::n4sid_oblique_projection(Y, U, i);

    Eigen::BDCSVD<Eigen::MatrixX<Scalar>> svd(O_i, Eigen::ComputeThinU | Eigen::ComputeThinV);

    return svd.singularValues();
}

template <std::size_t NX, typename Derived1, typename Derived2>
n4sid_result<typename Derived1::Scalar, NX, 1, 1> n4sid(const Eigen::MatrixBase<Derived1>& Y, const Eigen::MatrixBase<Derived2>& U, std::size_t block_rows = 0)

{
    using Scalar = typename Derived1::Scalar;
    static constexpr auto nx = static_cast<Eigen::Index>(NX);

    auto N = Y.cols();
    auto ny = Y.rows();

    auto i = static_cast<Eigen::Index>(block_rows);
    if(i == 0)
        i = std::min(static_cast<Eigen::Index>(N / 4), Eigen::Index{30});

    auto [O_i, Y_future, U_future, bi, bj] = detail::n4sid_oblique_projection(Y, U, i);

    Eigen::BDCSVD<Eigen::MatrixX<Scalar>> svd(O_i, Eigen::ComputeThinU | Eigen::ComputeThinV);
    auto sv = svd.singularValues();
    auto U_svd = svd.matrixU();

    Scalar cond = (sv.size() >= nx && sv(nx - 1) > Scalar{0}) ? sv(0) / sv(nx - 1) : std::numeric_limits<Scalar>::infinity();

    // Build observability matrix Gamma from SVD
    Eigen::Index rank = std::min(static_cast<Eigen::Index>(sv.size()), nx);
    Eigen::MatrixX<Scalar> Gamma(U_svd.rows(), nx);
    Gamma.setZero();
    for(Eigen::Index c = 0; c < rank; ++c)
        Gamma.col(c) = U_svd.col(c) * std::sqrt(sv(c));

    // Extract C from first ny rows
    Matrix<Scalar, 1, NX> C;
    for(Eigen::Index c = 0; c < nx; ++c)
        C(0, c) = Gamma(c);

    // Extract A via shift relation on observability matrix
    auto A = detail::extract_system_A<Scalar, NX>(Gamma, ny);

    // Recover B and D via least-squares
    auto [B, D] = detail::recover_BD<Scalar, NX>(A, C, Y, U);

    discrete_state_space<Scalar, NX, 1, 1> sys{.A = A, .B = B, .C = C, .D = D};

    auto metrics = detail::compute_n4sid_metrics<Scalar, NX>(A, B, C, D, Y, U);

    return {.system = sys, .singular_values = sv, .metrics = metrics, .condition_number = cond};
}

}

#endif
