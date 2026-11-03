#ifndef HPP_GUARD_CTRLPP_DETAIL_SCHUR_SWAP_STAGES_H
#define HPP_GUARD_CTRLPP_DETAIL_SCHUR_SWAP_STAGES_H

#include <Eigen/QR>
#include <Eigen/Dense>

#include <cmath>
#include <algorithm>

namespace ctrlpp::detail
{

// K = I kron A11 - A22^T kron I in a zero-padded fixed 4x4; vec(X) stacks columns, so X(i, j) is entry
// j * n1 + i.
template <typename Scalar, int N>
void assemble_sylvester_operator(
        const Eigen::Matrix<Scalar, N, N>& T, Eigen::Index p, Eigen::Index n1, Eigen::Index n2, Eigen::Matrix<Scalar, 4, 4>& K, Eigen::Matrix<Scalar, 4, 1>& rhs)
{
    for(Eigen::Index j = 0; j < n2; ++j)
    {
        for(Eigen::Index i = 0; i < n1; ++i)
        {
            const Eigen::Index row = j * n1 + i;
            for(Eigen::Index k = 0; k < n1; ++k)
                K(row, j * n1 + k) += T(p + i, p + k);
            for(Eigen::Index l = 0; l < n2; ++l)
                K(row, l * n1 + i) -= T(p + n1 + l, p + n1 + j);
            rhs(row) = T(p + i, p + n1 + j);
        }
    }
}

// Column pivoting moves the padding's unused columns to the back, so R(m - 1, m - 1) / R(0, 0) is the
// conditioning signal for every active size m = n1 * n2.
template <typename Scalar>
Scalar sylvester_pivot_ratio(const Eigen::ColPivHouseholderQR<Eigen::Matrix<Scalar, 4, 4>>& qr, Eigen::Index m)
{
    using std::abs;

    const auto R_mat = qr.matrixR();
    const Scalar r_first = abs(R_mat(0, 0));
    const Scalar r_last = abs(R_mat(m - 1, m - 1));
    return (r_first > Scalar{0}) ? (r_last / r_first) : Scalar{0};
}

// Q is the active w x w corner of the Householder Q of [-X; I], padded with the identity.
template <typename Scalar>
void build_swap_basis(const Eigen::Matrix<Scalar, 4, 1>& x_full, Eigen::Index n1, Eigen::Index n2, Eigen::Matrix<Scalar, 4, 4>& Q)
{
    const Eigen::Index w = n1 + n2;
    Eigen::Matrix<Scalar, 2, 2> X = Eigen::Matrix<Scalar, 2, 2>::Zero();
    for(Eigen::Index j = 0; j < n2; ++j)
        for(Eigen::Index i = 0; i < n1; ++i)
            X(i, j) = x_full(j * n1 + i);

    Eigen::Matrix<Scalar, 4, 2> G = Eigen::Matrix<Scalar, 4, 2>::Zero();
    for(Eigen::Index j = 0; j < n2; ++j)
    {
        for(Eigen::Index i = 0; i < n1; ++i)
            G(i, j) = -X(i, j);
        G(n1 + j, j) = Scalar{1};
    }

    Eigen::HouseholderQR<Eigen::Matrix<Scalar, 4, 2>> qr_G(G);
    const Eigen::Matrix<Scalar, 4, 4> Q_full = qr_G.householderQ();
    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
            Q(i, j) = Q_full(i, j);
}

// Explicit loops rather than slice expressions avoid a dynamic-size assignment into a fixed 4x4.
template <typename Scalar>
void transform_active_block(const Eigen::Matrix<Scalar, 4, 4>& Q, const Eigen::Matrix<Scalar, 4, 4>& A_block, Eigen::Index w, Eigen::Matrix<Scalar, 4, 4>& A_new)
{
    Eigen::Matrix<Scalar, 4, 4> QtA = Eigen::Matrix<Scalar, 4, 4>::Zero();
    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
        {
            Scalar acc = Scalar{0};
            for(Eigen::Index k = 0; k < w; ++k)
                acc += Q(k, i) * A_block(k, j);
            QtA(i, j) = acc;
        }
    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
        {
            Scalar acc = Scalar{0};
            for(Eigen::Index k = 0; k < w; ++k)
                acc += QtA(i, k) * Q(k, j);
            A_new(i, j) = acc;
        }
}

template <typename Scalar, int N>
void tentative_swap_block(const Eigen::Matrix<Scalar, N, N>& T, Eigen::Index p, Eigen::Index w, const Eigen::Matrix<Scalar, 4, 4>& Q, Eigen::Matrix<Scalar, 4, 4>& A_new)
{
    Eigen::Matrix<Scalar, 4, 4> A_block = Eigen::Matrix<Scalar, 4, 4>::Zero();
    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
            A_block(i, j) = T(p + i, p + j);
    transform_active_block<Scalar>(Q, A_block, w, A_new);
}

template <typename Scalar>
Scalar swap_residual_norm(const Eigen::Matrix<Scalar, 4, 4>& A_new, Eigen::Index n2, Eigen::Index w)
{
    using std::abs;

    Scalar residual_norm = Scalar{0};
    for(Eigen::Index i = n2; i < w; ++i)
        for(Eigen::Index j = 0; j < n2; ++j)
            residual_norm = std::max(residual_norm, abs(A_new(i, j)));
    return residual_norm;
}

// The new (2, 1) block is written as exact zeros.
template <typename Scalar, int N>
void write_swapped_block(Eigen::Matrix<Scalar, N, N>& T, Eigen::Index p, Eigen::Index n2, Eigen::Index w, const Eigen::Matrix<Scalar, 4, 4>& A_new)
{
    for(Eigen::Index i = 0; i < w; ++i)
        for(Eigen::Index j = 0; j < w; ++j)
            T(p + i, p + j) = A_new(i, j);
    for(Eigen::Index i = n2; i < w; ++i)
        for(Eigen::Index j = 0; j < n2; ++j)
            T(p + i, p + j) = Scalar{0};
}

// The active block already holds Q^T A Q, so only the rows and columns outside it take Q here.
template <typename Scalar, int N>
void apply_swap_outside_block(Eigen::Matrix<Scalar, N, N>& T, Eigen::Index p, Eigen::Index w, const Eigen::Matrix<Scalar, 4, 4>& Q)
{
    for(Eigen::Index j = 0; j < N; ++j)
    {
        if(j >= p && j < p + w)
            continue;
        Eigen::Matrix<Scalar, 4, 1> tmp = Eigen::Matrix<Scalar, 4, 1>::Zero();
        for(Eigen::Index i = 0; i < w; ++i)
            for(Eigen::Index k = 0; k < w; ++k)
                tmp(i) += Q(k, i) * T(p + k, j);
        for(Eigen::Index i = 0; i < w; ++i)
            T(p + i, j) = tmp(i);
    }
    for(Eigen::Index i = 0; i < N; ++i)
    {
        if(i >= p && i < p + w)
            continue;
        Eigen::Matrix<Scalar, 4, 1> tmp = Eigen::Matrix<Scalar, 4, 1>::Zero();
        for(Eigen::Index j = 0; j < w; ++j)
            for(Eigen::Index k = 0; k < w; ++k)
                tmp(j) += T(i, p + k) * Q(k, j);
        for(Eigen::Index j = 0; j < w; ++j)
            T(i, p + j) = tmp(j);
    }
}

template <typename Scalar, int N>
void accumulate_swap_basis(Eigen::Matrix<Scalar, N, N>& U, Eigen::Index p, Eigen::Index w, const Eigen::Matrix<Scalar, 4, 4>& Q)
{
    for(Eigen::Index i = 0; i < N; ++i)
    {
        Eigen::Matrix<Scalar, 4, 1> tmp = Eigen::Matrix<Scalar, 4, 1>::Zero();
        for(Eigen::Index j = 0; j < w; ++j)
            for(Eigen::Index k = 0; k < w; ++k)
                tmp(j) += U(i, p + k) * Q(k, j);
        for(Eigen::Index j = 0; j < w; ++j)
            U(i, p + j) = tmp(j);
    }
}

}

#endif
