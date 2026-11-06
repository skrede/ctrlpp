#ifndef HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_QP_COST_H
#define HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_QP_COST_H

#include "ctrlpp/mhe/detail/mhe_qp_layout.h"

#include "ctrlpp/types.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <span>
#include <vector>
#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

template <typename Scalar>
inline void add_upper_block(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index row_off, Eigen::Index col_off, const auto& block, Eigen::Index rows, Eigen::Index cols)
{
    for(Eigen::Index i = 0; i < rows; ++i)
    {
        Eigen::Index start_j = (row_off == col_off) ? i : 0;
        for(Eigen::Index j = start_j; j < cols; ++j)
            if(block(i, j) != Scalar{0})
                trips.emplace_back(static_cast<int>(row_off + i), static_cast<int>(col_off + j), block(i, j));
    }
}

// @cite rao2003 -- Hessian diagonal blocks from arrival cost, process noise, and measurement noise
template <typename Scalar, std::size_t NX>
inline void build_mhe_hessian_diagonal(std::vector<Eigen::Triplet<Scalar>>& trips, std::int32_t Ni, Scalar arrival_weight, const Matrix<Scalar, NX, NX>& P_arr_inv,
        const Matrix<Scalar, NX, NX>& AtQinvA, const Matrix<Scalar, NX, NX>& HtRinvH, const Matrix<Scalar, NX, NX>& Q_inv)
{
    constexpr int nx = static_cast<int>(NX);

    Matrix<Scalar, NX, NX> blk0 = arrival_weight * P_arr_inv + AtQinvA + HtRinvH;
    add_upper_block<Scalar>(trips, 0, 0, blk0, nx, nx);

    Matrix<Scalar, NX, NX> blk_mid = Q_inv + AtQinvA + HtRinvH;
    for(std::int32_t k = 1; k < Ni; ++k)
        add_upper_block<Scalar>(trips, k * nx, k * nx, blk_mid, nx, nx);

    Matrix<Scalar, NX, NX> blk_last = Q_inv + HtRinvH;
    add_upper_block<Scalar>(trips, Ni * nx, Ni * nx, blk_last, nx, nx);
}

// @cite rao2003 -- Hessian cross-term blocks -A^T Q^{-1} between consecutive states
template <typename Scalar, std::size_t NX>
inline void build_mhe_hessian_cross_terms(std::vector<Eigen::Triplet<Scalar>>& trips, std::int32_t Ni, const Matrix<Scalar, NX, NX>& AtQinv)
{
    constexpr int nx = static_cast<int>(NX);
    Matrix<Scalar, NX, NX> cross = (-AtQinv).eval();
    for(std::int32_t k = 0; k < Ni; ++k)
    {
        Eigen::Index row_off = k * nx;
        Eigen::Index col_off = (k + 1) * nx;
        for(Eigen::Index i = 0; i < nx; ++i)
            for(Eigen::Index j = 0; j < nx; ++j)
                if(cross(i, j) != Scalar{0})
                    trips.emplace_back(static_cast<int>(row_off + i), static_cast<int>(col_off + j), cross(i, j));
    }
}

// @cite rawlings2017 -- Rawlings, Mayne & Diehl, MPC, 2nd ed., 2017, Sec. 8.8.1 (block-banded stagewise structure)
template <typename Scalar, std::size_t NX, std::size_t NY>
Eigen::SparseMatrix<Scalar, Eigen::ColMajor> build_mhe_hessian(const mhe_qp_dims& dims, std::int32_t Ni, Scalar arrival_weight, const Matrix<Scalar, NX, NX>& P_arr_inv,
        const Matrix<Scalar, NX, NX>& Q_inv, const Matrix<Scalar, NY, NY>& R_inv, const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NY, NX>& H)
{
    constexpr int nx = static_cast<int>(NX);

    std::vector<Eigen::Triplet<Scalar>> p_trips;
    p_trips.reserve(static_cast<std::size_t>((Ni + 1) * nx * (nx + 1) / 2 * 3 + dims.n_slack));

    Matrix<Scalar, NX, NX> AtQinvA = (A.transpose() * Q_inv * A).eval();
    Matrix<Scalar, NX, NX> AtQinv = (A.transpose() * Q_inv).eval();
    Matrix<Scalar, NX, NX> HtRinvH = (H.transpose() * R_inv * H).eval();

    build_mhe_hessian_diagonal<Scalar, NX>(p_trips, Ni, arrival_weight, P_arr_inv, AtQinvA, HtRinvH, Q_inv);
    build_mhe_hessian_cross_terms<Scalar, NX>(p_trips, Ni, AtQinv);

    if(dims.n_slack > 0)
    {
        constexpr Scalar slack_reg{1e-6};
        Eigen::Index slack_off = dims.n_states;
        for(Eigen::Index i = 0; i < dims.n_slack; ++i)
            p_trips.emplace_back(static_cast<int>(slack_off + i), static_cast<int>(slack_off + i), slack_reg);
    }

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> P(dims.n_dec, dims.n_dec);
    P.setFromTriplets(p_trips.begin(), p_trips.end());
    return P;
}

// @cite rao2003 -- Linear cost from arrival cost, process noise cross-terms, and measurement terms
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
Eigen::VectorX<Scalar> build_mhe_linear_cost(const mhe_qp_dims& dims, std::int32_t Ni, Scalar arrival_weight, const Matrix<Scalar, NX, NX>& P_arr_inv,
        const Matrix<Scalar, NX, NX>& Q_inv, const Matrix<Scalar, NX, NX>& A_lin, const Matrix<Scalar, NX, NU>& B_lin, const Matrix<Scalar, NY, NX>& H_lin,
        const Matrix<Scalar, NY, NY>& R_inv, const Vector<Scalar, NX>& x_arrival, std::span<const Vector<Scalar, NU>> u_buf, std::span<const Vector<Scalar, NY>> z_buf)
{
    constexpr int nx = static_cast<int>(NX);

    Matrix<Scalar, NX, NX> AtQinv = (A_lin.transpose() * Q_inv).eval();
    Matrix<Scalar, NX, NY> HtRinv = (H_lin.transpose() * R_inv).eval();

    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(dims.n_dec);
    q.template segment<nx>(0) = -arrival_weight * P_arr_inv * x_arrival;

    for(std::int32_t k = 0; k < Ni; ++k)
    {
        Vector<Scalar, NX> Buk = B_lin * u_buf[static_cast<std::size_t>(k)];
        q.segment(k * nx, nx) += AtQinv * Buk;
        q.segment((k + 1) * nx, nx) -= Q_inv * Buk;
    }

    for(std::int32_t k = 0; k <= Ni; ++k)
        q.segment(k * nx, nx) -= HtRinv * z_buf[static_cast<std::size_t>(k)];

    return q;
}

}

#endif
