#ifndef HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_QP_BOUNDS_H
#define HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_QP_BOUNDS_H

#include "ctrlpp/types.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <span>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
inline void add_mhe_identity_rows(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, std::int32_t Ni, Eigen::Index col_off)
{
    constexpr int nx = static_cast<int>(NX);
    for(std::int32_t k = 0; k <= Ni; ++k)
    {
        Eigen::Index sk_off = col_off + k * nx;
        for(Eigen::Index i = 0; i < nx; ++i)
            trips.emplace_back(static_cast<int>(row + i), static_cast<int>(sk_off + i), Scalar{1});
        row += nx;
    }
}

template <typename Scalar, std::size_t NX>
inline void add_mhe_softened_rows(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, std::int32_t Ni, Eigen::Index slack_off, Scalar slack_sign)
{
    constexpr int nx = static_cast<int>(NX);
    for(std::int32_t k = 0; k <= Ni; ++k)
    {
        Eigen::Index xk_off = k * nx;
        Eigen::Index sk_off = slack_off + k * nx;
        for(Eigen::Index i = 0; i < nx; ++i)
        {
            trips.emplace_back(static_cast<int>(row + i), static_cast<int>(xk_off + i), Scalar{1});
            trips.emplace_back(static_cast<int>(row + i), static_cast<int>(sk_off + i), slack_sign);
        }
        row += nx;
    }
}

// A softened box is four row blocks, x + s_lo >= lb, x - s_up <= ub, s_lo >= 0 and s_up >= 0: one signed slack cannot
// soften both sides, and without the non-negativity rows the optimizer is rewarded for a negative slack.
template <typename Scalar, std::size_t NX>
inline void add_mhe_box_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, std::int32_t Ni, Eigen::Index n_states, bool need_slack)
{
    constexpr int nx = static_cast<int>(NX);

    if(!need_slack)
    {
        add_mhe_identity_rows<Scalar, NX>(trips, row, Ni, 0);
        return;
    }

    Eigen::Index n_state_vars = (Ni + 1) * nx;
    Eigen::Index s_lo_off = n_states;
    Eigen::Index s_up_off = n_states + n_state_vars;

    add_mhe_softened_rows<Scalar, NX>(trips, row, Ni, s_lo_off, Scalar{1});
    add_mhe_softened_rows<Scalar, NX>(trips, row, Ni, s_up_off, Scalar{-1});
    add_mhe_identity_rows<Scalar, NX>(trips, row, Ni, s_lo_off);
    add_mhe_identity_rows<Scalar, NX>(trips, row, Ni, s_up_off);
}

template <typename Scalar, std::size_t NX, std::size_t NY>
inline void add_mhe_residual_triplets(std::vector<Eigen::Triplet<Scalar>>& trips, Eigen::Index& row, std::int32_t Ni, const Matrix<Scalar, NY, NX>& H)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    for(std::int32_t k = 0; k <= Ni; ++k)
    {
        Eigen::Index xk_off = k * nx;
        for(Eigen::Index i = 0; i < ny; ++i)
            for(Eigen::Index j = 0; j < nx; ++j)
                if(H(i, j) != Scalar{0})
                    trips.emplace_back(static_cast<int>(row + i), static_cast<int>(xk_off + j), H(i, j));
        row += ny;
    }
}

template <typename Scalar, std::size_t NX>
inline void set_mhe_hard_box_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, std::int32_t Ni, const Vector<Scalar, NX>& lb, const Vector<Scalar, NX>& ub)
{
    constexpr int nx = static_cast<int>(NX);
    for(std::int32_t k = 0; k <= Ni; ++k)
    {
        l.segment(row, nx) = lb;
        u.segment(row, nx) = ub;
        row += nx;
    }
}

template <typename Scalar, std::size_t NX>
inline void set_mhe_soft_box_bounds(Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, std::int32_t Ni, const Vector<Scalar, NX>& lb, const Vector<Scalar, NX>& ub)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr auto inf = std::numeric_limits<Scalar>::infinity();
    for(std::int32_t k = 0; k <= Ni; ++k)
    {
        l.segment(row, nx) = lb;
        u.segment(row, nx).setConstant(inf);
        row += nx;
    }
    for(std::int32_t k = 0; k <= Ni; ++k)
    {
        l.segment(row, nx).setConstant(-inf);
        u.segment(row, nx) = ub;
        row += nx;
    }
    for(std::int32_t k = 0; k < 2 * (Ni + 1); ++k)
    {
        l.segment(row, nx).setZero();
        u.segment(row, nx).setConstant(inf);
        row += nx;
    }
}

template <typename Scalar, std::size_t NY>
inline void set_mhe_residual_bounds(
        Eigen::VectorX<Scalar>& l, Eigen::VectorX<Scalar>& u, Eigen::Index& row, std::int32_t Ni, std::span<const Vector<Scalar, NY>> z_buf, const Vector<Scalar, NY>& thresh)
{
    constexpr int ny = static_cast<int>(NY);
    for(std::int32_t k = 0; k <= Ni; ++k)
    {
        const auto& zk = z_buf[static_cast<std::size_t>(k)];
        for(Eigen::Index i = 0; i < ny; ++i)
        {
            l(row + i) = zk(i) - thresh(i);
            u(row + i) = zk(i) + thresh(i);
        }
        row += ny;
    }
}

}

#endif
