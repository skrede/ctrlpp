#ifndef HPP_GUARD_CTRLPP_MHE_MHE_QP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MHE_MHE_QP_FORMULATION_H

// @cite rao2003 -- Rao et al., "Constrained State Estimation for Nonlinear Discrete-Time Systems", 2003
// @cite rawlings2017 -- Rawlings, Mayne & Diehl, "Model Predictive Control: Theory, Computation, and Design", 2nd ed., 2017, Ch. 4 (State Estimation)
// @cite kuhl2011 -- Kuhl et al., "A Real-Time Algorithm for Moving Horizon State and Parameter Estimation", 2011

#include "ctrlpp/mhe/detail/mhe_qp_cost.h"
#include "ctrlpp/mhe/detail/mhe_qp_bounds.h"
#include "ctrlpp/mhe/detail/mhe_qp_layout.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/qp_types.h"

#include "ctrlpp/util/concepts.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <span>
#include <array>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <optional>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NY>
Eigen::SparseMatrix<Scalar, Eigen::ColMajor> build_mhe_constraint_matrix(const mhe_qp_dims& dims, std::int32_t Ni, const Matrix<Scalar, NX, NX>& A, const Matrix<Scalar, NY, NX>& H,
        bool has_box_bounds, bool has_soft_constraints, bool has_residual_bounds)
{
    constexpr int nx = static_cast<int>(NX);
    (void)A;
    std::vector<Eigen::Triplet<Scalar>> a_trips;
    a_trips.reserve(static_cast<std::size_t>(dims.n_box * 2 + dims.n_residual * nx));

    Eigen::Index row = 0;
    if(has_box_bounds)
        add_mhe_box_triplets<Scalar, NX>(a_trips, row, Ni, dims.n_states, has_soft_constraints);
    if(has_residual_bounds)
        add_mhe_residual_triplets<Scalar, NX, NY>(a_trips, row, Ni, H);

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A_con(dims.n_con, dims.n_dec);
    A_con.setFromTriplets(a_trips.begin(), a_trips.end());
    return A_con;
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
qp_problem<Scalar> build_mhe_qp_structure(std::size_t N, Scalar arrival_weight, const Matrix<Scalar, NX, NX>& P_arr_inv, const Matrix<Scalar, NX, NX>& Q_inv,
        const Matrix<Scalar, NY, NY>& R_inv, const std::array<Matrix<Scalar, NX, NX>, 1>& A_lin, const std::array<Matrix<Scalar, NY, NX>, 1>& H_lin, bool has_box_bounds,
        bool has_soft_constraints, Scalar soft_penalty, bool has_residual_bounds)
{
    const std::int32_t Ni = static_cast<std::int32_t>(N);
    auto dims = compute_mhe_dims<NX, NY>(N, has_box_bounds, has_soft_constraints, has_residual_bounds);

    auto P = build_mhe_hessian<Scalar, NX, NY>(dims, Ni, arrival_weight, P_arr_inv, Q_inv, R_inv, A_lin[0], H_lin[0]);
    auto A_con = build_mhe_constraint_matrix<Scalar, NX, NY>(dims, Ni, A_lin[0], H_lin[0], has_box_bounds, has_soft_constraints, has_residual_bounds);

    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(dims.n_dec);
    if(dims.n_slack > 0)
    {
        Eigen::Index slack_off = dims.n_states;
        for(Eigen::Index i = 0; i < dims.n_slack; ++i)
            q(slack_off + i) = soft_penalty;
    }

    Eigen::VectorX<Scalar> l = Eigen::VectorX<Scalar>::Zero(dims.n_con);
    Eigen::VectorX<Scalar> u = Eigen::VectorX<Scalar>::Zero(dims.n_con);

    return {std::move(P), std::move(q), std::move(A_con), std::move(l), std::move(u)};
}

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
std::pair<Eigen::VectorX<Scalar>, Eigen::VectorX<Scalar>> build_mhe_update_bounds(const mhe_qp_dims& dims, std::int32_t Ni, const Matrix<Scalar, NX, NU>& B_lin,
        std::span<const Vector<Scalar, NU>> u_buf, std::span<const Vector<Scalar, NY>> z_buf, bool has_box_bounds, const std::optional<Vector<Scalar, NX>>& x_min,
        const std::optional<Vector<Scalar, NX>>& x_max, const std::optional<Vector<Scalar, NY>>& residual_bound)
{
    constexpr auto inf = std::numeric_limits<Scalar>::infinity();
    (void)B_lin;
    (void)u_buf;

    Eigen::VectorX<Scalar> l(dims.n_con);
    Eigen::VectorX<Scalar> u(dims.n_con);
    Eigen::Index row = 0;

    if(has_box_bounds)
    {
        const bool need_slack = dims.n_slack > 0;
        Vector<Scalar, NX> lb = x_min.value_or(Vector<Scalar, NX>::Constant(-inf));
        Vector<Scalar, NX> ub = x_max.value_or(Vector<Scalar, NX>::Constant(inf));

        if(!need_slack)
            set_mhe_hard_box_bounds<Scalar, NX>(l, u, row, Ni, lb, ub);
        else
            set_mhe_soft_box_bounds<Scalar, NX>(l, u, row, Ni, lb, ub);
    }

    if(residual_bound.has_value())
        set_mhe_residual_bounds<Scalar, NY>(l, u, row, Ni, z_buf, *residual_bound);

    return {std::move(l), std::move(u)};
}

// @cite kuhl2011 -- Kuhl et al., "A Real-Time Algorithm for Moving Horizon State and Parameter Estimation", 2011 (real-time MHE reassembly with arrival-cost update)
template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU, std::size_t NY>
qp_update<Scalar> build_mhe_qp_update(std::size_t N, Scalar arrival_weight, const Matrix<Scalar, NX, NX>& P_arr_inv, const Matrix<Scalar, NX, NX>& Q_inv,
        const Matrix<Scalar, NY, NY>& R_inv, const Matrix<Scalar, NX, NX>& A_lin, const Matrix<Scalar, NX, NU>& B_lin, const Matrix<Scalar, NY, NX>& H_lin,
        const Vector<Scalar, NX>& x_arrival, std::span<const Vector<Scalar, NU>> u_buf, std::span<const Vector<Scalar, NY>> z_buf, bool has_box_bounds, bool has_soft_constraints,
        Scalar soft_penalty, const std::optional<Vector<Scalar, NX>>& x_min, const std::optional<Vector<Scalar, NX>>& x_max, const std::optional<Vector<Scalar, NY>>& residual_bound,
        const Eigen::VectorX<Scalar>& warm_x, const Eigen::VectorX<Scalar>& warm_y)
{
    const std::int32_t Ni = static_cast<std::int32_t>(N);
    (void)soft_penalty;
    // The residual block's height and its rows are decided by the same optional, so the update vectors cannot be sized
    // for a block the bounds loop declines to fill.
    auto dims = compute_mhe_dims<NX, NY>(N, has_box_bounds, has_soft_constraints, residual_bound.has_value());

    auto q = build_mhe_linear_cost<Scalar, NX, NU, NY>(dims, Ni, arrival_weight, P_arr_inv, Q_inv, A_lin, B_lin, H_lin, R_inv, x_arrival, u_buf, z_buf);
    auto [l, u] = build_mhe_update_bounds<Scalar, NX, NU, NY>(dims, Ni, B_lin, u_buf, z_buf, has_box_bounds, x_min, x_max, residual_bound);

    return {std::move(q), std::move(l), std::move(u), warm_x, warm_y};
}

}

#endif
