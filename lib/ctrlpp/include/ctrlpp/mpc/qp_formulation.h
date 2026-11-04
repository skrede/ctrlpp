#ifndef HPP_GUARD_CTRLPP_MPC_QP_FORMULATION_H
#define HPP_GUARD_CTRLPP_MPC_QP_FORMULATION_H

// @cite rawlings2017 -- Rawlings et al., "Model Predictive Control: Theory, Computation, and Design", 2017, Ch. 2 (sparse condensing)
// @cite borrelli2017 -- Borrelli, Bemporad & Morari, "Predictive Control for Linear and Hybrid Systems", 2017, Ch. 11 (linear MPC, soft-constraint slack penalties)

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/mpc/detail/qp_bounds.h"
#include "ctrlpp/mpc/detail/qp_layout.h"
#include "ctrlpp/mpc/detail/qp_cost_triplets.h"
#include "ctrlpp/mpc/detail/qp_constraint_triplets.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <span>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <optional>

namespace ctrlpp::detail
{

template <ctrlpp_floating_scalar Scalar, std::size_t NX, std::size_t NU>
Eigen::SparseMatrix<Scalar, Eigen::ColMajor> build_cost_matrix(std::int32_t N, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NU, NU>& R, const Matrix<Scalar, NX, NX>& Qf,
        bool has_soft_constraints, Scalar soft_penalty, const std::optional<Vector<Scalar, NX>>& soft_state_penalty = {})
{
    static_assert(NX >= 1 && NU >= 1, "QP formulation requires NX >= 1 and NU >= 1");
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    const auto decision = decision_layout<NX, NU>(N, has_soft_constraints);

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(static_cast<std::size_t>((N + 1) * nx * (nx + 1) / 2 + N * nu * (nu + 1) / 2 + decision.n_slack));

    for(std::int32_t k = 0; k < N; ++k)
        add_diagonal_block<Scalar, NX, NU>(triplets, k * nx, Q, nx);
    add_diagonal_block<Scalar, NX, NU>(triplets, N * nx, Qf, nx);
    for(std::int32_t k = 0; k < N; ++k)
        add_diagonal_block<Scalar, NX, NU>(triplets, decision.n_x_total + k * nu, R, nu);
    if(has_soft_constraints)
        add_soft_constraint_triplets<Scalar, NX, NU>(triplets, N, decision.n_x_total + decision.n_u_total, soft_penalty, soft_state_penalty);

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> P(decision.n_dec, decision.n_dec);
    P.setFromTriplets(triplets.begin(), triplets.end());
    return P;
}

// has_soft_constraints is a parameter rather than read off an optional: present state bounds may still be enforced
// hard, which is the caller's policy choice.
template <typename Scalar, std::size_t NX, std::size_t NU>
Eigen::SparseMatrix<Scalar, Eigen::ColMajor> build_constraint_matrix(std::int32_t N, const Matrix<Scalar, NX, NX>& A_sys, const Matrix<Scalar, NX, NU>& B_sys,
        const std::optional<Vector<Scalar, NX>>& x_min, const std::optional<Vector<Scalar, NX>>& x_max, bool has_soft_constraints, const std::optional<Vector<Scalar, NU>>& u_min,
        const std::optional<Vector<Scalar, NU>>& u_max, const std::optional<Vector<Scalar, NU>>& du_max, const std::optional<terminal_set<Scalar, NX>>& terminal_constraint = {})
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    const auto layout = constraint_layout<Scalar, NX, NU>(N, x_min, x_max, u_min, u_max, du_max, terminal_constraint);
    const auto decision = decision_layout<NX, NU>(N, has_soft_constraints);

    std::vector<Eigen::Triplet<Scalar>> triplets;
    triplets.reserve(static_cast<std::size_t>((N + 1) * nx + N * nx * nx + N * nx * nu + (layout.input_row - layout.state_row) + decision.n_slack +
            (layout.rate_row - layout.input_row) + 2 * (layout.terminal_row - layout.rate_row)));

    Eigen::Index row = 0;
    add_dynamics_triplets<Scalar, NX, NU>(triplets, row, N, decision.n_x_total, A_sys, B_sys);
    if(x_min.has_value() || x_max.has_value())
        add_state_bound_triplets<Scalar, NX, NU>(triplets, row, N, decision.n_x_total, decision.n_u_total, has_soft_constraints);
    if(u_min.has_value() || u_max.has_value())
        add_input_bound_triplets<Scalar, NX, NU>(triplets, row, N, decision.n_x_total);
    if(du_max.has_value())
        add_rate_bound_triplets<Scalar, NX, NU>(triplets, row, N, decision.n_x_total);
    if(terminal_constraint.has_value())
        add_terminal_set_triplets<Scalar, NX>(triplets, row, N * nx, *terminal_constraint);

    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A(layout.n_con, decision.n_dec);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
std::pair<Eigen::VectorX<Scalar>, Eigen::VectorX<Scalar>> build_bounds_vectors(std::int32_t N, const Vector<Scalar, NX>& x0, const std::optional<Vector<Scalar, NX>>& x_min,
        const std::optional<Vector<Scalar, NX>>& x_max, const std::optional<Vector<Scalar, NU>>& u_min, const std::optional<Vector<Scalar, NU>>& u_max,
        const std::optional<Vector<Scalar, NU>>& du_max, const std::optional<terminal_set<Scalar, NX>>& terminal_constraint = {})
{
    const auto layout = constraint_layout<Scalar, NX, NU>(N, x_min, x_max, u_min, u_max, du_max, terminal_constraint);

    Eigen::VectorX<Scalar> l(layout.n_con);
    Eigen::VectorX<Scalar> u(layout.n_con);

    Eigen::Index row = 0;
    set_dynamics_bounds<Scalar, NX>(l, u, row, N, x0);
    if(x_min.has_value() || x_max.has_value())
        set_state_bounds<Scalar, NX>(l, u, row, N, x_min, x_max);
    if(u_min.has_value() || u_max.has_value())
        set_input_bounds<Scalar, NU>(l, u, row, N, u_min, u_max);
    if(du_max.has_value())
        set_rate_bounds<Scalar, NU>(l, u, row, N, *du_max);
    if(terminal_constraint.has_value())
        set_terminal_bounds<Scalar, NX>(l, u, row, *terminal_constraint);

    return {std::move(l), std::move(u)};
}

template <typename Scalar, std::size_t NX, std::size_t NU>
Eigen::VectorX<Scalar> build_cost_vector(std::int32_t /*N*/, Eigen::Index n_dec, const Matrix<Scalar, NX, NX>& /*Q*/, const Matrix<Scalar, NX, NX>& /*Qf*/)
{
    return Eigen::VectorX<Scalar>::Zero(n_dec);
}

template <typename Scalar, std::size_t NX, std::size_t NU>
Eigen::VectorX<Scalar> build_cost_vector(std::int32_t N, Eigen::Index n_dec, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NX, NX>& Qf, const Vector<Scalar, NX>& x_ref)
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(n_dec);
    Vector<Scalar, NX> Qxr = Q * x_ref;
    for(std::int32_t k = 0; k < N; ++k)
        q.segment(k * nx, nx) = -Qxr;
    q.segment(N * nx, nx) = -(Qf * x_ref);
    return q;
}

template <typename Scalar, std::size_t NX, std::size_t NU>
Eigen::VectorX<Scalar> build_cost_vector(
        std::int32_t N, Eigen::Index n_dec, const Matrix<Scalar, NX, NX>& Q, const Matrix<Scalar, NX, NX>& Qf, std::span<const Vector<Scalar, NX>> x_ref)
{
    constexpr int nx = static_cast<int>(NX);
    Eigen::VectorX<Scalar> q = Eigen::VectorX<Scalar>::Zero(n_dec);
    for(std::int32_t k = 0; k < N; ++k)
        q.segment(k * nx, nx) = -(Q * x_ref[static_cast<std::size_t>(k)]);
    q.segment(N * nx, nx) = -(Qf * x_ref[static_cast<std::size_t>(N)]);
    return q;
}

}

#endif
