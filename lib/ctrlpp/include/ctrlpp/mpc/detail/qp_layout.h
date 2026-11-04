#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_QP_LAYOUT_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_QP_LAYOUT_H

#include "ctrlpp/mpc/terminal_set.h"

#include "ctrlpp/types.h"

#include <cstddef>
#include <cstdint>
#include <variant>
#include <optional>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX>
inline Eigen::Index terminal_constraint_rows(const std::optional<terminal_set<Scalar, NX>>& terminal_constraint)
{
    constexpr int nx = static_cast<int>(NX);
    if(!terminal_constraint.has_value())
        return 0;
    return std::visit(
            [](const auto& s) -> Eigen::Index
            {
                using T = std::decay_t<decltype(s)>;
                if constexpr(std::is_same_v<T, ellipsoidal_set<Scalar, NX>>)
                    return nx;
                else
                    return s.H.rows();
            },
            *terminal_constraint);
}

// Every block's presence is read off the optionals here and only here, so the constraint matrix, the bounds vectors
// and the per-solve rate-bound refresh cannot disagree on where a block starts or how tall the QP is.
struct qp_constraint_layout
{
    Eigen::Index state_row;
    Eigen::Index input_row;
    Eigen::Index rate_row;
    Eigen::Index terminal_row;
    Eigen::Index n_con;
};

template <typename Scalar, std::size_t NX, std::size_t NU>
inline qp_constraint_layout constraint_layout(std::int32_t N, const std::optional<Vector<Scalar, NX>>& x_min, const std::optional<Vector<Scalar, NX>>& x_max,
        const std::optional<Vector<Scalar, NU>>& u_min, const std::optional<Vector<Scalar, NU>>& u_max, const std::optional<Vector<Scalar, NU>>& du_max,
        const std::optional<terminal_set<Scalar, NX>>& terminal_constraint)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    const Eigen::Index state_row = (N + 1) * nx;
    const Eigen::Index input_row = state_row + ((x_min.has_value() || x_max.has_value()) ? N * nx : 0);
    const Eigen::Index rate_row = input_row + ((u_min.has_value() || u_max.has_value()) ? N * nu : 0);
    const Eigen::Index terminal_row = rate_row + (du_max.has_value() ? N * nu : 0);

    return {state_row, input_row, rate_row, terminal_row, terminal_row + terminal_constraint_rows<Scalar, NX>(terminal_constraint)};
}

// The cost matrix and the constraint matrix must agree on the decision width, so the slack count is decided here.
template <std::size_t NX>
inline Eigen::Index slack_columns(std::int32_t N, bool has_soft_constraints)
{
    constexpr int nx = static_cast<int>(NX);
    return has_soft_constraints ? N * nx : 0;
}

struct qp_decision_layout
{
    Eigen::Index n_x_total;
    Eigen::Index n_u_total;
    Eigen::Index n_slack;
    Eigen::Index n_dec;
};

template <std::size_t NX, std::size_t NU>
inline qp_decision_layout decision_layout(std::int32_t N, bool has_soft_constraints)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);

    const Eigen::Index n_x_total = (N + 1) * nx;
    const Eigen::Index n_u_total = N * nu;
    const Eigen::Index n_slack = slack_columns<NX>(N, has_soft_constraints);
    return {n_x_total, n_u_total, n_slack, n_x_total + n_u_total + n_slack};
}

}

#endif
