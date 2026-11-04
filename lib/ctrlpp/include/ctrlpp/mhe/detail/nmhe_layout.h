#ifndef HPP_GUARD_CTRLPP_MHE_DETAIL_NMHE_LAYOUT_H
#define HPP_GUARD_CTRLPP_MHE_DETAIL_NMHE_LAYOUT_H

#include "ctrlpp/mhe/mhe_config.h"

#include <Eigen/Core>

#include <cstddef>

namespace ctrlpp::detail
{

// The decision vector is [x_0 .. x_N, s_0 .. s_N] and the constraint rows are [path, residual pairs]; the path slack s_k,
// NC entries per stage, exists only for a soft path constraint.
template <std::size_t NX, std::size_t N>
struct nmhe_layout
{
    static constexpr Eigen::Index x_offset = 0;
    static constexpr Eigen::Index n_states = static_cast<Eigen::Index>((N + 1) * NX);
    static constexpr Eigen::Index path_slack_offset = n_states;
    static constexpr Eigen::Index path_con_start = 0;

    bool has_path_constraint;
    bool has_path_slack;
    bool has_residual;
    Eigen::Index n_path_slack;
    Eigen::Index n_vars;
    Eigen::Index n_path_con;
    Eigen::Index n_residual_con;
    Eigen::Index n_constraints;
    Eigen::Index residual_con_start;
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY, std::size_t N, std::size_t NC>
nmhe_layout<NX, N> make_nmhe_layout(const nmhe_config<Scalar, NX, NU, NY, N, NC>& config)
{
    constexpr int ny = static_cast<int>(NY);
    constexpr int nc = static_cast<int>(NC);
    constexpr int Ni = static_cast<int>(N);
    nmhe_layout<NX, N> layout;
    layout.has_path_constraint = NC > 0 && config.path_constraint_cb.has_value();
    layout.has_path_slack = config.soft_constraints && layout.has_path_constraint;
    layout.has_residual = config.residual_bound.has_value();
    layout.n_path_slack = layout.has_path_slack ? (Ni + 1) * nc : 0;
    layout.n_vars = layout.n_states + layout.n_path_slack;
    layout.n_path_con = layout.has_path_constraint ? (Ni + 1) * nc : 0;
    layout.n_residual_con = layout.has_residual ? (Ni + 1) * ny * 2 : 0;
    layout.n_constraints = layout.n_path_con + layout.n_residual_con;
    layout.residual_con_start = layout.path_con_start + layout.n_path_con;
    return layout;
}

}

#endif
