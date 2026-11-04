#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_LAYOUT_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_LAYOUT_H

#include "ctrlpp/mpc/nmpc_config.h"

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

// Decision vector [x_0..x_N, u_0..u_{N-1}, s_path_0..s_path_{N-1}, s_term]; constraint rows [initial state,
// continuity, rate pairs, path, terminal]. Slack exists only for a soft constraint of nonzero dimension.
struct nmpc_layout
{
    static constexpr Eigen::Index x_offset = 0;
    static constexpr Eigen::Index eq_start = 0;

    std::int32_t N;
    bool has_path_slack;
    bool has_term_slack;
    Eigen::Index n_path_slack;
    Eigen::Index n_term_slack;
    Eigen::Index n_vars;
    Eigen::Index u_offset;
    Eigen::Index path_slack_offset;
    Eigen::Index term_slack_offset;
    Eigen::Index n_rate;
    Eigen::Index n_path_con;
    Eigen::Index n_term_con;
    Eigen::Index n_constraints;
    Eigen::Index rate_start;
    Eigen::Index path_con_start;
    Eigen::Index term_con_start;
};

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC>
nmpc_layout make_nmpc_layout(const nmpc_config<Scalar, NX, NU, NC, NTC>& config)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int nu = static_cast<int>(NU);
    constexpr int nc = static_cast<int>(NC);
    constexpr int ntc = static_cast<int>(NTC);
    const std::int32_t N = config.horizon;
    nmpc_layout layout;
    layout.N = N;
    layout.has_path_slack = config.soft_constraints && NC > 0 && config.path_constraint_cb;
    layout.has_term_slack = config.soft_constraints && NTC > 0 && config.terminal_constraint_cb;
    layout.n_path_slack = layout.has_path_slack ? N * nc : 0;
    layout.n_term_slack = layout.has_term_slack ? ntc : 0;
    layout.n_vars = (N + 1) * nx + N * nu + layout.n_path_slack + layout.n_term_slack;
    layout.u_offset = (N + 1) * nx;
    layout.path_slack_offset = layout.u_offset + N * nu;
    layout.term_slack_offset = layout.path_slack_offset + layout.n_path_slack;
    layout.n_rate = config.du_max ? N * nu * 2 : 0;
    layout.n_path_con = config.path_constraint_cb ? N * nc : 0;
    layout.n_term_con = config.terminal_constraint_cb ? ntc : 0;
    layout.n_constraints = (N + 1) * nx + layout.n_rate + layout.n_path_con + layout.n_term_con;
    layout.rate_start = (N + 1) * nx;
    layout.path_con_start = layout.rate_start + layout.n_rate;
    layout.term_con_start = layout.path_con_start + layout.n_path_con;
    return layout;
}

}

#endif
