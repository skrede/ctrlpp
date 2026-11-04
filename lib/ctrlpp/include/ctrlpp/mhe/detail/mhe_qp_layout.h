#ifndef HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_QP_LAYOUT_H
#define HPP_GUARD_CTRLPP_MHE_DETAIL_MHE_QP_LAYOUT_H

#include <Eigen/Core>

#include <cstddef>
#include <cstdint>

namespace ctrlpp::detail
{

// z = [x_0..x_N, s_lo_0..s_lo_N, s_up_0..s_up_N], the slacks only for a softened box. The dynamics enter through
// the process-noise cost, not an equality block, which would force every residual to zero and void the weighting.
struct mhe_qp_dims
{
    Eigen::Index n_states;
    Eigen::Index n_slack;
    Eigen::Index n_dec;
    Eigen::Index n_dyn;
    Eigen::Index n_box;
    Eigen::Index n_residual;
    Eigen::Index n_con;
};

template <std::size_t NX, std::size_t NY>
constexpr mhe_qp_dims compute_mhe_dims(std::size_t N, bool has_box_bounds, bool has_soft_constraints, bool has_residual_bounds)
{
    constexpr int nx = static_cast<int>(NX);
    constexpr int ny = static_cast<int>(NY);
    const std::int32_t Ni = static_cast<std::int32_t>(N);

    Eigen::Index n_state_vars = (Ni + 1) * nx;
    bool need_slack = has_soft_constraints && has_box_bounds;
    Eigen::Index n_slack = need_slack ? 2 * n_state_vars : 0;
    Eigen::Index n_dec = n_state_vars + n_slack;

    Eigen::Index n_dyn = 0;
    Eigen::Index n_box = has_box_bounds ? (need_slack ? 4 * n_state_vars : n_state_vars) : 0;
    Eigen::Index n_residual = has_residual_bounds ? (Ni + 1) * ny : 0;
    Eigen::Index n_con = n_dyn + n_box + n_residual;

    return {n_state_vars, n_slack, n_dec, n_dyn, n_box, n_residual, n_con};
}

}

#endif
