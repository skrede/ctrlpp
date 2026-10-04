#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_QP_TERMINAL_RICCATI_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_QP_TERMINAL_RICCATI_H

#include "qp/sparse_mpc_sweep.h"
#include "qp/mpc_tracking_cells.h"

#include <Eigen/Dense>

#include <cstdint>

namespace ctrlpp::bench::problems::qp
{

struct terminal_riccati_settings
{
    double rho;
    std::int32_t max_iterations;
    double gain_tolerance;
};

inline constexpr terminal_riccati_settings terminal_riccati{1.0, 1000, 1e-5};

/// Port of the stationary Riccati iteration in TinyMPC's tiny_precompute_and_set_cache
/// (src/tinympc/tiny_api.cpp at commit 023f36bd5b27267e1a96ed47a774d120b5707f13), with the same
/// starting point, stopping rule and operation order, applied to Q + rho I and R + rho I.
inline auto stationary_cost_to_go(const mpc_sweep_plant& plant, const Eigen::MatrixXd& Q, const Eigen::MatrixXd& R)
    -> Eigen::MatrixXd
{
    const Eigen::MatrixXd& A = plant.transition;
    const Eigen::MatrixXd& B = plant.actuation;
    const Eigen::Index nx = A.rows();
    const Eigen::Index nu = B.cols();
    const Eigen::MatrixXd Q1 = Q + terminal_riccati.rho * Eigen::MatrixXd::Identity(nx, nx);
    const Eigen::MatrixXd R1 = R + terminal_riccati.rho * Eigen::MatrixXd::Identity(nu, nu);
    Eigen::MatrixXd K_next = Eigen::MatrixXd::Zero(nu, nx);
    Eigen::MatrixXd P_next = terminal_riccati.rho * Eigen::MatrixXd::Identity(nx, nx);
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(nx, nx);
    for(std::int32_t i = 0; i < terminal_riccati.max_iterations; ++i)
    {
        const Eigen::MatrixXd K = (R1 + B.transpose() * P_next * B).inverse() * B.transpose() * P_next * A;
        P = Q1 + A.transpose() * P_next * (A - B * K);
        if((K - K_next).cwiseAbs().maxCoeff() < terminal_riccati.gain_tolerance)
            break;
        K_next = K;
        P_next = P;
    }
    return P;
}

inline auto terminal_weight(const mpc_sweep_cell& cell, const mpc_sweep_plant& plant) -> Eigen::MatrixXd
{
    const Eigen::MatrixXd Q = stage_state_weight * Eigen::MatrixXd::Identity(cell.nx, cell.nx);
    if(!is_long_horizon(cell))
        return legacy_terminal_weight * Eigen::MatrixXd::Identity(cell.nx, cell.nx);
    return stationary_cost_to_go(plant, Q, stage_input_weight * Eigen::MatrixXd::Identity(cell.nu, cell.nu));
}

}

#endif
