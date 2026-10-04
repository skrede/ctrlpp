#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_QP_MPC_CONDENSING_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_QP_MPC_CONDENSING_H

#include "qp/sparse_mpc_sweep.h"
#include "qp/mpc_tracking_problem.h"

#include <Eigen/Dense>

#include <cstdint>

namespace ctrlpp::bench::problems::qp
{

/// The tracking problem over U = [u_0..u_{N-1}] alone, with X = [x_1..x_N] = Phi x0 + Gamma U substituted.
/// Only the gradient and the state rows' bounds depend on x0; everything else is posed once.
struct condensed_tracking
{
    Eigen::MatrixXd phi;
    Eigen::MatrixXd gamma;
    Eigen::MatrixXd gradient_map;
    Eigen::VectorXd stacked_target;
    Eigen::MatrixXd actuated_phi;
    Eigen::MatrixXd P;
    Eigen::MatrixXd A;
    Eigen::VectorXd l;
    Eigen::VectorXd u;
};

inline auto lifted_free_response(const mpc_sweep_cell& cell, const mpc_sweep_plant& plant) -> Eigen::MatrixXd
{
    Eigen::MatrixXd phi(cell.horizon * cell.nx, cell.nx);
    Eigen::MatrixXd power = plant.transition;
    for(std::int32_t k = 0; k < cell.horizon; ++k)
    {
        phi.middleRows(k * cell.nx, cell.nx) = power;
        power = plant.transition * power;
    }
    return phi;
}

/// Block (k, j) of Gamma is A^(k-j) B for j <= k, with k and j counted from zero.
inline auto lifted_input_response(const mpc_sweep_cell& cell, const mpc_sweep_plant& plant) -> Eigen::MatrixXd
{
    Eigen::MatrixXd gamma = Eigen::MatrixXd::Zero(cell.horizon * cell.nx, cell.horizon * cell.nu);
    Eigen::MatrixXd block = plant.actuation;
    for(std::int32_t d = 0; d < cell.horizon; ++d)
    {
        for(std::int32_t j = 0; j + d < cell.horizon; ++j)
            gamma.block((j + d) * cell.nx, j * cell.nu, cell.nx, cell.nu) = block;
        block = plant.transition * block;
    }
    return gamma;
}

inline auto stacked_state_weight(const mpc_tracking_spec& spec) -> Eigen::MatrixXd
{
    const std::int32_t states = spec.cell.horizon * spec.cell.nx;
    Eigen::MatrixXd weight = stage_state_weight * Eigen::MatrixXd::Identity(states, states);
    weight.bottomRightCorner(spec.cell.nx, spec.cell.nx) = spec.terminal_weight;
    return weight;
}

/// Selects the actuated states i < nu of stages 1..N, in the state rows' order.
inline auto actuated_rows(const mpc_sweep_cell& cell, const Eigen::MatrixXd& lifted) -> Eigen::MatrixXd
{
    Eigen::MatrixXd rows(cell.horizon * cell.nu, lifted.cols());
    for(std::int32_t k = 0; k < cell.horizon; ++k)
        rows.middleRows(k * cell.nu, cell.nu) = lifted.middleRows(k * cell.nx, cell.nu);
    return rows;
}

inline void condense_rows(const mpc_tracking_spec& spec, condensed_tracking& out)
{
    const std::int32_t inputs = spec.cell.horizon * spec.cell.nu;
    const std::int32_t state_rows = tracking_state_rows(spec);
    out.A = Eigen::MatrixXd::Zero(inputs + state_rows, inputs);
    out.A.topRows(inputs).setIdentity();
    out.l = Eigen::VectorXd::Constant(inputs + state_rows, -spec.variant.input_limit);
    out.u = Eigen::VectorXd::Constant(inputs + state_rows, spec.variant.input_limit);
    if(state_rows == 0)
        return;
    out.A.bottomRows(state_rows) = actuated_rows(spec.cell, out.gamma);
    out.actuated_phi = actuated_rows(spec.cell, out.phi);
    out.l.tail(state_rows).setConstant(-*spec.variant.state_limit);
    out.u.tail(state_rows).setConstant(*spec.variant.state_limit);
}

inline auto condense_tracking(const mpc_tracking_spec& spec) -> condensed_tracking
{
    const mpc_sweep_plant plant = sweep_plant(spec.cell);
    condensed_tracking out;
    out.phi = lifted_free_response(spec.cell, plant);
    out.gamma = lifted_input_response(spec.cell, plant);
    const Eigen::MatrixXd weight = stacked_state_weight(spec);
    out.gradient_map = out.gamma.transpose() * weight;
    out.stacked_target = tracking_target(spec).replicate(spec.cell.horizon, 1);
    const std::int32_t inputs = spec.cell.horizon * spec.cell.nu;
    out.P = out.gradient_map * out.gamma
          + stage_input_weight * Eigen::MatrixXd::Identity(inputs, inputs);
    condense_rows(spec, out);
    return out;
}

inline void condense_initial_state(const mpc_tracking_spec& spec, const condensed_tracking& posed,
                                   const Eigen::VectorXd& x0, Eigen::VectorXd& q, Eigen::VectorXd& l,
                                   Eigen::VectorXd& u)
{
    q.noalias() = posed.gradient_map * (posed.phi * x0 - posed.stacked_target);
    const std::int32_t state_rows = tracking_state_rows(spec);
    if(state_rows == 0)
        return;
    const Eigen::VectorXd shift = posed.actuated_phi * x0;
    l.tail(state_rows) = Eigen::VectorXd::Constant(state_rows, -*spec.variant.state_limit) - shift;
    u.tail(state_rows) = Eigen::VectorXd::Constant(state_rows, *spec.variant.state_limit) - shift;
}

}

#endif
