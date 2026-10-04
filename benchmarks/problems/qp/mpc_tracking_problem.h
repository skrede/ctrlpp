#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_QP_MPC_TRACKING_PROBLEM_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_QP_MPC_TRACKING_PROBLEM_H

#include "qp/sparse_mpc_sweep.h"
#include "qp/mpc_tracking_cells.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>
#include <cstdint>

namespace ctrlpp::bench::problems::qp
{

/// The sweep's plant and decision layout generalized over the variant's input and state boxes, its
/// tracking target and the terminal weight. Rows come in the order initial state, dynamics, input
/// boxes, then (with a state box) the actuated states' boxes over stages 1..N, all stage-major.
struct mpc_tracking_spec
{
    mpc_sweep_cell cell;
    mpc_variant variant;
    Eigen::MatrixXd terminal_weight;
};

struct mpc_tracking_problem
{
    mpc_sweep_layout layout;
    std::int32_t state_rows;
    Eigen::SparseMatrix<double> P;
    Eigen::VectorXd q;
    Eigen::SparseMatrix<double> A;
    Eigen::VectorXd l;
    Eigen::VectorXd u;
};

inline auto tracking_target(const mpc_tracking_spec& spec) -> Eigen::VectorXd
{
    Eigen::VectorXd target = Eigen::VectorXd::Zero(spec.cell.nx);
    target.head(spec.cell.nu).setConstant(spec.variant.tracking_target);
    return target;
}

inline auto tracking_state_rows(const mpc_tracking_spec& spec) -> std::int32_t
{
    return spec.variant.state_limit ? spec.cell.horizon * spec.cell.nu : 0;
}

inline auto tracking_cost(const mpc_tracking_spec& spec, const mpc_sweep_layout& layout)
    -> Eigen::SparseMatrix<double>
{
    const mpc_sweep_cell& cell = spec.cell;
    const std::int32_t terminal = cell.horizon * cell.nx;
    std::vector<Eigen::Triplet<double>> entries;
    for(std::int32_t i = 0; i < terminal; ++i)
        entries.emplace_back(i, i, stage_state_weight);
    for(std::int32_t r = 0; r < cell.nx; ++r)
        for(std::int32_t c = 0; c < cell.nx; ++c)
            entries.emplace_back(terminal + r, terminal + c, spec.terminal_weight(r, c));
    for(std::int32_t i = layout.n_x; i < layout.n_dec; ++i)
        entries.emplace_back(i, i, stage_input_weight);
    Eigen::SparseMatrix<double> P(layout.n_dec, layout.n_dec);
    P.setFromTriplets(entries.begin(), entries.end());
    P.makeCompressed();
    return P;
}

inline auto tracking_gradient(const mpc_tracking_spec& spec, const mpc_sweep_layout& layout) -> Eigen::VectorXd
{
    const mpc_sweep_cell& cell = spec.cell;
    const Eigen::VectorXd target = tracking_target(spec);
    Eigen::VectorXd q = Eigen::VectorXd::Zero(layout.n_dec);
    for(std::int32_t k = 0; k < cell.horizon; ++k)
        q.segment(k * cell.nx, cell.nx) = -stage_state_weight * target;
    q.segment(cell.horizon * cell.nx, cell.nx) = -spec.terminal_weight * target;
    return q;
}

inline void append_box_rows(const mpc_tracking_spec& spec, const mpc_sweep_layout& layout,
                            std::vector<Eigen::Triplet<double>>& entries)
{
    const mpc_sweep_cell& cell = spec.cell;
    const std::int32_t inputs = cell.horizon * cell.nu;
    for(std::int32_t j = 0; j < inputs; ++j)
        entries.emplace_back(layout.n_dyn + j, layout.n_x + j, 1.0);
    for(std::int32_t k = 1; k <= cell.horizon && spec.variant.state_limit; ++k)
        for(std::int32_t i = 0; i < cell.nu; ++i)
            entries.emplace_back(layout.n_dyn + inputs + (k - 1) * cell.nu + i, k * cell.nx + i, 1.0);
}

inline void fill_box_bounds(const mpc_tracking_spec& spec, const mpc_sweep_layout& layout, Eigen::VectorXd& l,
                            Eigen::VectorXd& u)
{
    const std::int32_t inputs = spec.cell.horizon * spec.cell.nu;
    l.segment(layout.n_dyn, inputs).setConstant(-spec.variant.input_limit);
    u.segment(layout.n_dyn, inputs).setConstant(spec.variant.input_limit);
    if(!spec.variant.state_limit)
        return;
    l.tail(tracking_state_rows(spec)).setConstant(-*spec.variant.state_limit);
    u.tail(tracking_state_rows(spec)).setConstant(*spec.variant.state_limit);
}

/// Posed with a zero initial state; set_initial_state() writes each tick's.
inline auto build_tracking_problem(const mpc_tracking_spec& spec) -> mpc_tracking_problem
{
    const mpc_sweep_layout layout = sweep_layout(spec.cell);
    const std::int32_t state_rows = tracking_state_rows(spec);
    const std::int32_t rows = layout.n_con + state_rows;
    std::vector<Eigen::Triplet<double>> entries;
    append_dynamics_rows(spec.cell, layout, sweep_plant(spec.cell), entries);
    append_box_rows(spec, layout, entries);
    Eigen::SparseMatrix<double> A(rows, layout.n_dec);
    A.setFromTriplets(entries.begin(), entries.end());
    A.makeCompressed();
    Eigen::VectorXd l = Eigen::VectorXd::Zero(rows);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(rows);
    fill_box_bounds(spec, layout, l, u);
    return {layout, state_rows, tracking_cost(spec, layout), tracking_gradient(spec, layout), A, l, u};
}

inline void set_initial_state(const Eigen::VectorXd& x0, Eigen::VectorXd& l, Eigen::VectorXd& u)
{
    l.head(x0.size()) = x0;
    u.head(x0.size()) = x0;
}

}

#endif
