#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_QP_SPARSE_MPC_SWEEP_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_QP_SPARSE_MPC_SWEEP_H

#include "qp/qp_accuracy.h"

#include "ctrlpp/mpc/qp_types.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

namespace ctrlpp::bench::problems::qp
{

struct mpc_sweep_cell
{
    int nx;
    int nu;
    int horizon;
};

struct mpc_sweep_layout
{
    int n_x;
    int n_con;
    int n_dec;
    int n_dyn;
};

struct mpc_sweep_plant
{
    Eigen::MatrixXd transition;
    Eigen::MatrixXd actuation;
};

struct mpc_sweep_bounds
{
    Eigen::VectorXd l;
    Eigen::VectorXd u;
    Eigen::SparseMatrix<double> A;
};

struct mpc_sweep_problem
{
    Eigen::VectorXd q;
    Eigen::VectorXd l;
    Eigen::VectorXd u;
    ctrlpp::qp_problem<double> problem;
};

/// Tiny-embedded through mid-scale, matching the cells the nonlinear
/// predictive-control sweeps in the same family already use.
inline auto mpc_sweep_cells() -> std::vector<mpc_sweep_cell>
{
    return {{2, 1, 10}, {4, 2, 10}, {4, 2, 20}, {8, 3, 20}, {8, 3, 30}, {12, 4, 30}};
}

inline auto sweep_layout(const mpc_sweep_cell& cell) -> mpc_sweep_layout
{
    const int n_x = (cell.horizon + 1) * cell.nx;
    return {n_x, n_x + cell.horizon * cell.nu, n_x + cell.horizon * cell.nu, n_x};
}

/// A mildly-coupled stable plant: spectral radius below one, the inputs acting
/// on the leading states.
inline auto sweep_plant(const mpc_sweep_cell& cell) -> mpc_sweep_plant
{
    Eigen::MatrixXd transition = 0.9 * Eigen::MatrixXd::Identity(cell.nx, cell.nx);
    for(int i = 0; i + 1 < cell.nx; ++i)
        transition(i, i + 1) = 0.1;
    Eigen::MatrixXd actuation = Eigen::MatrixXd::Zero(cell.nx, cell.nu);
    for(int i = 0; i < cell.nu; ++i)
        actuation(i, i) = 1.0;
    return {transition, actuation};
}

inline auto sweep_cost(const mpc_sweep_cell& cell, const mpc_sweep_layout& layout) -> Eigen::SparseMatrix<double>
{
    constexpr double state_weight = 1.0;
    constexpr double input_weight = 0.1;
    constexpr double terminal_weight = 10.0;
    std::vector<Eigen::Triplet<double>> entries;
    for(int k = 0; k <= cell.horizon; ++k)
        for(int i = 0; i < cell.nx; ++i)
            entries.emplace_back(k * cell.nx + i, k * cell.nx + i,
                                 k == cell.horizon ? terminal_weight : state_weight);
    for(int k = 0; k < cell.horizon; ++k)
        for(int i = 0; i < cell.nu; ++i)
            entries.emplace_back(layout.n_x + k * cell.nu + i, layout.n_x + k * cell.nu + i, input_weight);
    Eigen::SparseMatrix<double> P(layout.n_dec, layout.n_dec);
    P.setFromTriplets(entries.begin(), entries.end());
    P.makeCompressed();
    return P;
}

/// The initial condition as an equality on the first state block, then one
/// transition equality per step.
inline void append_dynamics_rows(const mpc_sweep_cell& cell, const mpc_sweep_layout& layout,
                                 const mpc_sweep_plant& plant, std::vector<Eigen::Triplet<double>>& entries)
{
    for(int i = 0; i < cell.nx; ++i)
        entries.emplace_back(i, i, 1.0);
    for(int k = 0; k < cell.horizon; ++k)
        for(int i = 0; i < cell.nx; ++i)
        {
            const int row = cell.nx + k * cell.nx + i;
            entries.emplace_back(row, (k + 1) * cell.nx + i, 1.0);
            for(int j = 0; j < cell.nx; ++j)
                entries.emplace_back(row, k * cell.nx + j, -plant.transition(i, j));
            for(int j = 0; j < cell.nu; ++j)
                entries.emplace_back(row, layout.n_x + k * cell.nu + j, -plant.actuation(i, j));
        }
}

inline auto sweep_bounds(const mpc_sweep_cell& cell, const mpc_sweep_layout& layout, const mpc_sweep_plant& plant)
    -> mpc_sweep_bounds
{
    constexpr double initial_state = 1.0;
    constexpr double input_limit = 10.0;
    std::vector<Eigen::Triplet<double>> entries;
    append_dynamics_rows(cell, layout, plant, entries);
    Eigen::VectorXd l = Eigen::VectorXd::Zero(layout.n_con);
    Eigen::VectorXd u = Eigen::VectorXd::Zero(layout.n_con);
    l.head(cell.nx).setConstant(initial_state);
    u.head(cell.nx).setConstant(initial_state);
    for(int k = 0; k < cell.horizon; ++k)
        for(int i = 0; i < cell.nu; ++i)
        {
            const int row = layout.n_dyn + k * cell.nu + i;
            entries.emplace_back(row, layout.n_x + k * cell.nu + i, 1.0);
            l(row) = -input_limit;
            u(row) = input_limit;
        }
    Eigen::SparseMatrix<double> A(layout.n_con, layout.n_dec);
    A.setFromTriplets(entries.begin(), entries.end());
    A.makeCompressed();
    return {l, u, A};
}

/// Regulation of the plant to the origin from a nonzero initial state, so every
/// cell is feasible and both solvers converge to the same optimum.
inline auto build_mpc_sweep_problem(const mpc_sweep_cell& cell) -> mpc_sweep_problem
{
    const mpc_sweep_layout layout = sweep_layout(cell);
    const mpc_sweep_bounds bounds = sweep_bounds(cell, layout, sweep_plant(cell));
    const Eigen::VectorXd q = Eigen::VectorXd::Zero(layout.n_dec);
    return {q, bounds.l, bounds.u,
            {.P = sweep_cost(cell, layout), .q = q, .A = bounds.A, .l = bounds.l, .u = bounds.u}};
}

inline auto dense_form(const mpc_sweep_problem& posed) -> dense_program
{
    return {Eigen::MatrixXd(posed.problem.P), posed.q, Eigen::MatrixXd(posed.problem.A), posed.l, posed.u};
}

}

#endif
