#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_INSTANCE_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_INSTANCE_H

#include "qp/qp_accuracy.h"
#include "qp/mpc_tracking_cells.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cstdint>
#include <optional>

namespace ctrlpp::bench::qp_solvers
{

struct ocp_structure
{
    problems::qp::mpc_sweep_cell cell;
    Eigen::MatrixXd transition;
    Eigen::MatrixXd actuation;
    Eigen::MatrixXd terminal_weight;
    Eigen::VectorXd target;
    double input_limit;
    std::optional<double> state_limit;
};

/// The data that stays fixed over one sequence of ticks. A row is an equality when its bounds coincide
/// on the posing tick, and stays one on every later tick. The leading box_rows rows are identity rows
/// on the leading variables, which the solvers that separate variable bounds take as such.
struct qp_instance
{
    Eigen::MatrixXd P;
    Eigen::MatrixXd A;
    Eigen::SparseMatrix<double> P_sparse;
    Eigen::SparseMatrix<double> A_sparse;
    std::int32_t box_rows;
    std::optional<ocp_structure> ocp;
};

struct qp_tick
{
    Eigen::VectorXd q;
    Eigen::VectorXd l;
    Eigen::VectorXd u;
    Eigen::VectorXd x0;
};

inline auto program_of(const qp_instance& instance, const qp_tick& tick) -> problems::qp::dense_program
{
    return {instance.P, tick.q, instance.A, tick.l, tick.u};
}

}

#endif
