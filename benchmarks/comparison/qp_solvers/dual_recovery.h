#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_DUAL_RECOVERY_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_DUAL_RECOVERY_H

#include "qp/qp_accuracy.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace ctrlpp::bench::qp_solvers
{

/// Multipliers for a primal answer that came without any: the least-squares solution of
/// P x + q + A_act^T y_act = 0 over the equality rows and the rows within threshold (1 + |bound|) of a
/// bound, zero on every other row. A sign that contradicts the active bound is left in place for the
/// judge's complementarity term to catch.
auto recover_duals(const problems::qp::dense_program& program, const Eigen::SparseMatrix<double>& A,
                   const Eigen::VectorXd& x, double threshold) -> Eigen::VectorXd;

}

#endif
