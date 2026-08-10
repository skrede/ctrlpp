#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_QP_DENSE_MPC_SHAPED_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_QP_DENSE_MPC_SHAPED_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>

namespace ctrlpp::bench::problems::qp
{

inline constexpr int dense_mpc_n_vars = 10;
inline constexpr int dense_mpc_n_cons = 5;

inline auto make_dense_mpc_hessian() -> Eigen::SparseMatrix<double>
{
    Eigen::SparseMatrix<double> H(dense_mpc_n_vars, dense_mpc_n_vars);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < dense_mpc_n_vars; ++i)
    {
        triplets.emplace_back(i, i, 2.0 + 0.1 * i);
        if(i + 1 < dense_mpc_n_vars)
        {
            triplets.emplace_back(i, i + 1, 0.1);
            triplets.emplace_back(i + 1, i, 0.1);
        }
    }
    H.setFromTriplets(triplets.begin(), triplets.end());
    return H;
}

inline auto make_dense_mpc_constraint_matrix() -> Eigen::SparseMatrix<double>
{
    Eigen::SparseMatrix<double> A(dense_mpc_n_cons, dense_mpc_n_vars);
    std::vector<Eigen::Triplet<double>> triplets;
    for(int i = 0; i < dense_mpc_n_cons; ++i)
    {
        triplets.emplace_back(i, 2 * i, 1.0);
        triplets.emplace_back(i, 2 * i + 1, 0.5);
    }
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

inline auto make_dense_mpc_gradient() -> Eigen::VectorXd
{
    return Eigen::VectorXd::LinSpaced(dense_mpc_n_vars, -1.0, 1.0);
}

/// How much each constraint row is asked for by the minimizer of the objective
/// alone. It is the only scale the posed program itself offers for how large a
/// bound has to be to matter, so the bounds below are derived from it rather
/// than chosen.
inline auto dense_mpc_row_demand() -> Eigen::VectorXd
{
    const Eigen::MatrixXd P(make_dense_mpc_hessian());
    const Eigen::MatrixXd A(make_dense_mpc_constraint_matrix());
    return (A * P.ldlt().solve(-make_dense_mpc_gradient())).cwiseAbs();
}

/// Bounds sit at the mean demand, so the rows that ask for more than average
/// bind and the rest stay slack. A corpus whose rows are all slack races the
/// solvers on an effectively unconstrained dense program; one whose rows are all
/// tight poses an equality-constrained problem with no active set to find.
inline auto dense_mpc_bound_magnitude() -> double
{
    return dense_mpc_row_demand().mean();
}

inline auto make_dense_mpc_lower_bound(double magnitude = dense_mpc_bound_magnitude()) -> Eigen::VectorXd
{
    return Eigen::VectorXd::Constant(dense_mpc_n_cons, -magnitude);
}

inline auto make_dense_mpc_upper_bound(double magnitude = dense_mpc_bound_magnitude()) -> Eigen::VectorXd
{
    return Eigen::VectorXd::Constant(dense_mpc_n_cons, magnitude);
}

}

#endif
