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

inline auto make_dense_mpc_lower_bound() -> Eigen::VectorXd
{
    return Eigen::VectorXd::Constant(dense_mpc_n_cons, -2.0);
}

inline auto make_dense_mpc_upper_bound() -> Eigen::VectorXd
{
    return Eigen::VectorXd::Constant(dense_mpc_n_cons, 2.0);
}

}

#endif
