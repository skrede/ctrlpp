#ifndef HPP_GUARD_CTRLPP_MPC_QP_TYPES_H
#define HPP_GUARD_CTRLPP_MPC_QP_TYPES_H

#include <Eigen/Sparse>

#include <cstdint>

namespace ctrlpp
{

enum class solve_status : std::uint8_t
{
    optimal,
    solved_inaccurate,
    infeasible,
    unbounded,
    max_iterations,
    time_limit,
    non_convex,
    error
};

/// @brief Structured failure modes for `osqp_solver::try_setup`.
///
///  * settings_allocation_failed : `OSQPSettings_new` returned null; the settings
///                                 block could not be allocated.
///  * setup_failed               : `osqp_setup` returned a nonzero exit flag; the
///                                 problem data or settings were rejected, or the
///                                 workspace could not be allocated.
enum class osqp_setup_error : std::uint8_t
{
    settings_allocation_failed,
    setup_failed
};

template <typename Scalar>
struct qp_problem
{
    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> P;
    Eigen::VectorX<Scalar> q;
    Eigen::SparseMatrix<Scalar, Eigen::ColMajor> A;
    Eigen::VectorX<Scalar> l;
    Eigen::VectorX<Scalar> u;
};

template <typename Scalar>
struct qp_update
{
    Eigen::VectorX<Scalar> q;
    Eigen::VectorX<Scalar> l;
    Eigen::VectorX<Scalar> u;
    Eigen::VectorX<Scalar> warm_x;
    Eigen::VectorX<Scalar> warm_y;
};

template <typename Scalar>
struct qp_result
{
    solve_status status{solve_status::error};
    Eigen::VectorX<Scalar> x{};
    Eigen::VectorX<Scalar> y{};
    Scalar objective{};
    Scalar solve_time{};
    int iterations{};
    Scalar primal_residual{};
    Scalar dual_residual{};
};

}

#endif
