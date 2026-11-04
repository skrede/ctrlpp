#ifndef HPP_GUARD_CTRLPP_MPC_QP_TYPES_H
#define HPP_GUARD_CTRLPP_MPC_QP_TYPES_H

#include "ctrlpp/types.h"

#include <Eigen/Sparse>

#include <cstdint>

namespace ctrlpp
{

// No backend reports invalid_backend_result: the consumer that finds a reported result short or non-finite sets it,
// so the diagnostics do not call a discarded answer a good solve.
enum class solve_status : std::uint8_t
{
    optimal,
    solved_inaccurate,
    infeasible,
    unbounded,
    max_iterations,
    time_limit,
    non_convex,
    error,
    invalid_backend_result
};

enum class qp_preset : std::uint8_t
{
    accuracy,
    speed
};

enum class solve_result_status : std::uint8_t
{
    converged,
    solved_inaccurate,
    budget_exhausted
};

enum class solver_error : std::uint8_t
{
    infeasible,
    invalid_problem,
    setup_incomplete,
    invalid_backend_result
};

enum class controller_construction_error : std::uint8_t
{
    non_positive_horizon,
    horizon_overflow,
    invalid_input_bounds,
    invalid_state_bounds,
    invalid_rate_bound,
    non_finite_cost_weight,
    asymmetric_cost_weight,
    indefinite_cost_weight,
    invalid_constraint_penalty,
    singular_input_weight,
    non_finite_system,
    non_finite_state_space_weight,
    mismatched_terminal_set,
    non_finite_terminal_set,
    asymmetric_terminal_set,
    indefinite_terminal_set,
    negative_terminal_set_level
};

// No conversion to Vector<Scalar, NU>: a caller reaches the input through ->input, so the status is never dropped
// unseen.
template <typename Scalar, std::size_t NU>
struct solve_output
{
    Vector<Scalar, NU> input;
    solve_result_status status;
};

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
    qp_result()
        : status{solve_status::error}
        , x{}
        , y{}
        , objective{}
        , solve_time{}
        , iterations{}
        , primal_residual{}
        , dual_residual{}
    {
    }

    solve_status status;
    Eigen::VectorX<Scalar> x;
    Eigen::VectorX<Scalar> y;
    Scalar objective;
    Scalar solve_time;
    std::int32_t iterations;
    Scalar primal_residual;
    Scalar dual_residual;
};

}

#endif
