#ifndef HPP_GUARD_CTRLPP_MPC_QP_TYPES_H
#define HPP_GUARD_CTRLPP_MPC_QP_TYPES_H

#include "ctrlpp/types.h"

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

/// @brief Soft result status carried on the SUCCESS branch of a controller solve.
///
/// A controller `solve()` that returns a value always produces a usable control
/// input; this status refines how much to trust it. It is the caller-facing soft
/// vocabulary shared by both `mpc` and `nmpc`.
///
///  * converged          : the solver reached its convergence tolerances.
///  * solved_inaccurate  : a usable iterate was returned but tolerances were only
///                         partially met; treat with caution.
///  * budget_exhausted   : the iteration or time budget ran out before
///                         convergence, but the best iterate found is returned so
///                         the caller may still command it knowingly.
enum class solve_result_status : std::uint8_t
{
    converged,
    solved_inaccurate,
    budget_exhausted
};

/// @brief Hard failure modes carried on the ERROR branch of a controller solve.
///
/// These are the un-ignorable failures: no usable control input exists. The
/// caller cannot extract an input without first confronting the error, which is
/// exactly the property the `ctrlpp::expected` return channel provides.
///
///  * infeasible        : the problem as posed has no feasible point.
///  * invalid_problem   : the problem is unbounded, non-convex, or the solver
///                        reported an internal error; the data is not solvable.
///  * setup_incomplete  : the controller's one-time solver setup failed, so no
///                        solve can run at all.
enum class solver_error : std::uint8_t
{
    infeasible,
    invalid_problem,
    setup_incomplete
};

/// @brief Success payload of a controller `solve()`.
///
/// Aggregates the applied control `input` with the soft `status` describing how
/// the solve terminated. There is deliberately NO implicit conversion to
/// `Vector<Scalar, NU>`: call sites must reach the input explicitly through
/// `->input`, so the accompanying `status` can never be silently dropped.
template <typename Scalar, std::size_t NU>
struct solve_output
{
    Vector<Scalar, NU> input;
    solve_result_status status;
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
