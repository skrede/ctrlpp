#ifndef HPP_GUARD_CTRLPP_MPC_QP_TYPES_H
#define HPP_GUARD_CTRLPP_MPC_QP_TYPES_H

#include "ctrlpp/types.h"

#include <Eigen/Sparse>

#include <cstdint>

namespace ctrlpp
{

/// @brief Outcome of a single solve, as reported to the caller.
///
/// Every value except the last is a backend's own termination status, mapped
/// into this shared vocabulary by the backend bridge.
///
///  * invalid_backend_result : the backend reported a termination status of its
///                             own but returned a result whose dimensions do not
///                             match the posed problem. No backend ever reports
///                             this; it is set by the consumer that detected the
///                             mismatch, so a reader of the diagnostics is not
///                             told the solve went well when its answer was
///                             discarded. It is the diagnostics-channel twin of
///                             solver_error::invalid_backend_result.
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

/// @brief Backend-agnostic QP tuning preset, selecting the accuracy/speed
/// tradeoff for any solver modeling the `qp_solver` concept.
///
/// The presets differ in a single knob: solution polishing. Polishing runs an
/// extra active-set refinement (re-solving a reduced KKT system) that drives the
/// operator-splitting iterate from its ~tolerance-level answer to near machine
/// precision, at a per-solve cost that can approach that of the ADMM loop itself.
///
///  * accuracy : polish on. The refined, high-precision iterate. Use when the
///               QP solution feeds something tolerance-sensitive (e.g. an SQP
///               inner solve) or when tight constraint feasibility is required.
///  * speed    : polish off. The raw operator-splitting iterate, converged to the
///               solver's stopping tolerance. For warm-resolve linear MPC the
///               per-step polish refinement is unnecessary -- the unpolished
///               iterate already meets the control tolerance -- so this is the
///               better default there.
enum class qp_preset : std::uint8_t
{
    accuracy,
    speed
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
///  * infeasible             : the problem as posed has no feasible point.
///  * invalid_problem        : the problem is unbounded, non-convex, or the
///                             solver reported an internal error; the data is
///                             not solvable.
///  * setup_incomplete       : the controller's one-time solver setup failed, so
///                             no solve can run at all.
///  * invalid_backend_result : the solver reported a status the controller
///                             accepts and then returned a primal or dual whose
///                             length does not cover the dimensions of the posed
///                             problem, so the fixed-width slices the extraction
///                             takes out of it would read past its end. This is
///                             deliberately NOT folded into invalid_problem:
///                             there the caller must fix the problem it posed,
///                             here the problem is well formed and the backend's
///                             answer is not, so the two demand different
///                             remedies and must stay distinguishable.
enum class solver_error : std::uint8_t
{
    infeasible,
    invalid_problem,
    setup_incomplete,
    invalid_backend_result
};

/// @brief Structured failure modes for the controller construction factories.
///
/// A runtime prediction horizon is a caller-supplied signed value that scales
/// every derived dimension of the posed problem, so it is validated once, at
/// construction, before any dimension product is formed or any storage is
/// reserved. The horizon is deliberately kept signed: a mistaken negative value
/// stays representable as negative and is therefore rejectable, whereas an
/// unsigned field would silently turn the same mistake into an enormous
/// allocation.
///
///  * non_positive_horizon : the prediction horizon is zero or negative. A
///                           horizon of zero poses no input to optimize over and
///                           leaves the first-input extraction reading past the
///                           end of the decision vector; a negative horizon
///                           drives every derived dimension negative.
///  * horizon_overflow     : the horizon is large enough that the derived
///                           decision or constraint dimension would not be
///                           representable in the horizon's own type, so the
///                           products that size the problem would wrap.
enum class controller_construction_error : std::uint8_t
{
    non_positive_horizon,
    horizon_overflow
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
