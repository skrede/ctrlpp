#ifndef HPP_GUARD_CTRLPP_MPC_NLP_TYPES_H
#define HPP_GUARD_CTRLPP_MPC_NLP_TYPES_H

/// @brief Setup-error types for the NLP solver backends.

#include <cstdint>

namespace ctrlpp
{

/// @brief Structured failure modes for `nlopt_solver::setup`.
///
///  * incompatible_equality_constraints : the selected algorithm (raw MMA or raw
///    CCSAQ) cannot handle equality constraints. The auglag-wrapped variants
///    (`auglag_mma`, `auglag_ccsaq`) absorb equality constraints into the outer
///    augmented-Lagrangian penalty and accept the same problem.
///  * backend_rejected : NLopt or the adapter raised an error other than an
///    allocation failure while building or configuring the optimizer (an
///    argument NLopt rejects, such as a bound vector of the wrong length, or any
///    other error, standard or not), so no optimizer exists.
///  * allocation_failed : an allocation failed while the optimizer was built or
///    configured (`std::bad_alloc`, which NLopt also raises for its own
///    out-of-memory), so no optimizer exists.
enum class nlopt_setup_error : std::uint8_t
{
    incompatible_equality_constraints,
    backend_rejected,
    allocation_failed
};

/// @brief Structured failure modes for the compile-time-horizon NLP
/// formulation factory `build_nmpc_problem_static`.
///
/// The compile-time decision dimension NV = (NH+1)*NX + NH*NU pins the shape of
/// the posed problem, so a configuration that would produce a different runtime
/// dimension is rejected before any of the dependent problem data is built.
///
///  * horizon_mismatch    : `config.horizon` disagrees with the compile-time
///                          horizon NH, so every derived offset into the
///                          decision vector would be computed against the wrong
///                          dimension.
///  * slack_not_supported : a soft-constraint configuration would add slack
///                          decision variables, growing the decision dimension
///                          past NV. Only the hard-constraint (slack-free) cut
///                          is supported on the compile-time-dimension path.
enum class nlp_formulation_error : std::uint8_t
{
    horizon_mismatch,
    slack_not_supported
};

/// @brief Which stopping test ended an NLP solve, for backends that report it.
///
/// `nlp_result::status` folds every successful stop into `solve_status::optimal`,
/// which is the right summary for a caller acting on the answer but cannot say
/// what the answer is certified to. A stop on `stationarity` means the backend's
/// first-order optimality residual fell below its threshold, which is what a
/// bound on the distance to the exact optimum can be built from; a stop on the
/// change in the objective or in the iterate certifies no such distance.
///
///  * unreported       : the backend does not report its stopping test.
///  * stationarity     : the first-order optimality (KKT) residual fell below its
///                       threshold (`argmin_settings::kkt_tol`).
///  * objective_change : the relative change in the objective fell below its
///                       threshold.
///  * step_size        : the relative change in the iterate fell below its
///                       threshold.
///  * stalled          : the iterate or the objective stopped making progress,
///                       including a line search that accepted no step.
///  * iteration_limit  : the iteration or evaluation budget ran out.
///  * time_limit       : the wall-clock budget ran out.
///  * other            : any other stop, including a failed solve.
enum class nlp_stop_criterion : std::uint8_t
{
    unreported,
    stationarity,
    objective_change,
    step_size,
    stalled,
    iteration_limit,
    time_limit,
    other
};

/// @brief Structured failure modes for `argmin_solver::setup`, listed in the
/// order setup checks them, so a problem that fails more than one check is
/// refused with the first.
///
///  * constraints_unsupported : the problem poses general constraint rows
///    (`n_constraints > 0`) on the unconstrained bridge (Constrained = false),
///    which handles variable bounds only, so it is refused rather than solved
///    without them.
///  * dimension_mismatch : the fixed-dimension bridge is sized for the
///    compile-time decision dimension NV, and the problem declares another
///    `n_vars`.
///  * constraint_cap_exceeded : the problem's partitioned constraint row count
///    (equality rows plus one per finite inequality side, so a two-sided row
///    counts twice) exceeds the compile-time cap MaxM, so argmin's inline
///    result-multiplier storage could not hold the multipliers.
///  * equality_rows_exceed_dimension : the problem poses more equality rows than
///    decision variables. argmin's active-set QP needs at most n equality rows.
enum class argmin_setup_error : std::uint8_t
{
    constraints_unsupported,
    constraint_cap_exceeded,
    dimension_mismatch,
    equality_rows_exceed_dimension
};

}

#endif
