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
enum class nlopt_setup_error : std::uint8_t
{
    incompatible_equality_constraints
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
///  * stalled          : the iterate or the objective stopped making progress.
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

/// @brief Structured failure modes for `argmin_solver::setup`. The argmin
/// adapter has no setup failure mode, so this carries no enumerators; the
/// fallible signature is what the solver concept requires of every backend.
enum class argmin_setup_error : std::uint8_t
{
};

}

#endif
