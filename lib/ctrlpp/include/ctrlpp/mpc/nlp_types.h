#ifndef HPP_GUARD_CTRLPP_MPC_NLP_TYPES_H
#define HPP_GUARD_CTRLPP_MPC_NLP_TYPES_H

/// @brief Setup-error types for the NLP solver backends.

#include <cstdint>

namespace ctrlpp
{

/// @brief Structured failure modes for `nlopt_solver::try_setup`.
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

/// @brief Structured failure modes for `argmin_solver::try_setup`. The argmin
/// adapter has no setup failure mode, so this carries no enumerators; the
/// fallible signature is retained for parity with the other solver backends.
enum class argmin_setup_error : std::uint8_t
{
};

}

#endif
