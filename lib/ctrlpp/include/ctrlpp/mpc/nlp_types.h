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

/// @brief Structured failure modes for `argmin_solver::try_setup`.
///
///  * incompatible_equality_constraints : the selected policy is a raw
///    MMA-family algorithm (`argmin_mma` / `argmin_gcmma`) which cannot
///    represent equality constraints. The auglag-wrapped variants
///    (`argmin_auglag<argmin_mma>`, `argmin_auglag<argmin_gcmma>`) absorb
///    equality constraints into the outer augmented-Lagrangian penalty and
///    accept the same problem.
enum class argmin_setup_error : std::uint8_t
{
    incompatible_equality_constraints
};

}

#endif
