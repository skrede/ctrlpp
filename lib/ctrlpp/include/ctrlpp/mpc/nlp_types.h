#ifndef HPP_GUARD_CTRLPP_MPC_NLP_TYPES_H
#define HPP_GUARD_CTRLPP_MPC_NLP_TYPES_H

#include <cstdint>

namespace ctrlpp
{

enum class nlopt_setup_error : std::uint8_t
{
    incompatible_equality_constraints,
    backend_rejected,
    allocation_failed
};

enum class nlp_formulation_error : std::uint8_t
{
    horizon_mismatch,
    slack_not_supported
};

// Only a stop on stationarity bounds the distance to the exact optimum; a stop on the change in the objective
// or in the iterate certifies no such distance.
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

// The enumerators are not in the order setup checks them: setup refuses with constraints_unsupported, then
// dimension_mismatch, then constraint_cap_exceeded, then equality_rows_exceed_dimension.
enum class argmin_setup_error : std::uint8_t
{
    constraints_unsupported,
    constraint_cap_exceeded,
    dimension_mismatch,
    equality_rows_exceed_dimension
};

}

#endif
