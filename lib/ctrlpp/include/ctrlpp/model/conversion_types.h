#ifndef HPP_GUARD_CTRLPP_MODEL_CONVERSION_TYPES_H
#define HPP_GUARD_CTRLPP_MODEL_CONVERSION_TYPES_H

namespace ctrlpp
{

enum class conversion_error
{
    non_finite_input,
    zero_leading_coefficient,
    eigenvalues_not_converged,
    unrepresentable_result,
};

}

#endif
