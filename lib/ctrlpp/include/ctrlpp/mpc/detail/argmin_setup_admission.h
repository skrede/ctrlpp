#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_ARGMIN_SETUP_ADMISSION_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_ARGMIN_SETUP_ADMISSION_H

#include "ctrlpp/mpc/nlp_types.h"

#include <Eigen/Core>

#include <cmath>
#include <cstdint>
#include <optional>

namespace ctrlpp::detail
{

template <typename Scalar>
bool is_equality_row(Scalar lower, Scalar upper)
{
    return lower == upper;
}

template <typename Scalar>
bool is_finite_side(Scalar bound)
{
    return std::isfinite(static_cast<double>(bound));
}

// 64-bit, so counting two sides per row cannot overflow where Eigen::Index is 32-bit.
struct constraint_row_count
{
    std::int64_t equality;
    std::int64_t inequality;
};

template <typename Problem>
constraint_row_count count_constraint_rows(const Problem& problem)
{
    constraint_row_count count{0, 0};
    for(Eigen::Index i = 0; i < problem.n_constraints; ++i)
    {
        if(is_equality_row(problem.c_lower[i], problem.c_upper[i]))
            ++count.equality;
        else
            count.inequality += std::int64_t{is_finite_side(problem.c_upper[i])} + std::int64_t{is_finite_side(problem.c_lower[i])};
    }
    return count;
}

// The order is the bridge's kind, then the fixed shapes, then the row counts.
template <bool Constrained, int NV, int MaxM, typename Problem>
std::optional<argmin_setup_error> argmin_setup_refusal(const Problem& problem)
{
    if(!Constrained && problem.n_constraints > 0)
        return argmin_setup_error::constraints_unsupported;
    if(NV != Eigen::Dynamic && problem.n_vars != NV)
        return argmin_setup_error::dimension_mismatch;
    if constexpr(Constrained)
    {
        const constraint_row_count rows = count_constraint_rows(problem);
        if(MaxM != Eigen::Dynamic && rows.equality + rows.inequality > MaxM)
            return argmin_setup_error::constraint_cap_exceeded;
        if(rows.equality > problem.n_vars)
            return argmin_setup_error::equality_rows_exceed_dimension;
    }
    return std::nullopt;
}

}

#endif
