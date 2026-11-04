#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_ARGMIN_PROBLEM_ROWS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_ARGMIN_PROBLEM_ROWS_H

#include "ctrlpp/mpc/detail/argmin_setup_admission.h"

#include <Eigen/Core>

#include <vector>

namespace ctrlpp::detail
{

// argmin reads a compile-time constraint cap from the presence of constraint_count (has_constraint_count in
// argmin/formulation/concepts.h) and asserts runtime_m <= cap for any value present, so the dynamic case omits the
// member rather than carrying a sentinel.
template <int MaxM>
struct argmin_constraint_bound
{
    static constexpr int constraint_count = MaxM;
};

template <>
struct argmin_constraint_bound<Eigen::Dynamic>
{
};

// The same classification as count_constraint_rows, which setup admits against MaxM.
template <typename Problem>
inline void classify_constraint_rows(
        const Problem& prob, std::vector<Eigen::Index>& eq_indices, std::vector<Eigen::Index>& ineq_upper_indices, std::vector<Eigen::Index>& ineq_lower_indices)
{
    eq_indices.clear();
    ineq_upper_indices.clear();
    ineq_lower_indices.clear();

    for(Eigen::Index i = 0; i < prob.n_constraints; ++i)
    {
        if(is_equality_row(prob.c_lower[i], prob.c_upper[i]))
        {
            eq_indices.push_back(i);
            continue;
        }
        if(is_finite_side(prob.c_upper[i]))
            ineq_upper_indices.push_back(i);
        if(is_finite_side(prob.c_lower[i]))
            ineq_lower_indices.push_back(i);
    }
}

}

#endif
