#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_QP_REFUSAL_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_ARGMIN_QP_REFUSAL_H

#include <argmin/qp/qp_types.h>

#include <cstdio>
#include <optional>

namespace ctrlpp::bench
{

inline char const* qp_error_name(argmin::qp_error error)
{
    switch(error)
    {
        case argmin::qp_error::dimension_mismatch:
            return "dimension_mismatch";
        case argmin::qp_error::invalid_bounds:
            return "invalid_bounds";
        case argmin::qp_error::non_finite_input:
            return "non_finite_input";
        case argmin::qp_error::capacity_exceeded:
            return "capacity_exceeded";
        case argmin::qp_error::infeasible_start:
            return "infeasible_start";
        case argmin::qp_error::invalid_problem:
            return "invalid_problem";
        case argmin::qp_error::invalid_options:
            return "invalid_options";
        case argmin::qp_error::factorization_failed:
            return "factorization_failed";
    }
    return "an unnamed qp_error";
}

inline bool admm_refused(char const* call, const std::optional<argmin::qp_error>& error)
{
    if(!error)
        return false;
    std::fprintf(stderr, "argmin::sparse_admm_qp_solver::%s refused with %s; a timing measured against a solver that did not solve is meaningless\n", call,
                 qp_error_name(*error));
    return true;
}

}

#endif
