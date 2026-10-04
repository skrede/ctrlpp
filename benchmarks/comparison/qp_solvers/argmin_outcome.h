#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ARGMIN_OUTCOME_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ARGMIN_OUTCOME_H

#include "arm.h"
#include "registry.h"

#include <argmin/qp/qp_types.h>

#include <array>
#include <optional>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

inline auto status_token(argmin::qp_solve_status status) -> std::string_view
{
    constexpr std::array<std::string_view, 5> names{"solved", "solved_inaccurate", "max_iterations", "primal_infeasible",
                                                    "dual_infeasible"};
    return names[static_cast<std::size_t>(status)];
}

inline auto error_token(argmin::qp_error error) -> std::string_view
{
    constexpr std::array<std::string_view, 8> names{"qp_error_dimension_mismatch", "qp_error_invalid_bounds",
                                                    "qp_error_non_finite_input", "qp_error_capacity_exceeded",
                                                    "qp_error_infeasible_start", "qp_error_invalid_problem",
                                                    "qp_error_invalid_options", "qp_error_factorization_failed"};
    return names[static_cast<std::size_t>(error)];
}

inline auto claim_of(argmin::qp_solve_status status) -> solver_claim
{
    switch(status)
    {
        case argmin::qp_solve_status::solved:
        case argmin::qp_solve_status::solved_inaccurate:
            return solver_claim::solved;
        case argmin::qp_solve_status::max_iterations:
            return solver_claim::iteration_cap;
        default:
            return solver_claim::infeasible;
    }
}

template <typename Result>
void fill_answer(const Result& result, const std::optional<argmin::qp_error>& error, arm_answer& out)
{
    out.x = error ? Eigen::VectorXd() : Eigen::VectorXd(result.x);
    out.y = error ? Eigen::VectorXd() : Eigen::VectorXd(result.y);
    out.iterations = error ? -1 : result.iterations;
    out.status = error ? error_token(*error) : status_token(result.status);
    out.claim = error ? solver_claim::error : claim_of(result.status);
}

template <typename Options>
auto admm_options(tier level, bool forced_polish) -> Options
{
    Options options;
    options.eps_abs = solver_tolerance(level);
    options.eps_rel = solver_tolerance(level);
    options.max_iterations = iteration_cap;
    options.warm_start = true;
    options.polish = level == tier::high || forced_polish;
    return options;
}

}

#endif
