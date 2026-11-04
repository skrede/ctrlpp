#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_ARGMIN_OPTIONS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_ARGMIN_OPTIONS_H

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <argmin/result/status.h>

#include <argmin/solver/options.h>
#include <argmin/solver/time_budget_options.h>

#include <chrono>
#include <cstdint>

namespace ctrlpp::detail
{

template <bool Constrained, typename Scalar>
inline argmin::solver_options<argmin_ctrlpp_convergence> make_argmin_solver_options(const argmin_settings<Scalar>& s)
{
    argmin::solver_options<argmin_ctrlpp_convergence> opts;

    opts.max_iterations = static_cast<std::uint32_t>(s.max_eval);

    if constexpr(Constrained)
    {
        if(s.constraint_tol > Scalar{0})
        {
            opts.constraint_tolerance = static_cast<double>(s.constraint_tol);
            opts.feasibility_tolerance = static_cast<double>(s.constraint_tol);
        }
    }

    opts.set_objective_threshold_rel(static_cast<double>(s.ftol_rel));
    opts.set_step_threshold_rel(static_cast<double>(s.xtol_rel));
    opts.set_gradient_threshold(static_cast<double>(s.kkt_tol));

    return opts;
}

template <bool Constrained, typename Scalar>
inline argmin::time_budget_options<argmin_ctrlpp_convergence> make_argmin_time_budget_options(const argmin_timed_settings<Scalar>& s)
{
    argmin::time_budget_options<argmin_ctrlpp_convergence> opts;
    opts.core = make_argmin_solver_options<Constrained>(s);

    opts.max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::duration<double>(static_cast<double>(s.max_time)));

    return opts;
}

constexpr solve_status map_argmin_status(argmin::solver_status s)
{
    switch(s)
    {
        case argmin::solver_status::converged:
        case argmin::solver_status::ftol_reached:
        case argmin::solver_status::xtol_reached:
            return solve_status::optimal;
        case argmin::solver_status::max_iterations:
        case argmin::solver_status::budget_exhausted:
        case argmin::solver_status::maxeval_reached:
            return solve_status::max_iterations;
        case argmin::solver_status::time_limit_reached:
            return solve_status::time_limit;
        case argmin::solver_status::stalled:
        case argmin::solver_status::roundoff_limited:
        case argmin::solver_status::objective_stalled:
        case argmin::solver_status::line_search_failed:
        case argmin::solver_status::trust_region_step_rejected:
            return solve_status::solved_inaccurate;
        case argmin::solver_status::invalid_problem:
        case argmin::solver_status::diverged:
        case argmin::solver_status::aborted:
        case argmin::solver_status::running:
            return solve_status::error;
    }
    return solve_status::error;
}

constexpr nlp_stop_criterion map_argmin_stop_criterion(argmin::solver_status s)
{
    switch(s)
    {
        case argmin::solver_status::converged:
            return nlp_stop_criterion::stationarity;
        case argmin::solver_status::ftol_reached:
            return nlp_stop_criterion::objective_change;
        case argmin::solver_status::xtol_reached:
            return nlp_stop_criterion::step_size;
        case argmin::solver_status::stalled:
        case argmin::solver_status::objective_stalled:
        case argmin::solver_status::line_search_failed:
            return nlp_stop_criterion::stalled;
        case argmin::solver_status::max_iterations:
        case argmin::solver_status::budget_exhausted:
        case argmin::solver_status::maxeval_reached:
            return nlp_stop_criterion::iteration_limit;
        case argmin::solver_status::time_limit_reached:
            return nlp_stop_criterion::time_limit;
        case argmin::solver_status::roundoff_limited:
        case argmin::solver_status::trust_region_step_rejected:
        case argmin::solver_status::invalid_problem:
        case argmin::solver_status::diverged:
        case argmin::solver_status::aborted:
        case argmin::solver_status::running:
            return nlp_stop_criterion::other;
    }
    return nlp_stop_criterion::other;
}

}

#endif
