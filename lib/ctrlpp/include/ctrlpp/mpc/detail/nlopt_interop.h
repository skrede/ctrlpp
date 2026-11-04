#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NLOPT_INTEROP_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NLOPT_INTEROP_H

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nlopt_types.h"

#include <nlopt.hpp>

#include <new>
#include <vector>
#include <cstddef>
#include <utility>
#include <exception>

#if defined(__GLIBCXX__)
    #include <cxxabi.h>
#endif

namespace ctrlpp::detail
{

inline nlopt::algorithm to_nlopt_algorithm(nlopt_algorithm alg)
{
    switch(alg)
    {
        case nlopt_algorithm::slsqp:
            return nlopt::LD_SLSQP;
        case nlopt_algorithm::mma:
            return nlopt::LD_MMA;
        case nlopt_algorithm::cobyla:
            return nlopt::LN_COBYLA;
        case nlopt_algorithm::isres:
            return nlopt::GN_ISRES;
        case nlopt_algorithm::auglag_mma:
            return nlopt::AUGLAG_EQ;
        case nlopt_algorithm::auglag_ccsaq:
            return nlopt::AUGLAG_EQ;
        case nlopt_algorithm::ccsaq:
            return nlopt::LD_CCSAQ;
    }
    return nlopt::LD_SLSQP;
}

inline nlopt::algorithm to_nlopt_inner(nlopt_algorithm alg)
{
    switch(alg)
    {
        case nlopt_algorithm::auglag_mma:
            return nlopt::LD_MMA;
        case nlopt_algorithm::auglag_ccsaq:
            return nlopt::LD_CCSAQ;
        default:
            return nlopt::LD_SLSQP;
    }
}

constexpr solve_status map_nlopt_result(nlopt::result rc)
{
    switch(rc)
    {
        case nlopt::SUCCESS:
            [[fallthrough]];
        case nlopt::STOPVAL_REACHED:
            [[fallthrough]];
        case nlopt::FTOL_REACHED:
            [[fallthrough]];
        case nlopt::XTOL_REACHED:
            return solve_status::optimal;
        case nlopt::MAXEVAL_REACHED:
            return solve_status::max_iterations;
        case nlopt::MAXTIME_REACHED:
            return solve_status::time_limit;
        case nlopt::ROUNDOFF_LIMITED:
            return solve_status::solved_inaccurate;
        default:
            return solve_status::error;
    }
}

// A forced unwind (glibc thread cancellation) is re-raised; one raised inside optimize is still ended by NLopt's callback wrapper.
// https://gcc.gnu.org/onlinedocs/libstdc++/manual/using_exceptions.html#intro.using.exception.compat
template <typename Result, typename Fn>
inline Result nlopt_contained(Fn&& fn, Result on_allocation_failure, Result on_other_failure)
{
    try
    {
        return std::forward<Fn>(fn)();
    }
#if defined(__GLIBCXX__)
    catch(const abi::__forced_unwind&)
    {
        throw;
    }
#endif
    catch(const std::bad_alloc&)
    {
        return on_allocation_failure;
    }
    catch(const std::exception&)
    {
        return on_other_failure;
    }
    catch(...)
    {
        return on_other_failure;
    }
}

template <typename Scalar>
inline std::vector<double> to_nlopt_vector(const Eigen::VectorX<Scalar>& v)
{
    std::vector<double> out(static_cast<std::size_t>(v.size()));
    for(Eigen::Index i = 0; i < v.size(); ++i)
    {
        out[static_cast<std::size_t>(i)] = static_cast<double>(v[i]);
    }
    return out;
}

template <typename Scalar>
inline void set_nlopt_variable_bounds(nlopt::opt& optimizer, const nlp_problem<Scalar>& problem)
{
    if(problem.x_lower.size() > 0)
        optimizer.set_lower_bounds(to_nlopt_vector(problem.x_lower));
    if(problem.x_upper.size() > 0)
        optimizer.set_upper_bounds(to_nlopt_vector(problem.x_upper));
}

template <typename Scalar>
inline void set_nlopt_stopping_criteria(nlopt::opt& optimizer, const nlopt_settings<Scalar>& settings)
{
    optimizer.set_ftol_rel(static_cast<double>(settings.ftol_rel));
    optimizer.set_xtol_rel(static_cast<double>(settings.xtol_rel));
    optimizer.set_maxeval(settings.max_eval);
    if(settings.max_time > Scalar{0})
        optimizer.set_maxtime(static_cast<double>(settings.max_time));
}

// NLopt refuses a start outside the bounds, so a valid box that excludes the previous
// iterate would otherwise be refused. The comparisons leave inverted or NaN bounds
// unchanged for NLopt to refuse, which std::clamp's lower <= upper precondition would not.
template <typename Scalar>
inline void project_onto_nlopt_bounds(std::vector<double>& x, const nlp_problem<Scalar>& problem)
{
    const auto n = static_cast<Eigen::Index>(x.size());
    const bool has_lower = problem.x_lower.size() == n;
    const bool has_upper = problem.x_upper.size() == n;
    for(Eigen::Index i = 0; i < n; ++i)
    {
        double& xi = x[static_cast<std::size_t>(i)];
        if(has_upper && xi > static_cast<double>(problem.x_upper[i]))
            xi = static_cast<double>(problem.x_upper[i]);
        if(has_lower && xi < static_cast<double>(problem.x_lower[i]))
            xi = static_cast<double>(problem.x_lower[i]);
    }
}

template <typename Scalar>
inline void set_nlopt_local_optimizer(nlopt::opt& optimizer, const nlopt_settings<Scalar>& settings, unsigned n_vars)
{
    nlopt::opt local(to_nlopt_inner(settings.algorithm), n_vars);
    set_nlopt_stopping_criteria(local, settings);
    optimizer.set_local_optimizer(local);
}

template <typename Scalar>
inline nlp_result<Scalar> failed_nlopt_result()
{
    nlp_result<Scalar> result{};
    result.status = solve_status::error;
    return result;
}

}

#endif
