#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_SOLVER_H

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/argmin_problem.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <nablapp/result/status.h>
#include <nablapp/solver/options.h>
#include <nablapp/solver/basic_solver.h>

#include <cmath>
#include <chrono>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace ctrlpp
{

template <typename Scalar, typename Policy>
class argmin_solver
{
public:
    using scalar_type = Scalar;
    using nablapp_policy = typename Policy::algorithm;
    using solver_type = nablapp::basic_solver<nablapp_policy, Eigen::Dynamic, argmin_problem<Scalar>>;

    explicit argmin_solver(argmin_settings<Scalar> settings = {})
        : settings_{settings}
    {}

    void setup(const nlp_problem<Scalar>& problem)
    {
        problem_ = &problem;
        solver_ = std::nullopt;
        bridge_.partition(problem);

        if constexpr(std::is_same_v<Policy, argmin_mma>)
        {
            if(bridge_.num_equality() > 0)
                throw std::invalid_argument("MMA algorithm does not support equality constraints");
        }
    }

    auto solve(const nlp_update<Scalar>& update) -> nlp_result<Scalar>
    {
        prepare_solver(update.x0);
        auto result = solver_->solve();
        return translate_result(result);
    }

    auto step(const nlp_update<Scalar>& update, int max_steps) -> nlp_result<Scalar>
    {
        prepare_solver(update.x0);
        auto result = solver_->step_n(static_cast<std::uint32_t>(max_steps));
        return translate_result(result);
    }

private:
    void prepare_solver(const Eigen::VectorX<Scalar>& x0)
    {
        if(!solver_)
        {
            solver_.emplace(nablapp_policy{}, bridge_, x0, make_solver_options());
        }
        else
        {
            if(settings_.warm_start == warm_start_mode::curvature)
                solver_->reset(x0);
            else
                solver_->reset_clear(x0);
        }
    }

    auto make_solver_options() const -> nablapp::solver_options<>
    {
        nablapp::solver_options<> opts;
        opts.max_iterations = static_cast<std::uint32_t>(settings_.max_eval);

        if(settings_.max_time > Scalar{0})
        {
            opts.max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::duration<double>(static_cast<double>(settings_.max_time)));
        }

        if(settings_.constraint_tol > Scalar{0})
            opts.constraint_tolerance = static_cast<double>(settings_.constraint_tol);

        opts.set_objective_threshold(static_cast<double>(settings_.ftol_rel));
        opts.set_step_threshold(static_cast<double>(settings_.xtol_rel));

        return opts;
    }

    static constexpr auto map_status(nablapp::solver_status s) -> solve_status
    {
        switch(s)
        {
        case nablapp::solver_status::converged:
        case nablapp::solver_status::ftol_reached:
        case nablapp::solver_status::xtol_reached:
            return solve_status::optimal;
        case nablapp::solver_status::max_iterations:
        case nablapp::solver_status::budget_exhausted:
        case nablapp::solver_status::maxeval_reached:
            return solve_status::max_iterations;
        case nablapp::solver_status::time_limit_reached:
            return solve_status::time_limit;
        case nablapp::solver_status::stalled:
        case nablapp::solver_status::roundoff_limited:
        case nablapp::solver_status::objective_stalled:
            return solve_status::solved_inaccurate;
        case nablapp::solver_status::diverged:
        case nablapp::solver_status::aborted:
        case nablapp::solver_status::running:
            return solve_status::error;
        }
        return solve_status::error;
    }

    auto translate_result(const nablapp::solve_result<Scalar, Eigen::Dynamic>& r) const -> nlp_result<Scalar>
    {
        return nlp_result<Scalar>{
            .status = map_status(r.status),
            .x = r.x,
            .objective = r.objective_value,
            .solve_time = static_cast<Scalar>(std::chrono::duration<double>(r.wall_time).count()),
            .iterations = static_cast<int>(r.iterations),
            .primal_residual = r.constraint_violation,
        };
    }

    argmin_settings<Scalar> settings_;
    const nlp_problem<Scalar>* problem_{nullptr};
    argmin_problem<Scalar> bridge_;
    std::optional<solver_type> solver_;
};

}

#endif
