#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_SOLVER_H

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/argmin_problem.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <argmin/result/status.h>

#include <argmin/solver/options.h>
#include <argmin/solver/basic_solver.h>

#include <chrono>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <type_traits>

namespace ctrlpp
{

template <typename Scalar, typename Policy, bool Constrained = true>
class argmin_solver
{
public:
    using scalar_type = Scalar;
    using argmin_policy = typename Policy::algorithm;
    using bridge_type = std::conditional_t<Constrained,
        argmin_constrained_problem<Scalar>,
        argmin_problem<Scalar>>;
    using solver_type = argmin::basic_solver<argmin_policy, Eigen::Dynamic, bridge_type>;
    using settings_type = std::conditional_t<
        is_mma_family_v<Policy>,
        argmin_mma_settings<Scalar>,
        argmin_settings<Scalar>>;

    explicit argmin_solver(settings_type settings = {})
        : settings_{settings}
    {}

    void setup(const nlp_problem<Scalar>& problem)
    {
        problem_ = &problem;
        solver_ = std::nullopt;

        if constexpr(Constrained)
        {
            bridge_.partition(problem);

        }
        else
        {
            bridge_.bind(problem);
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
            if constexpr(is_mma_family_v<Policy>)
            {
                solver_.emplace(argmin_policy{}, bridge_, x0,
                                make_solver_options(),
                                make_mma_policy_opts());
            }
            else
            {
                solver_.emplace(argmin_policy{}, bridge_, x0, make_solver_options());
            }
        }
        else
        {
            const auto ws = [&]() -> warm_start_mode
            {
                if constexpr(is_mma_family_v<Policy>)
                    return settings_.base.warm_start;
                else
                    return settings_.warm_start;
            }();
            if(ws == warm_start_mode::curvature)
                solver_->reset(x0);
            else
                solver_->reset_clear(x0);
        }
    }

    auto make_solver_options() const -> argmin::solver_options<>
    {
        argmin::solver_options<> opts;

        auto const& s = [&]() -> auto const&
        {
            if constexpr(is_mma_family_v<Policy>)
                return settings_.base;
            else
                return settings_;
        }();

        opts.max_iterations = static_cast<std::uint32_t>(s.max_eval);

        if(s.max_time > Scalar{0})
        {
            opts.max_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::duration<double>(static_cast<double>(s.max_time)));
        }

        if constexpr(Constrained)
        {
            if(s.constraint_tol > Scalar{0})
                opts.constraint_tolerance = static_cast<double>(s.constraint_tol);
        }

        opts.set_objective_threshold(static_cast<double>(s.ftol_rel));
        opts.set_step_threshold(static_cast<double>(s.xtol_rel));

        return opts;
    }

    auto make_mma_policy_opts() const -> typename argmin_policy::options_type
    {
        typename argmin_policy::options_type po{};
        if constexpr(std::is_same_v<Policy, argmin_mma>)
        {
            po.asymptote_init = static_cast<double>(settings_.asymptote_init);
            po.asymptote_expand = static_cast<double>(settings_.asymptote_incr);
            po.asymptote_contract = static_cast<double>(settings_.asymptote_decr);
        }
        else if constexpr(std::is_same_v<Policy, argmin_gcmma>)
        {
            po.asymptote_init = static_cast<double>(settings_.asymptote_init);
            po.asymptote_expand = static_cast<double>(settings_.asymptote_incr);
            po.asymptote_contract = static_cast<double>(settings_.asymptote_decr);
            po.max_inner_iterations =
                static_cast<std::uint16_t>(settings_.gcmma_inner_max);
        }
        else
        {
            // argmin_auglag<Inner> with is_mma_family_v<Inner>.
            if constexpr(std::is_same_v<Policy, argmin_auglag<argmin_mma>>)
            {
                po.inner_opts.asymptote_init =
                    static_cast<double>(settings_.asymptote_init);
                po.inner_opts.asymptote_expand =
                    static_cast<double>(settings_.asymptote_incr);
                po.inner_opts.asymptote_contract =
                    static_cast<double>(settings_.asymptote_decr);
            }
            else if constexpr(std::is_same_v<Policy, argmin_auglag<argmin_gcmma>>)
            {
                po.inner_opts.asymptote_init =
                    static_cast<double>(settings_.asymptote_init);
                po.inner_opts.asymptote_expand =
                    static_cast<double>(settings_.asymptote_incr);
                po.inner_opts.asymptote_contract =
                    static_cast<double>(settings_.asymptote_decr);
                po.inner_opts.max_inner_iterations =
                    static_cast<std::uint16_t>(settings_.gcmma_inner_max);
            }
        }
        return po;
    }

    static constexpr auto map_status(argmin::solver_status s) -> solve_status
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
            return solve_status::solved_inaccurate;
        case argmin::solver_status::diverged:
        case argmin::solver_status::aborted:
        case argmin::solver_status::running:
            return solve_status::error;
        }
        return solve_status::error;
    }

    auto translate_result(const argmin::solve_result<Scalar, Eigen::Dynamic>& r) const -> nlp_result<Scalar>
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

    settings_type settings_;
    const nlp_problem<Scalar>* problem_{nullptr};
    bridge_type bridge_;
    std::optional<solver_type> solver_;
};

}

#endif
