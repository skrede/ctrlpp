#ifndef HPP_GUARD_CTRLPP_MPC_ARGMIN_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_ARGMIN_SOLVER_H

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/argmin_problem.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include "ctrlpp/mpc/detail/argmin_options.h"
#include "ctrlpp/mpc/detail/argmin_solver_state.h"
#include "ctrlpp/mpc/detail/argmin_setup_admission.h"

#include "ctrlpp/expected.h"

#include <argmin/result/status.h>

#include <chrono>
#include <cstdint>
#include <variant>
#include <stdexcept>
#include <type_traits>

namespace ctrlpp
{

// A fixed NV rebinds the policy and the bridge to fixed-size decision storage and binds an nlp_problem_static; a fixed
// MaxM also makes argmin's result multipliers inline. Binding both gives the allocation-free steady-state solve.
template <typename Scalar, typename Policy, bool Constrained = true, int NV = Eigen::Dynamic, int MaxM = Eigen::Dynamic, argmin_budget Budget = argmin_step_budget>
class argmin_solver : private detail::argmin_solver_state<Scalar, Policy, Constrained, NV, MaxM, Budget>
{
    using base = detail::argmin_solver_state<Scalar, Policy, Constrained, NV, MaxM, Budget>;

public:
    using scalar_type = Scalar;
    using argmin_policy = typename base::argmin_policy;
    using bridge_type = typename base::bridge_type;
    using problem_type = typename base::problem_type;
    using step_solver_type = typename base::step_solver_type;
    using timed_solver_type = typename base::timed_solver_type;
    using solver_storage_type = typename base::solver_storage_type;
    using settings_type = typename base::settings_type;

    // Inherited rather than forwarded: a forwarding constructor leaves a copy of the settings in the creating frame.
    using base::base;

    explicit argmin_solver()
        : base{settings_type{}}
    {
    }

    ctrlpp::expected<void, argmin_setup_error> setup(const problem_type& problem)
    {
        if(const auto refusal = detail::argmin_setup_refusal<Constrained, NV, MaxM>(problem))
            return ctrlpp::unexpected(*refusal);

        this->m_problem = &problem;
        this->m_solver.template emplace<std::monostate>();

        if constexpr(Constrained)
            this->m_bridge->partition(problem);
        else
            this->m_bridge->bind(problem);

        return {};
    }

    nlp_result<Scalar> solve(const nlp_update<Scalar>& update)
    {
        this->prepare_solver(update.x0);
        return with_solver(
                [&](auto& solver)
                {
                    auto result = solver.solve();
                    return translate_result(result);
                });
    }

    // Writes the primal into the caller's solution, which must already have the problem's dimension, and leaves the
    // result's x empty, so the static path allocates no decision vector.
    nlp_result<Scalar> solve_into(const nlp_update<Scalar>& update, Eigen::Ref<Eigen::VectorX<Scalar>> solution)
    {
        this->prepare_solver(update.x0);
        return with_solver(
                [&](auto& solver)
                {
                    auto result = solver.solve();
                    solution = result.x;
                    return translate_meta(result);
                });
    }

    nlp_result<Scalar> step(const nlp_update<Scalar>& update, std::int32_t max_steps)
    {
        this->prepare_solver(update.x0);
        return with_solver(
                [&](auto& solver)
                {
                    auto result = solver.step_n(static_cast<std::uint32_t>(max_steps));
                    return translate_result(result);
                });
    }

private:
    template <typename Result>
    nlp_result<Scalar> translate_meta(const Result& r) const
    {
        nlp_result<Scalar> out;
        out.status = detail::map_argmin_status(r.status);
        out.objective = r.objective_value;
        out.solve_time = result_wall_time(r);
        out.iterations = static_cast<std::int32_t>(r.iterations);
        out.primal_residual = r.constraint_violation;
        out.stop_criterion = detail::map_argmin_stop_criterion(r.status);
        return out;
    }

    template <typename Result>
    nlp_result<Scalar> translate_result(const Result& r) const
    {
        auto out = translate_meta(r);
        out.x = r.x;
        return out;
    }

    template <typename Result>
    static Scalar result_wall_time(const Result& r)
    {
        if constexpr(requires { r.wall_time; })
        {
            return static_cast<Scalar>(std::chrono::duration<double>(r.wall_time).count());
        }
        else
        {
            return Scalar{0};
        }
    }

    template <typename SolverFn>
    nlp_result<Scalar> with_solver(SolverFn&& fn)
    {
        return std::visit(
                [&](auto& solver) -> nlp_result<Scalar>
                {
                    using solver_t = std::decay_t<decltype(solver)>;
                    if constexpr(std::is_same_v<solver_t, std::monostate>)
                    {
                        nlp_result<Scalar> failed;
                        failed.status = solve_status::error;
                        return failed;
                    }
                    else
                    {
                        return fn(solver);
                    }
                },
                this->m_solver);
    }
};

}

#endif
