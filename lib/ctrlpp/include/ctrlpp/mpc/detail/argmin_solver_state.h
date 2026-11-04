#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_ARGMIN_SOLVER_STATE_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_ARGMIN_SOLVER_STATE_H

#include "ctrlpp/mpc/detail/argmin_options.h"

#include "ctrlpp/mpc/argmin_problem.h"
#include "ctrlpp/mpc/argmin_policies.h"

#include <argmin/solver/step_budget_solver.h>
#include <argmin/solver/step_and_time_budget_solver.h>

#include <memory>
#include <variant>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename Scalar, typename Policy, bool Constrained, int NV, int MaxM, argmin_budget Budget>
class argmin_solver_state
{
protected:
    static constexpr bool time_budgeted = std::is_same_v<Budget, argmin_step_and_time_budget>;

public:
    using argmin_policy = rebind_argmin_algorithm_t<typename Policy::algorithm, NV>;
    using bridge_type = std::conditional_t<Constrained, argmin_constrained_problem<Scalar, NV, MaxM>, argmin_problem<Scalar, NV, MaxM>>;
    using problem_type = typename bridge_type::problem_type;
    using step_solver_type = argmin::step_budget_solver<argmin_policy, NV, bridge_type, argmin_ctrlpp_convergence>;
    using timed_solver_type = argmin::step_and_time_budget_solver<argmin_policy, NV, bridge_type, argmin_ctrlpp_convergence>;
    using solver_storage_type = std::conditional_t<time_budgeted, std::variant<std::monostate, step_solver_type, timed_solver_type>, std::variant<std::monostate, step_solver_type>>;
    using settings_type = typename Budget::template settings_type<Scalar>;

    explicit argmin_solver_state(settings_type settings)
        : m_settings{settings}
        , m_problem{nullptr}
        , m_bridge{std::make_unique<bridge_type>()}
    {
    }

    // The bridge sits behind a unique_ptr, so a defaulted move keeps valid the references to it that argmin's solver
    // core and policy state hold.
    argmin_solver_state(argmin_solver_state&&) = default;
    argmin_solver_state& operator=(argmin_solver_state&&) = default;

    // argmin's solver core is not copyable, so a copy binds a fresh bridge to the same problem and leaves the solver
    // empty; the next solve or step emplaces it.
    argmin_solver_state(const argmin_solver_state& other)
        : m_settings{other.m_settings}
        , m_problem{other.m_problem}
        , m_bridge{std::make_unique<bridge_type>()}
        , m_solver{}
    {
        rebind_bridge();
    }

    argmin_solver_state& operator=(const argmin_solver_state& other)
    {
        if(this != &other)
        {
            auto fresh = std::make_unique<bridge_type>();
            m_settings = other.m_settings;
            m_problem = other.m_problem;
            // The stale solver caches the old bridge's address, so it goes before the bridge does.
            m_solver.template emplace<std::monostate>();
            m_bridge = std::move(fresh);
            rebind_bridge();
        }
        return *this;
    }

protected:
    settings_type m_settings;
    // Null until setup(). An optional reference would widen the solver, which grows the solve frames of the controllers that hold it.
    const problem_type* m_problem;
    // Declared before m_solver, which caches the bridge's address, so the solver is destroyed first.
    std::unique_ptr<bridge_type> m_bridge;
    solver_storage_type m_solver;

    void rebind_bridge()
    {
        if(m_problem == nullptr)
            return;

        if constexpr(Constrained)
            m_bridge->partition(*m_problem);
        else
            m_bridge->bind(*m_problem);
    }

    void prepare_solver(const Eigen::VectorX<Scalar>& x0)
    {
        if(std::holds_alternative<std::monostate>(m_solver))
        {
            emplace_solver(x0);
            return;
        }

        const auto ws = m_settings.warm_start;

        std::visit(
                [&](auto& solver)
                {
                    using solver_t = std::decay_t<decltype(solver)>;
                    if constexpr(!std::is_same_v<solver_t, std::monostate>)
                    {
                        if(ws == warm_start_mode::curvature)
                            solver.reset(x0);
                        else
                            solver.reset_clear(x0);
                    }
                },
                m_solver);
    }

    void emplace_solver(const Eigen::VectorX<Scalar>& x0)
    {
        if constexpr(time_budgeted)
        {
            if(m_settings.max_time > Scalar{0})
            {
                emplace_timed_solver(x0);
                return;
            }
        }
        emplace_step_solver(x0);
    }

    void emplace_step_solver(const Eigen::VectorX<Scalar>& x0)
    {
        m_solver.template emplace<step_solver_type>(argmin_policy{}, *m_bridge, x0, make_argmin_solver_options<Constrained>(m_settings));
    }

    void emplace_timed_solver(const Eigen::VectorX<Scalar>& x0)
    {
        m_solver.template emplace<timed_solver_type>(argmin_policy{}, *m_bridge, x0, make_argmin_time_budget_options<Constrained>(m_settings));
    }
};

}

#endif
