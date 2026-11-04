#ifndef HPP_GUARD_CTRLPP_MPC_NLOPT_SOLVER_H
#define HPP_GUARD_CTRLPP_MPC_NLOPT_SOLVER_H

// NLopt's C++ API throws by upstream design. setup and solve contain every exception, so this opt-in backend needs
// exception support and is exempt from the library's no-exceptions floor.

#include "ctrlpp/mpc/nlp_types.h"
#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nlopt_types.h"

#include "ctrlpp/mpc/detail/nlopt_interop.h"
#include "ctrlpp/mpc/detail/nlopt_constraint_rows.h"

#include "ctrlpp/expected.h"

#include <nlopt.hpp>

#include <span>
#include <chrono>
#include <vector>
#include <cstddef>
#include <utility>
#include <type_traits>

namespace ctrlpp
{

template <typename Scalar>
class nlopt_solver
{
    static_assert(std::is_same_v<Scalar, double>, "NLopt operates in double precision only");

public:
    using scalar_type = Scalar;

    explicit nlopt_solver(nlopt_settings<Scalar> settings = {})
        : m_rows{}
        , m_problem{nullptr}
        , m_settings{settings}
        , m_optimizer{}
        , m_ready{false}
    {
    }

    // NLopt holds every callback by the address of this solver or of its constraint rows, so a copy or a move re-runs
    // setup against the same referenced problem, which the caller keeps at a stable address for the solver's lifetime.
    nlopt_solver(const nlopt_solver& other)
        : m_rows{}
        , m_problem{other.m_problem}
        , m_settings{other.m_settings}
        , m_optimizer{}
        , m_ready{false}
    {
        rebind_callbacks();
    }

    nlopt_solver(nlopt_solver&& other)
        : nlopt_solver(std::as_const(other))
    {
    }

    nlopt_solver& operator=(const nlopt_solver& other)
    {
        if(this != &other)
        {
            m_settings = other.m_settings;
            m_problem = other.m_problem;
            rebind_callbacks();
        }
        return *this;
    }

    nlopt_solver& operator=(nlopt_solver&& other) { return *this = std::as_const(other); }

    ~nlopt_solver() = default;

    ctrlpp::expected<void, nlopt_setup_error> setup(const nlp_problem<Scalar>& problem)
    {
        using outcome = ctrlpp::expected<void, nlopt_setup_error>;
        m_ready = false;
        m_problem = &problem;
        const outcome configured = detail::nlopt_contained<outcome>(
                [&] { return configure(problem); }, ctrlpp::unexpected(nlopt_setup_error::allocation_failed), ctrlpp::unexpected(nlopt_setup_error::backend_rejected));
        m_ready = configured.has_value();
        return configured;
    }

    nlp_result<Scalar> solve(const nlp_update<Scalar>& update)
    {
        const auto t0 = std::chrono::steady_clock::now();
        nlp_result<Scalar> result = m_ready ? guarded_solve(update) : detail::failed_nlopt_result<Scalar>();
        result.solve_time = static_cast<Scalar>(std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        return result;
    }

private:
    detail::nlopt_constraint_rows<Scalar> m_rows;
    // Null until setup(). An optional reference would widen the solver and grow the frame that builds an estimator around it.
    const nlp_problem<Scalar>* m_problem;
    nlopt_settings<Scalar> m_settings;
    nlopt::opt m_optimizer;
    bool m_ready;

    ctrlpp::expected<void, nlopt_setup_error> configure(const nlp_problem<Scalar>& problem)
    {
        const bool rejects_equalities = m_settings.algorithm == nlopt_algorithm::mma || m_settings.algorithm == nlopt_algorithm::ccsaq;
        if(auto partitioned = m_rows.partition(problem, rejects_equalities); !partitioned.has_value())
            return partitioned;

        m_optimizer = nlopt::opt(detail::to_nlopt_algorithm(m_settings.algorithm), static_cast<unsigned>(problem.n_vars));
        detail::set_nlopt_variable_bounds(m_optimizer, problem);
        m_optimizer.set_min_objective(objective_callback, this);
        m_rows.register_callbacks(m_optimizer, static_cast<double>(m_settings.constraint_tol));
        detail::set_nlopt_stopping_criteria(m_optimizer, m_settings);
        if(m_settings.algorithm == nlopt_algorithm::auglag_mma || m_settings.algorithm == nlopt_algorithm::auglag_ccsaq)
            detail::set_nlopt_local_optimizer(m_optimizer, m_settings, static_cast<unsigned>(problem.n_vars));
        return {};
    }

    // solve reports every failure through nlp_result::status, so no exception may cross it; run_optimizer's
    // round-off cost re-evaluation and the violation evaluation in solved_result run inside this guard.
    nlp_result<Scalar> guarded_solve(const nlp_update<Scalar>& update)
    {
        return detail::nlopt_contained<nlp_result<Scalar>>(
                [&]
                {
                    auto x = detail::to_nlopt_vector(update.x0);
                    detail::project_onto_nlopt_bounds(x, *m_problem);
                    double opt_f{};
                    const solve_status status = run_optimizer(x, opt_f);
                    return solved_result(x, opt_f, status);
                },
                detail::failed_nlopt_result<Scalar>(), detail::failed_nlopt_result<Scalar>());
    }

    solve_status run_optimizer(std::vector<double>& x, double& opt_f)
    {
        try
        {
            return detail::map_nlopt_result(m_optimizer.optimize(x, opt_f));
        }
        catch(const nlopt::roundoff_limited&)
        {
            opt_f = m_problem->cost_cb({x.data(), x.size()});
            return solve_status::solved_inaccurate;
        }
    }

    nlp_result<Scalar> solved_result(const std::vector<double>& x, double opt_f, solve_status status) const
    {
        nlp_result<Scalar> result{};
        result.status = status;
        result.x.resize(static_cast<Eigen::Index>(x.size()));
        for(Eigen::Index i = 0; i < result.x.size(); ++i)
            result.x[i] = static_cast<Scalar>(x[static_cast<std::size_t>(i)]);
        result.objective = static_cast<Scalar>(opt_f);
        result.iterations = m_optimizer.get_numevals();
        result.primal_residual = m_rows.max_violation(x);
        return result;
    }

    // The re-setup fails when the referenced problem changed since the source's setup. A copy or move cannot report
    // that, so it stays unready and its next solve reports an error, as does a solver whose source was never set up.
    void rebind_callbacks()
    {
        m_ready = false;
        if(m_problem != nullptr)
            static_cast<void>(setup(*m_problem));
    }

    static double objective_callback(const std::vector<double>& x, std::vector<double>& grad, void* data)
    {
        auto* self = static_cast<nlopt_solver*>(data);
        std::span<const double> x_span{x.data(), x.size()};

        if(!grad.empty() && self->m_problem->gradient_cb)
        {
            std::span<double> grad_span{grad.data(), grad.size()};
            self->m_problem->gradient_cb(x_span, grad_span);
        }

        return self->m_problem->cost_cb(x_span);
    }
};

}

#endif
