#ifndef HPP_GUARD_CTRLPP_NMPC_H
#define HPP_GUARD_CTRLPP_NMPC_H

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/model/dynamics_model.h"
#include "ctrlpp/mpc/nlp_formulation.h"

#include <Eigen/Dense>

#include <span>
#include <cmath>
#include <memory>
#include <vector>
#include <cstddef>
#include <utility>
#include <optional>
#include <algorithm>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX, std::size_t NU, nlp_solver Solver, dynamics_model<Scalar, NX, NU> Dynamics, std::size_t NC = 0, std::size_t NTC = 0>
class nmpc
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int ntc = static_cast<int>(NTC);

public:
    nmpc(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config)
        : nmpc{std::move(dynamics), config, Solver{}}
    {}

    nmpc(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, Solver solver)
        : m_dynamics{std::move(dynamics)}
        , m_config{config}
        , m_N{config.horizon}
        , m_has_path_slack{config.soft_constraints && NC > 0 && config.path_constraint.has_value()}
        , m_has_term_slack{config.soft_constraints && NTC > 0 && config.terminal_constraint.has_value()}
        , m_num_path_slack{m_has_path_slack ? m_N * nc : 0}
        , m_num_term_slack{m_has_term_slack ? ntc : 0}
        , m_num_vars{(m_N + 1) * nx + m_N * nu + m_num_path_slack + m_num_term_slack}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>()}
        , m_solver{std::move(solver)}
    {
        m_state->x_ref.resize(static_cast<std::size_t>(m_N + 1), Vector<Scalar, NX>::Zero());
        m_problem = std::make_unique<nlp_problem<Scalar>>(
            detail::build_nmpc_problem<Scalar, NX, NU, NC, NTC>(m_dynamics, m_config, m_state));
        m_setup_failed = !detail::setup_nlp_solver(m_solver, *m_problem);
        m_warm_z = Eigen::VectorX<Scalar>::Zero(m_num_vars);
    }

    // Move is correct-by-default: `m_problem` lives behind a `unique_ptr` (a
    // stable heap address), so moving relocates only the owning pointer while the
    // pointed-to problem stays put. The solver's bridge caches `&m_problem` by
    // reference (argmin_problem.h:24-27), so that cached pointer stays valid
    // across a defaulted move. `m_state` is already a heap-stable shared_ptr.
    nmpc(nmpc&&) = default;
    nmpc& operator=(nmpc&&) = default;
    ~nmpc() = default;

    // Copy is an independent fork (the solver's copy is copy-deleted-by-argmin
    // and re-emplaced lazily, so aliasing would be unsafe anyway): deep-copy the
    // formulation state, rebuild the problem bound to that fresh state, then
    // re-setup the forked solver against the new problem. The two controllers
    // share no mutable state and solve independently.
    nmpc(const nmpc& other)
        : m_dynamics{other.m_dynamics}
        , m_config{other.m_config}
        , m_N{other.m_N}
        , m_has_path_slack{other.m_has_path_slack}
        , m_has_term_slack{other.m_has_term_slack}
        , m_num_path_slack{other.m_num_path_slack}
        , m_num_term_slack{other.m_num_term_slack}
        , m_num_vars{other.m_num_vars}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>(*other.m_state)}
        , m_problem{std::make_unique<nlp_problem<Scalar>>(
              detail::build_nmpc_problem<Scalar, NX, NU, NC, NTC>(m_dynamics, m_config, m_state))}
        , m_solver{other.m_solver}
        , m_warm_z{other.m_warm_z}
        , m_last_solution{other.m_last_solution}
        , m_last_diagnostics{other.m_last_diagnostics}
        , m_u_prev{other.m_u_prev}
    {
        m_setup_failed = !detail::setup_nlp_solver(m_solver, *m_problem);
    }

    nmpc& operator=(const nmpc& other)
    {
        nmpc tmp{other};
        *this = std::move(tmp);
        return *this;
    }

    std::optional<Vector<Scalar, NU>> solve(const Vector<Scalar, NX>& x0)
    {
        for(auto& ref : m_state->x_ref)
            ref.setZero();
        return solve_impl(x0);
    }

    std::optional<Vector<Scalar, NU>> solve(const Vector<Scalar, NX>& x0, const Vector<Scalar, NX>& x_ref)
    {
        for(auto& ref : m_state->x_ref)
            ref = x_ref;
        return solve_impl(x0);
    }

    std::optional<Vector<Scalar, NU>> solve(const Vector<Scalar, NX>& x0, std::span<const Vector<Scalar, NX>> x_ref)
    {
        const auto len = std::min(x_ref.size(), m_state->x_ref.size());
        for(std::size_t k = 0; k < len; ++k)
            m_state->x_ref[k] = x_ref[k];
        if(!x_ref.empty())
            for(std::size_t k = len; k < m_state->x_ref.size(); ++k)
                m_state->x_ref[k] = x_ref.back();
        return solve_impl(x0);
    }

    std::pair<std::vector<Vector<Scalar, NX>>, std::vector<Vector<Scalar, NU>>> trajectory() const
    {
        std::vector<Vector<Scalar, NX>> states;
        std::vector<Vector<Scalar, NU>> inputs;
        states.reserve(static_cast<std::size_t>(m_N + 1));
        inputs.reserve(static_cast<std::size_t>(m_N));

        for(int k = 0; k <= m_N; ++k)
            states.push_back(m_last_solution.segment(k * nx, nx));
        const int u_offset = (m_N + 1) * nx;
        for(int k = 0; k < m_N; ++k)
            inputs.push_back(m_last_solution.segment(u_offset + k * nu, nu));

        return {std::move(states), std::move(inputs)};
    }

    mpc_diagnostics<Scalar> diagnostics() const { return m_last_diagnostics; }
    const nlp_problem<Scalar>& problem() const { return *m_problem; }
    const Eigen::VectorX<Scalar>& last_solution() const { return m_last_solution; }

private:
    std::optional<Vector<Scalar, NU>> solve_impl(const Vector<Scalar, NX>& x0)
    {
        if(m_setup_failed)
            return std::nullopt;

        m_state->x0 = x0;
        m_state->u_prev = m_u_prev;

        nlp_update<Scalar> update;
        update.x0 = m_warm_z;

        auto result = m_solver.solve(update);
        populate_diagnostics(result);

        if(result.status != solve_status::optimal && result.status != solve_status::solved_inaccurate)
            return std::nullopt;

        m_last_solution = result.x;
        populate_constraint_diagnostics(result.x);
        shift_warm_start(result.x);
        return extract_first_input(result.x);
    }

    void populate_diagnostics(const nlp_result<Scalar>& result)
    {
        m_last_diagnostics = mpc_diagnostics<Scalar>{.status = result.status,
                                                     .iterations = result.iterations,
                                                     .solve_time = result.solve_time,
                                                     .cost = result.objective,
                                                     .primal_residual = result.primal_residual,
                                                     .dual_residual = Scalar{0},
                                                     .max_constraint_violation = result.primal_residual,
                                                     .max_path_constraint_violation = Scalar{0},
                                                     .max_terminal_constraint_violation = Scalar{0},
                                                     .total_slack = Scalar{0}};
    }

    void populate_constraint_diagnostics(const Eigen::VectorX<Scalar>& z)
    {
        const int u_offset = (m_N + 1) * nx;
        compute_path_constraint_violation(z, u_offset);
        compute_terminal_constraint_violation(z);
        compute_total_slack(z, u_offset);
    }

    void compute_path_constraint_violation(const Eigen::VectorX<Scalar>& z, int u_offset)
    {
        if constexpr(NC > 0)
        {
            if(m_config.path_constraint)
            {
                const auto& g = *m_config.path_constraint;
                Scalar max_viol{0};

                for(int k = 0; k < m_N; ++k)
                {
                    Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + k * nx);
                    Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);
                    Vector<Scalar, NC> gk = g(xk, uk);

                    for(int i = 0; i < nc; ++i)
                        max_viol = std::max(max_viol, gk[i]);
                }

                m_last_diagnostics.max_path_constraint_violation = max_viol;
            }
        }
    }

    void compute_terminal_constraint_violation(const Eigen::VectorX<Scalar>& z)
    {
        if constexpr(NTC > 0)
        {
            if(m_config.terminal_constraint)
            {
                const auto& h = *m_config.terminal_constraint;
                Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + m_N * nx);
                Vector<Scalar, NTC> hN = h(xN);

                Scalar max_viol{0};
                for(int i = 0; i < ntc; ++i)
                    max_viol = std::max(max_viol, hN[i]);

                m_last_diagnostics.max_terminal_constraint_violation = max_viol;
            }
        }
    }

    void compute_total_slack(const Eigen::VectorX<Scalar>& z, int u_offset)
    {
        Scalar total{0};
        if(m_has_path_slack)
        {
            const int slack_off = u_offset + m_N * nu;
            for(int k = 0; k < m_N; ++k)
            {
                Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + slack_off + k * nc);
                total += sk.sum();
            }
        }
        if(m_has_term_slack)
        {
            const int term_off = u_offset + m_N * nu + m_num_path_slack;
            Eigen::Map<const Vector<Scalar, NTC>> st(z.data() + term_off);
            total += st.sum();
        }
        m_last_diagnostics.total_slack = total;
    }

    void shift_warm_start(const Eigen::VectorX<Scalar>& sol)
    {
        m_warm_z = sol;
        shift_warm_states(sol);
        shift_warm_inputs(sol);
        shift_warm_path_slack(sol);
        zero_warm_terminal_slack();
    }

    void shift_warm_states(const Eigen::VectorX<Scalar>& sol)
    {
        for(int k = 0; k < m_N; ++k)
            m_warm_z.segment(k * nx, nx) = sol.segment((k + 1) * nx, nx);
    }

    void shift_warm_inputs(const Eigen::VectorX<Scalar>& sol)
    {
        const int u_offset = (m_N + 1) * nx;
        for(int k = 0; k < m_N - 1; ++k)
            m_warm_z.segment(u_offset + k * nu, nu) = sol.segment(u_offset + (k + 1) * nu, nu);
    }

    void shift_warm_path_slack(const Eigen::VectorX<Scalar>& sol)
    {
        if(!m_has_path_slack)
            return;
        const int u_offset = (m_N + 1) * nx;
        const int slack_off = u_offset + m_N * nu;
        for(int k = 0; k < m_N - 1; ++k)
            m_warm_z.segment(slack_off + k * nc, nc) = sol.segment(slack_off + (k + 1) * nc, nc);
        m_warm_z.segment(slack_off + (m_N - 1) * nc, nc).setZero();
    }

    void zero_warm_terminal_slack()
    {
        if(!m_has_term_slack)
            return;
        const int u_offset = (m_N + 1) * nx;
        const int term_off = u_offset + m_N * nu + m_num_path_slack;
        m_warm_z.segment(term_off, ntc).setZero();
    }

    auto extract_first_input(const Eigen::VectorX<Scalar>& sol) -> Vector<Scalar, NU>
    {
        const int u_offset = (m_N + 1) * nx;
        Vector<Scalar, NU> u0 = sol.segment(u_offset, nu);
        m_u_prev = u0;
        return u0;
    }

    Dynamics m_dynamics;
    nmpc_config<Scalar, NX, NU, NC, NTC> m_config;
    int m_N;
    bool m_has_path_slack;
    bool m_has_term_slack;
    int m_num_path_slack;
    int m_num_term_slack;
    int m_num_vars;

    std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> m_state;
    // Held behind a stable heap address so a defaulted move does not relocate the
    // object the solver's bridge points at by reference (argmin_problem.h:24-27).
    std::unique_ptr<nlp_problem<Scalar>> m_problem;
    Solver m_solver{};

    Eigen::VectorX<Scalar> m_warm_z;
    Eigen::VectorX<Scalar> m_last_solution;
    mpc_diagnostics<Scalar> m_last_diagnostics{};
    Vector<Scalar, NU> m_u_prev{Vector<Scalar, NU>::Zero()};
    bool m_setup_failed{false};
};

}

#endif
