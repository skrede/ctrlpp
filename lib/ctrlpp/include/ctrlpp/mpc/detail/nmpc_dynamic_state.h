#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_DYNAMIC_STATE_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_DYNAMIC_STATE_H

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/mpc/nlp_formulation.h"

#include <Eigen/Dense>

#include <limits>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, typename Solver, typename Dynamics, std::size_t NC, std::size_t NTC>
class nmpc_dynamic_state
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int ntc = static_cast<int>(NTC);

protected:
    Dynamics m_dynamics;
    nmpc_config<Scalar, NX, NU, NC, NTC> m_config;
    std::int32_t m_N;
    bool m_has_path_slack;
    bool m_has_term_slack;
    Eigen::Index m_num_path_slack;
    Eigen::Index m_num_term_slack;
    Eigen::Index m_num_vars;

    std::shared_ptr<nmpc_formulation_state<Scalar, NX, NU>> m_state;
    std::unique_ptr<nlp_problem<Scalar>> m_problem;
    Solver m_solver;

    Eigen::VectorX<Scalar> m_warm_z;
    Eigen::VectorX<Scalar> m_last_solution;
    mpc_diagnostics<Scalar> m_last_diagnostics;
    Vector<Scalar, NU> m_u_prev;
    bool m_setup_failed;
    bool m_has_solution;

    nmpc_dynamic_state(Dynamics dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, Solver&& solver)
        : m_dynamics{std::move(dynamics)}
        , m_config{config}
        , m_N{config.horizon}
        , m_has_path_slack{config.soft_constraints && NC > 0 && config.path_constraint_cb.has_value()}
        , m_has_term_slack{config.soft_constraints && NTC > 0 && config.terminal_constraint_cb.has_value()}
        , m_num_path_slack{m_has_path_slack ? m_N * nc : 0}
        , m_num_term_slack{m_has_term_slack ? ntc : 0}
        , m_num_vars{(m_N + 1) * nx + m_N * nu + m_num_path_slack + m_num_term_slack}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>()}
        , m_solver{std::move(solver)}
        , m_last_diagnostics{}
        , m_u_prev{Vector<Scalar, NU>::Zero()}
        , m_setup_failed{false}
        , m_has_solution{false}
    {
        m_state->x_ref.resize(static_cast<std::size_t>(m_N + 1), Vector<Scalar, NX>::Zero());
        m_problem = std::make_unique<nlp_problem<Scalar>>(detail::build_nmpc_problem<Scalar, NX, NU, NC, NTC>(m_dynamics, m_config, m_state));
        m_warm_z = Eigen::VectorX<Scalar>::Zero(m_num_vars);
    }

    // A defaulted move is correct: the problem sits behind a unique_ptr, so the reference the solver's bridge holds
    // to it survives the move.
    nmpc_dynamic_state(nmpc_dynamic_state&&) = default;
    nmpc_dynamic_state& operator=(nmpc_dynamic_state&&) = default;
    ~nmpc_dynamic_state() = default;

    // A copy is an independent fork: it deep-copies the formulation state, poses its own problem against that copy
    // and sets its forked solver up against it, so the two controllers share no mutable state.
    nmpc_dynamic_state(const nmpc_dynamic_state& other)
        : m_dynamics{other.m_dynamics}
        , m_config{other.m_config}
        , m_N{other.m_N}
        , m_has_path_slack{other.m_has_path_slack}
        , m_has_term_slack{other.m_has_term_slack}
        , m_num_path_slack{other.m_num_path_slack}
        , m_num_term_slack{other.m_num_term_slack}
        , m_num_vars{other.m_num_vars}
        , m_state{std::make_shared<nmpc_formulation_state<Scalar, NX, NU>>(*other.m_state)}
        , m_problem{std::make_unique<nlp_problem<Scalar>>(detail::build_nmpc_problem<Scalar, NX, NU, NC, NTC>(m_dynamics, m_config, m_state))}
        , m_solver{other.m_solver}
        , m_warm_z{other.m_warm_z}
        , m_last_solution{other.m_last_solution}
        , m_last_diagnostics{other.m_last_diagnostics}
        , m_u_prev{other.m_u_prev}
        , m_setup_failed{false}
        , m_has_solution{other.m_has_solution}
    {
        m_setup_failed = !detail::setup_nlp_solver(m_solver, *m_problem).has_value();
    }

    nmpc_dynamic_state& operator=(const nmpc_dynamic_state& other)
    {
        nmpc_dynamic_state tmp{other};
        *this = std::move(tmp);
        return *this;
    }

    // Both derived dimensions, N (nx + nu + nc) + nx + ntc decision variables and N (nx + 2 nu + nc) + nx + ntc
    // constraint rows, stay representable in the horizon's own type up to this horizon.
    static std::int32_t horizon_bound()
    {
        constexpr std::int32_t per_step = nx + 2 * nu + nc;
        constexpr std::int32_t constant_dimensions = nx + ntc;
        return (std::numeric_limits<std::int32_t>::max() - constant_dimensions) / per_step;
    }
};

}

#endif
