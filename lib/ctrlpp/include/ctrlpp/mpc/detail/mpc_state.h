#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_MPC_STATE_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_MPC_STATE_H

#include "ctrlpp/mpc/detail/qp_layout.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/qp_types.h"
#include "ctrlpp/mpc/diagnostics.h"
#include "ctrlpp/mpc/terminal_cost.h"
#include "ctrlpp/mpc/qp_formulation.h"

#include "ctrlpp/model/state_space.h"

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>
#include <utility>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, typename Solver, std::size_t NY, typename Config>
class mpc_state
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);

protected:
    Solver m_solver;
    Config m_config;
    discrete_state_space<Scalar, NX, NU, NY> m_system;
    Matrix<Scalar, NX, NX> m_Q_state;
    Matrix<Scalar, NX, NX> m_Qf_state;
    Matrix<Scalar, NX, NY> m_CtQ;
    Matrix<Scalar, NX, NY> m_Qf_linear;

    qp_update<Scalar> m_update;
    Eigen::VectorX<Scalar> m_last_primal;
    Eigen::VectorX<Scalar> m_last_dual;
    mpc_diagnostics<Scalar> m_last_diagnostics;
    Vector<Scalar, NU> m_u_prev;
    bool m_has_solution;

    Eigen::Index m_n_dec;
    Eigen::Index m_n_con;
    Eigen::Index m_n_x_total;
    Eigen::Index m_rate_row;

    mpc_state(const discrete_state_space<Scalar, NX, NU, NY>& system, const Config& config, Matrix<Scalar, NX, NX>&& Q_state, mpc_terminal_weights<Scalar, NX, NY>&& terminal,
            Solver&& solver)
        : m_solver{std::move(solver)}
        , m_config{config}
        , m_system{system}
        , m_Q_state{std::move(Q_state)}
        , m_Qf_state{std::move(terminal.state)}
        , m_CtQ{system.C.transpose() * config.Q}
        , m_Qf_linear{std::move(terminal.linear)}
        , m_update{}
        , m_last_primal{}
        , m_last_dual{}
        , m_last_diagnostics{}
        , m_u_prev{Vector<Scalar, NU>::Zero()}
        , m_has_solution{false}
        , m_n_dec{}
        , m_n_con{}
        , m_n_x_total{}
        , m_rate_row{}
    {
        compute_dimensions();
        allocate_update_vectors();
    }

    // One derivation sizes the constraint matrix, the bounds vectors and the rate row, so the artifacts handed to
    // the backend have the height of the problem this controller believes it posed.
    void compute_dimensions()
    {
        std::int32_t N = m_config.horizon;
        m_n_x_total = (N + 1) * nx;
        m_n_dec = m_n_x_total + N * nu + detail::slack_columns<NX>(N, uses_soft_state_constraints());

        auto const layout = constraint_layout();
        m_n_con = layout.n_con;
        m_rate_row = layout.rate_row;
    }

    // Present state bounds are necessary but not sufficient: enforcing them hard is the caller's policy choice.
    bool uses_soft_state_constraints() const { return (m_config.x_min.has_value() || m_config.x_max.has_value()) && !m_config.hard_state_constraints; }

    detail::qp_constraint_layout constraint_layout() const
    {
        return detail::constraint_layout<Scalar, NX, NU>(
                m_config.horizon, m_config.x_min, m_config.x_max, m_config.u_min, m_config.u_max, m_config.du_max, m_config.terminal_constraint_set);
    }

    qp_problem<Scalar> initial_problem() const
    {
        std::int32_t N = m_config.horizon;
        bool use_soft = uses_soft_state_constraints();

        // Both builders read the same optionals, so the matrix and the bounds
        // vectors cannot come out of this function with different row counts.
        auto P = detail::build_cost_matrix<Scalar, NX, NU>(N, m_Q_state, m_config.R, m_Qf_state, use_soft, m_config.soft_penalty, m_config.soft_state_penalty);
        auto A = detail::build_constraint_matrix<Scalar, NX, NU>(
                N, m_system.A, m_system.B, m_config.x_min, m_config.x_max, use_soft, m_config.u_min, m_config.u_max, m_config.du_max, m_config.terminal_constraint_set);

        Vector<Scalar, NX> x0_dummy = Vector<Scalar, NX>::Zero();
        auto [l, u] = detail::build_bounds_vectors<Scalar, NX, NU>(
                N, x0_dummy, m_config.x_min, m_config.x_max, m_config.u_min, m_config.u_max, m_config.du_max, m_config.terminal_constraint_set);

        auto q = detail::build_cost_vector<Scalar, NX, NU>(N, m_n_dec, m_Q_state, m_Qf_state);

        return qp_problem<Scalar>{.P = std::move(P), .q = std::move(q), .A = std::move(A), .l = std::move(l), .u = std::move(u)};
    }

    void allocate_update_vectors()
    {
        m_update.q.resize(m_n_dec);
        m_update.q.setZero();
        m_update.l.resize(m_n_con);
        m_update.l.setZero();
        m_update.u.resize(m_n_con);
        m_update.u.setZero();
    }
};

}

#endif
