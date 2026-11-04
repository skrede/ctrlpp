#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_DYNAMIC_RESULT_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_DYNAMIC_RESULT_H

#include "ctrlpp/mpc/detail/nmpc_dynamic_state.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/nlp_solver.h"
#include "ctrlpp/mpc/diagnostics.h"

#include <Eigen/Dense>

#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, typename Solver, typename Dynamics, std::size_t NC, std::size_t NTC>
class nmpc_dynamic_result : protected nmpc_dynamic_state<Scalar, NX, NU, Solver, Dynamics, NC, NTC>
{
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int ntc = static_cast<int>(NTC);

    using base = nmpc_dynamic_state<Scalar, NX, NU, Solver, Dynamics, NC, NTC>;

protected:
    using base::base;

    void populate_diagnostics(const nlp_result<Scalar>& result)
    {
        mpc_diagnostics<Scalar> d;
        d.status = result.status;
        d.iterations = result.iterations;
        d.solve_time = result.solve_time;
        d.cost = result.objective;
        d.primal_residual = result.primal_residual;
        d.dual_residual = Scalar{0};
        d.max_constraint_violation = result.primal_residual;
        d.max_path_constraint_violation = Scalar{0};
        d.max_terminal_constraint_violation = Scalar{0};
        d.total_slack = Scalar{0};
        d.stop_criterion = result.stop_criterion;
        this->m_last_diagnostics = d;
    }

    void populate_constraint_diagnostics(const Eigen::VectorX<Scalar>& z)
    {
        const Eigen::Index u_offset = (this->m_N + 1) * nx;
        compute_path_constraint_violation(z, u_offset);
        compute_terminal_constraint_violation(z);
        compute_total_slack(z, u_offset);
    }

    void compute_path_constraint_violation(const Eigen::VectorX<Scalar>& z, Eigen::Index u_offset)
    {
        if constexpr(NC > 0)
        {
            if(this->m_config.path_constraint_cb)
            {
                const auto& g = *this->m_config.path_constraint_cb;
                Scalar max_viol{0};

                for(std::int32_t k = 0; k < this->m_N; ++k)
                {
                    Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + k * nx);
                    Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + u_offset + k * nu);
                    Vector<Scalar, NC> gk = g(xk, uk);

                    for(Eigen::Index i = 0; i < nc; ++i)
                        max_viol = std::max(max_viol, gk[i]);
                }

                this->m_last_diagnostics.max_path_constraint_violation = max_viol;
            }
        }
    }

    void compute_terminal_constraint_violation(const Eigen::VectorX<Scalar>& z)
    {
        if constexpr(NTC > 0)
        {
            if(this->m_config.terminal_constraint_cb)
            {
                const auto& h = *this->m_config.terminal_constraint_cb;
                Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + this->m_N * nx);
                Vector<Scalar, NTC> hN = h(xN);

                Scalar max_viol{0};
                for(Eigen::Index i = 0; i < ntc; ++i)
                    max_viol = std::max(max_viol, hN[i]);

                this->m_last_diagnostics.max_terminal_constraint_violation = max_viol;
            }
        }
    }

    void compute_total_slack(const Eigen::VectorX<Scalar>& z, Eigen::Index u_offset)
    {
        Scalar total{0};
        if(this->m_has_path_slack)
        {
            const Eigen::Index slack_off = u_offset + this->m_N * nu;
            for(std::int32_t k = 0; k < this->m_N; ++k)
            {
                Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + slack_off + k * nc);
                total += sk.sum();
            }
        }
        if(this->m_has_term_slack)
        {
            const Eigen::Index term_off = u_offset + this->m_N * nu + this->m_num_path_slack;
            Eigen::Map<const Vector<Scalar, NTC>> st(z.data() + term_off);
            total += st.sum();
        }
        this->m_last_diagnostics.total_slack = total;
    }

    void shift_warm_start(const Eigen::VectorX<Scalar>& sol)
    {
        this->m_warm_z = sol;
        shift_warm_states(sol);
        shift_warm_inputs(sol);
        shift_warm_path_slack(sol);
        zero_warm_terminal_slack();
    }

    void shift_warm_states(const Eigen::VectorX<Scalar>& sol)
    {
        for(std::int32_t k = 0; k < this->m_N; ++k)
            this->m_warm_z.segment(k * nx, nx) = sol.segment((k + 1) * nx, nx);
    }

    void shift_warm_inputs(const Eigen::VectorX<Scalar>& sol)
    {
        const Eigen::Index u_offset = (this->m_N + 1) * nx;
        for(std::int32_t k = 0; k < this->m_N - 1; ++k)
            this->m_warm_z.segment(u_offset + k * nu, nu) = sol.segment(u_offset + (k + 1) * nu, nu);
    }

    void shift_warm_path_slack(const Eigen::VectorX<Scalar>& sol)
    {
        if(!this->m_has_path_slack)
            return;
        const Eigen::Index u_offset = (this->m_N + 1) * nx;
        const Eigen::Index slack_off = u_offset + this->m_N * nu;
        for(std::int32_t k = 0; k < this->m_N - 1; ++k)
            this->m_warm_z.segment(slack_off + k * nc, nc) = sol.segment(slack_off + (k + 1) * nc, nc);
        this->m_warm_z.segment(slack_off + (this->m_N - 1) * nc, nc).setZero();
    }

    void zero_warm_terminal_slack()
    {
        if(!this->m_has_term_slack)
            return;
        const Eigen::Index u_offset = (this->m_N + 1) * nx;
        const Eigen::Index term_off = u_offset + this->m_N * nu + this->m_num_path_slack;
        this->m_warm_z.segment(term_off, ntc).setZero();
    }

    Vector<Scalar, NU> extract_first_input(const Eigen::VectorX<Scalar>& sol)
    {
        const Eigen::Index u_offset = (this->m_N + 1) * nx;
        Vector<Scalar, NU> u0 = sol.segment(u_offset, nu);
        this->m_u_prev = u0;
        return u0;
    }
};

}

#endif
