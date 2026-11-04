#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_COST_TERMS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_COST_TERMS_H

#include "ctrlpp/mpc/detail/nmpc_layout.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/nmpc_config.h"

#include <Eigen/Dense>

#include <span>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC, typename State>
class nmpc_cost
{
public:
    nmpc_cost(const nmpc_config<Scalar, NX, NU, NC, NTC>& config, std::shared_ptr<State> state, const nmpc_layout& layout)
        : m_N{layout.N}
        , m_has_path_slack{layout.has_path_slack}
        , m_has_term_slack{layout.has_term_slack}
        , m_u_offset{layout.u_offset}
        , m_state{std::move(state)}
        , m_path_slack_offset{layout.path_slack_offset}
        , m_term_slack_offset{layout.term_slack_offset}
        , m_config{config}
    {
    }

    Scalar operator()(std::span<const Scalar> z) const
    {
        Scalar total = stage_costs(z);
        add_terminal_cost(total, z);
        add_slack_penalty(total, z);
        return total;
    }

private:
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int nc = static_cast<int>(NC);

    std::int32_t m_N;
    bool m_has_path_slack;
    bool m_has_term_slack;
    Eigen::Index m_u_offset;
    std::shared_ptr<State> m_state;
    Eigen::Index m_path_slack_offset;
    Eigen::Index m_term_slack_offset;
    // const, as a lambda's by-copy capture would be: moving the closure into std::function then copies it, which
    // measured smaller construction frames than moving it.
    const nmpc_config<Scalar, NX, NU, NC, NTC> m_config;

    Scalar stage_costs(std::span<const Scalar> z) const
    {
        Scalar total{0};
        for(std::int32_t k = 0; k < m_N; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + nmpc_layout::x_offset + k * nx);
            Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + m_u_offset + k * nu);

            if(m_config.stage_cost_cb)
            {
                total += (*m_config.stage_cost_cb)(xk, uk);
            }
            else
            {
                const auto& x_ref_k = (static_cast<std::size_t>(k) < m_state->x_ref.size()) ? m_state->x_ref[static_cast<std::size_t>(k)] : Vector<Scalar, NX>::Zero();
                Vector<Scalar, NX> dx = xk - x_ref_k;
                total += Scalar{0.5} * dx.dot(m_config.Q * dx) + Scalar{0.5} * uk.dot(m_config.R * uk);
            }
        }
        return total;
    }

    void add_terminal_cost(Scalar& total, std::span<const Scalar> z) const
    {
        Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + nmpc_layout::x_offset + m_N * nx);
        if(m_config.terminal_cost_cb)
        {
            total += (*m_config.terminal_cost_cb)(xN);
        }
        else
        {
            const auto Qf = m_config.Qf.value_or(m_config.Q);
            const auto& x_ref_N = (static_cast<std::size_t>(m_N) < m_state->x_ref.size()) ? m_state->x_ref[static_cast<std::size_t>(m_N)] : Vector<Scalar, NX>::Zero();
            Vector<Scalar, NX> dx = xN - x_ref_N;
            total += Scalar{0.5} * dx.dot(Qf * dx);
        }
    }

    void add_slack_penalty(Scalar& total, std::span<const Scalar> z) const
    {
        if(m_has_path_slack)
        {
            for(std::int32_t k = 0; k < m_N; ++k)
            {
                Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + m_path_slack_offset + k * nc);
                total += m_config.path_penalty.dot(sk);
            }
        }
        if(m_has_term_slack)
        {
            Eigen::Map<const Vector<Scalar, NTC>> st(z.data() + m_term_slack_offset);
            total += m_config.terminal_penalty.dot(st);
        }
    }
};

}

#endif
