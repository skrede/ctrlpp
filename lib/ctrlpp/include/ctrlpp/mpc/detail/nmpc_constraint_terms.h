#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_CONSTRAINT_TERMS_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_CONSTRAINT_TERMS_H

#include "ctrlpp/mpc/detail/nmpc_layout.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/nmpc_config.h"

#include <Eigen/Dense>

#include <span>
#include <memory>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC, typename Dynamics, typename State>
class nmpc_constraints
{
    // A plain function is held by reference and any other callable as a const copy, as a lambda's by-copy capture
    // holds them.
    using dynamics_holder = std::conditional_t<std::is_function_v<Dynamics>, Dynamics&, const Dynamics>;

public:
    nmpc_constraints(const Dynamics& dynamics, const nmpc_config<Scalar, NX, NU, NC, NTC>& config, std::shared_ptr<State> state, const nmpc_layout& layout)
        : m_N{layout.N}
        , m_has_path_slack{layout.has_path_slack}
        , m_has_term_slack{layout.has_term_slack}
        , m_u_offset{layout.u_offset}
        , m_rate_start{layout.rate_start}
        , m_dynamics{dynamics}
        , m_path_con_start{layout.path_con_start}
        , m_term_con_start{layout.term_con_start}
        , m_state{std::move(state)}
        , m_path_slack_offset{layout.path_slack_offset}
        , m_term_slack_offset{layout.term_slack_offset}
        , m_config{config}
    {
    }

    void operator()(std::span<const Scalar> z, std::span<Scalar> c) const
    {
        set_dynamics_residuals(z, c);
        if(m_config.du_max)
        {
            set_rate_residuals(z, c);
        }
        if(m_config.path_constraint_cb)
        {
            const auto& g = *m_config.path_constraint_cb;
            for(std::int32_t k = 0; k < m_N; ++k)
            {
                set_path_stage(g, z, c, k);
            }
        }
        if(m_config.terminal_constraint_cb)
        {
            set_terminal_residuals(z, c);
        }
    }

private:
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int ntc = static_cast<int>(NTC);

    std::int32_t m_N;
    bool m_has_path_slack;
    bool m_has_term_slack;
    Eigen::Index m_u_offset;
    Eigen::Index m_rate_start;
    dynamics_holder m_dynamics;
    Eigen::Index m_path_con_start;
    Eigen::Index m_term_con_start;
    std::shared_ptr<State> m_state;
    Eigen::Index m_path_slack_offset;
    Eigen::Index m_term_slack_offset;
    // const, as a lambda's by-copy capture would be: moving the closure into std::function then copies it, which
    // measured smaller construction frames than moving it.
    const nmpc_config<Scalar, NX, NU, NC, NTC> m_config;

    void set_dynamics_residuals(std::span<const Scalar> z, std::span<Scalar> c) const
    {
        for(Eigen::Index i = 0; i < nx; ++i)
        {
            c[static_cast<std::size_t>(nmpc_layout::eq_start + i)] = z[static_cast<std::size_t>(nmpc_layout::x_offset + i)] - m_state->x0[i];
        }

        for(std::int32_t k = 0; k < m_N; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + nmpc_layout::x_offset + k * nx);
            Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + m_u_offset + k * nu);

            Vector<Scalar, NX> x_next = m_dynamics(xk, uk);

            for(Eigen::Index i = 0; i < nx; ++i)
            {
                c[static_cast<std::size_t>(nmpc_layout::eq_start + (k + 1) * nx + i)] = z[static_cast<std::size_t>(nmpc_layout::x_offset + (k + 1) * nx + i)] - x_next[i];
            }
        }
    }

    void set_rate_residuals(std::span<const Scalar> z, std::span<Scalar> c) const
    {
        const auto& du_max = *m_config.du_max;
        for(std::int32_t k = 0; k < m_N; ++k)
        {
            Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + m_u_offset + k * nu);

            Vector<Scalar, NU> uk_prev;
            if(k == 0)
            {
                uk_prev = m_state->u_prev;
            }
            else
            {
                uk_prev = Eigen::Map<const Vector<Scalar, NU>>(z.data() + m_u_offset + (k - 1) * nu);
            }

            for(Eigen::Index j = 0; j < nu; ++j)
            {
                Scalar du = uk[j] - uk_prev[j];
                c[static_cast<std::size_t>(m_rate_start + k * nu * 2 + j * 2)] = du - du_max[j];
                c[static_cast<std::size_t>(m_rate_start + k * nu * 2 + j * 2 + 1)] = -du - du_max[j];
            }
        }
    }

    template <typename G>
    void set_path_stage(const G& g, std::span<const Scalar> z, std::span<Scalar> c, std::int32_t k) const
    {
        Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + nmpc_layout::x_offset + k * nx);
        Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + m_u_offset + k * nu);

        Vector<Scalar, NC> gk = g(xk, uk);

        if(m_has_path_slack)
        {
            Eigen::Map<const Vector<Scalar, NC>> sk(z.data() + m_path_slack_offset + k * nc);
            gk -= sk;
        }

        for(Eigen::Index i = 0; i < nc; ++i)
        {
            c[static_cast<std::size_t>(m_path_con_start + k * nc + i)] = gk[i];
        }
    }

    void set_terminal_residuals(std::span<const Scalar> z, std::span<Scalar> c) const
    {
        const auto& h = *m_config.terminal_constraint_cb;
        Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + nmpc_layout::x_offset + m_N * nx);

        Vector<Scalar, NTC> hN = h(xN);

        if(m_has_term_slack)
        {
            Eigen::Map<const Vector<Scalar, NTC>> st(z.data() + m_term_slack_offset);
            hN -= st;
        }

        for(Eigen::Index i = 0; i < ntc; ++i)
        {
            c[static_cast<std::size_t>(m_term_con_start + i)] = hN[i];
        }
    }
};

}

#endif
