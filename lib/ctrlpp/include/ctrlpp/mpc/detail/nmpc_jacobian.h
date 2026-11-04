#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_JACOBIAN_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_JACOBIAN_H

#include "ctrlpp/mpc/detail/nmpc_layout.h"

#include "ctrlpp/mpc/nmpc_config.h"

#include <Eigen/Dense>

#include <span>
#include <cmath>
#include <limits>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <functional>

namespace ctrlpp::detail
{

// The dynamics and constraint callables expose no Jacobian of their own, so the whole Jacobian is central-differenced
// (the eps^(1/3) magnitude-scaled step of finite_diff_*) and its structural entries (identity rows, rate +/-1, slack -1)
// are then overwritten exactly. argmin reads the Jacobian column-major: entry (i, j) is at i + j * n_constraints.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC>
class nmpc_constraint_jacobian
{
public:
    nmpc_constraint_jacobian(
            const nmpc_config<Scalar, NX, NU, NC, NTC>& config, const nmpc_layout& layout, const std::function<void(std::span<const Scalar>, std::span<Scalar>)>& constraints)
        : m_N{layout.N}
        , m_n_vars{layout.n_vars}
        , m_has_path_slack{layout.has_path_slack}
        , m_has_term_slack{layout.has_term_slack}
        , m_u_offset{layout.u_offset}
        , m_rate_start{layout.rate_start}
        , m_c_plus(static_cast<std::size_t>(layout.n_constraints))
        , m_n_constraints{layout.n_constraints}
        , m_c_minus(static_cast<std::size_t>(layout.n_constraints))
        , m_path_con_start{layout.path_con_start}
        , m_term_con_start{layout.term_con_start}
        , m_z_scratch(static_cast<std::size_t>(layout.n_vars))
        , m_path_slack_offset{layout.path_slack_offset}
        , m_term_slack_offset{layout.term_slack_offset}
        , m_config{config}
        , m_constraints_cb{constraints}
    {
    }

    void operator()(std::span<const Scalar> z, std::span<Scalar> jac)
    {
        central_difference(z, jac);
        set_dynamics_rows(jac);
        if(m_config.du_max)
        {
            set_rate_rows(jac);
        }
        if(m_config.path_constraint_cb && m_has_path_slack)
        {
            set_path_slack_columns(jac);
        }
        if(m_config.terminal_constraint_cb && m_has_term_slack)
        {
            set_terminal_slack_columns(jac);
        }
    }

private:
    static constexpr int nx = static_cast<int>(NX);
    static constexpr int nu = static_cast<int>(NU);
    static constexpr int nc = static_cast<int>(NC);
    static constexpr int ntc = static_cast<int>(NTC);

    std::int32_t m_N;
    Eigen::Index m_n_vars;
    bool m_has_path_slack;
    bool m_has_term_slack;
    Eigen::Index m_u_offset;
    Eigen::Index m_rate_start;
    std::vector<Scalar> m_c_plus;
    Eigen::Index m_n_constraints;
    std::vector<Scalar> m_c_minus;
    Eigen::Index m_path_con_start;
    Eigen::Index m_term_con_start;
    std::vector<Scalar> m_z_scratch;
    Eigen::Index m_path_slack_offset;
    Eigen::Index m_term_slack_offset;
    // const, as a lambda's by-copy capture would be: moving the closure into std::function then copies it, which
    // measured smaller construction frames than moving it.
    const nmpc_config<Scalar, NX, NU, NC, NTC> m_config;
    std::function<void(std::span<const Scalar>, std::span<Scalar>)> m_constraints_cb;

    Scalar& at(std::span<Scalar> J, Eigen::Index i, Eigen::Index j) const
    {
        return J[static_cast<std::size_t>(i) + static_cast<std::size_t>(j) * static_cast<std::size_t>(m_n_constraints)];
    }

    void central_difference(std::span<const Scalar> z, std::span<Scalar> jac)
    {
        const auto step_scale = std::cbrt(std::numeric_limits<Scalar>::epsilon());
        std::copy(z.begin(), z.end(), m_z_scratch.begin());
        for(Eigen::Index j = 0; j < m_n_vars; ++j)
        {
            const auto jz = static_cast<std::size_t>(j);
            const Scalar h_raw = step_scale * std::max(Scalar{1}, std::abs(z[jz]));
            const Scalar temp = z[jz] + h_raw;
            const Scalar h = temp - z[jz];
            const Scalar orig = m_z_scratch[jz];

            m_z_scratch[jz] = orig + h;
            m_constraints_cb(std::span<const Scalar>{m_z_scratch.data(), m_z_scratch.size()}, std::span<Scalar>{m_c_plus.data(), m_c_plus.size()});
            m_z_scratch[jz] = orig - h;
            m_constraints_cb(std::span<const Scalar>{m_z_scratch.data(), m_z_scratch.size()}, std::span<Scalar>{m_c_minus.data(), m_c_minus.size()});

            for(Eigen::Index i = 0; i < m_n_constraints; ++i)
            {
                at(jac, i, j) = (m_c_plus[static_cast<std::size_t>(i)] - m_c_minus[static_cast<std::size_t>(i)]) / (Scalar{2} * h);
            }
            m_z_scratch[jz] = orig;
        }
    }

    void set_dynamics_rows(std::span<Scalar> jac) const
    {
        for(Eigen::Index i = 0; i < nx; ++i)
        {
            for(Eigen::Index j = 0; j < m_n_vars; ++j)
            {
                at(jac, nmpc_layout::eq_start + i, j) = Scalar{0};
            }
            at(jac, nmpc_layout::eq_start + i, nmpc_layout::x_offset + i) = Scalar{1};
        }

        for(std::int32_t k = 0; k < m_N; ++k)
        {
            for(Eigen::Index i = 0; i < nx; ++i)
            {
                at(jac, nmpc_layout::eq_start + (k + 1) * nx + i, nmpc_layout::x_offset + (k + 1) * nx + i) = Scalar{1};
            }
        }
    }

    void set_rate_rows(std::span<Scalar> jac) const
    {
        for(std::int32_t k = 0; k < m_N; ++k)
        {
            for(Eigen::Index j = 0; j < nu; ++j)
            {
                const Eigen::Index row_p = m_rate_start + k * nu * 2 + j * 2;
                const Eigen::Index row_m = row_p + 1;
                for(Eigen::Index col = 0; col < m_n_vars; ++col)
                {
                    at(jac, row_p, col) = Scalar{0};
                    at(jac, row_m, col) = Scalar{0};
                }

                at(jac, row_p, m_u_offset + k * nu + j) = Scalar{1};
                at(jac, row_m, m_u_offset + k * nu + j) = Scalar{-1};
                if(k > 0)
                {
                    at(jac, row_p, m_u_offset + (k - 1) * nu + j) = Scalar{-1};
                    at(jac, row_m, m_u_offset + (k - 1) * nu + j) = Scalar{1};
                }
            }
        }
    }

    void set_path_slack_columns(std::span<Scalar> jac) const
    {
        for(std::int32_t k = 0; k < m_N; ++k)
        {
            for(Eigen::Index i = 0; i < nc; ++i)
            {
                at(jac, m_path_con_start + k * nc + i, m_path_slack_offset + k * nc + i) = Scalar{-1};
            }
        }
    }

    void set_terminal_slack_columns(std::span<Scalar> jac) const
    {
        for(Eigen::Index i = 0; i < ntc; ++i)
        {
            at(jac, m_term_con_start + i, m_term_slack_offset + i) = Scalar{-1};
        }
    }
};

}

#endif
