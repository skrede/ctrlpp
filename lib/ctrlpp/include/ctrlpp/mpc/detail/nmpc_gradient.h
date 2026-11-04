#ifndef HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_GRADIENT_H
#define HPP_GUARD_CTRLPP_MPC_DETAIL_NMPC_GRADIENT_H

#include "ctrlpp/mpc/detail/nmpc_layout.h"

#include "ctrlpp/types.h"

#include "ctrlpp/detail/numerical_diff.h"

#include "ctrlpp/mpc/nmpc_config.h"

#include <Eigen/Dense>

#include <span>
#include <memory>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <functional>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC, typename State>
class nmpc_quadratic_gradient
{
public:
    nmpc_quadratic_gradient(const nmpc_config<Scalar, NX, NU, NC, NTC>& config, std::shared_ptr<State> state, const nmpc_layout& layout)
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

    void operator()(std::span<const Scalar> z, std::span<Scalar> grad) const
    {
        const auto Qf = m_config.Qf.value_or(m_config.Q);
        set_stage_gradient(z, grad);
        set_terminal_gradient(Qf, z, grad);
        set_slack_gradient(grad);
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

    void set_stage_gradient(std::span<const Scalar> z, std::span<Scalar> grad) const
    {
        for(std::int32_t k = 0; k < m_N; ++k)
        {
            Eigen::Map<const Vector<Scalar, NX>> xk(z.data() + nmpc_layout::x_offset + k * nx);
            Eigen::Map<const Vector<Scalar, NU>> uk(z.data() + m_u_offset + k * nu);

            const auto& x_ref_k = (static_cast<std::size_t>(k) < m_state->x_ref.size()) ? m_state->x_ref[static_cast<std::size_t>(k)] : Vector<Scalar, NX>::Zero();

            Eigen::Map<Vector<Scalar, NX>> gx(grad.data() + nmpc_layout::x_offset + k * nx);
            gx = m_config.Q * (xk - x_ref_k);

            Eigen::Map<Vector<Scalar, NU>> gu(grad.data() + m_u_offset + k * nu);
            gu = m_config.R * uk;
        }
    }

    void set_terminal_gradient(const Matrix<Scalar, NX, NX>& Qf, std::span<const Scalar> z, std::span<Scalar> grad) const
    {
        Eigen::Map<const Vector<Scalar, NX>> xN(z.data() + nmpc_layout::x_offset + m_N * nx);
        const auto& x_ref_N = (static_cast<std::size_t>(m_N) < m_state->x_ref.size()) ? m_state->x_ref[static_cast<std::size_t>(m_N)] : Vector<Scalar, NX>::Zero();
        Eigen::Map<Vector<Scalar, NX>> gxN(grad.data() + nmpc_layout::x_offset + m_N * nx);
        gxN = Qf * (xN - x_ref_N);
    }

    void set_slack_gradient(std::span<Scalar> grad) const
    {
        if(m_has_path_slack)
        {
            for(std::int32_t k = 0; k < m_N; ++k)
            {
                Eigen::Map<Vector<Scalar, NC>> gs(grad.data() + m_path_slack_offset + k * nc);
                gs = m_config.path_penalty;
            }
        }
        if(m_has_term_slack)
        {
            Eigen::Map<Vector<Scalar, NTC>> gst(grad.data() + m_term_slack_offset);
            gst = m_config.terminal_penalty;
        }
    }
};

// The quadratic tracking cost has the analytic gradient H z + q, evaluated block by block without forming H;
// a custom stage or terminal cost falls back to a finite-difference gradient whose scratch is sized once.
template <typename Scalar, std::size_t NX, std::size_t NU, std::size_t NC, std::size_t NTC, typename State>
std::function<void(std::span<const Scalar>, std::span<Scalar>)> make_nmpc_gradient(
        const nmpc_config<Scalar, NX, NU, NC, NTC>& config, std::shared_ptr<State> state, const nmpc_layout& layout, const std::function<Scalar(std::span<const Scalar>)>& cost)
{
    if(!config.stage_cost_cb && !config.terminal_cost_cb)
    {
        return nmpc_quadratic_gradient<Scalar, NX, NU, NC, NTC, State>(config, std::move(state), layout);
    }
    return [cost, scratch = std::vector<Scalar>(static_cast<std::size_t>(layout.n_vars))](std::span<const Scalar> z, std::span<Scalar> grad) mutable
    { finite_diff_gradient<Scalar>(cost, z, grad, std::span<Scalar>{scratch.data(), scratch.size()}); };
}

}

#endif
