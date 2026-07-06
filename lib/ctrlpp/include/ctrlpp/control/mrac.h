#ifndef HPP_GUARD_CTRLPP_CONTROL_MRAC_H
#define HPP_GUARD_CTRLPP_CONTROL_MRAC_H

/// @brief Stateful MRAC controller with Lyapunov-based adaptation and compile-time robustification.
///
/// The adaptation law projects the tracking error through B^T only, which is the
/// Lyapunov gradient B^T P e specialized to P = I. This is exact when the reference
/// model is chosen so that A_m^T + A_m is negative definite (P = I solves the
/// Lyapunov equation A_m^T P + P A_m = -Q); for a general stable A_m a different P
/// would be required, and no configurable P weighting is provided here. The update
/// also carries no explicit dt factor, so the adaptation gains gamma_x and gamma_r
/// absorb the sample time: rescale them proportionally if the sample rate changes.
///
/// @cite slotine1991 -- Slotine & Li, "Applied Nonlinear Control", 1991, Ch. 8

#include "ctrlpp/types.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/control/mrac_config.h"
#include "ctrlpp/control/mrac_policies.h"

#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"

#include <cmath>
#include <cstddef>
#include <type_traits>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NX = 1, std::size_t NU = 1,
          typename Robustification = no_robustification>
class mrac_controller
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

public:
    using config_type = mrac_config<Scalar, NX, NU, Robustification>;
    using state_type = Vector<Scalar, NX>;
    using input_type = Vector<Scalar, NU>;
    using theta_x_type = Matrix<Scalar, NU, NX>;
    using theta_r_type = Matrix<Scalar, NU, NU>;

    explicit mrac_controller(const config_type& cfg)
        : m_cfg{cfg}
        , m_x_model{cfg.x_model_0}
        , m_tracking_error{state_type::Zero()}
        , m_theta_x{cfg.theta_x_0}
        , m_theta_r{cfg.theta_r_0}
    {
    }

    auto evaluate(const state_type& x, const input_type& r) -> input_type
    {
        m_x_model = propagate(m_cfg.reference_model, m_x_model, r);

        m_tracking_error = x - m_x_model;

        input_type u = m_theta_x * x + m_theta_r * r;

        auto e_proj = (m_cfg.reference_model.B.transpose() * m_tracking_error).eval();

        if constexpr(std::is_same_v<Robustification, no_robustification>)
        {
            m_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x;
            m_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r;
        }
        else if constexpr(std::is_same_v<Robustification, dead_zone>)
        {
            if(compute_error_norm(m_tracking_error) > m_cfg.robustification.threshold)
            {
                m_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x;
                m_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r;
            }
        }
        else if constexpr(std::is_same_v<Robustification, sigma_modification>)
        {
            m_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x
                                   + m_cfg.robustification.sigma * m_theta_x;
            m_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r
                                   + m_cfg.robustification.sigma * m_theta_r;
        }
        else if constexpr(std::is_same_v<Robustification, e_modification>)
        {
            auto e_norm = compute_error_norm(m_tracking_error);
            m_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x
                                   + m_cfg.robustification.delta * e_norm * m_theta_x;
            m_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r
                                   + m_cfg.robustification.delta * e_norm * m_theta_r;
        }

        return u;
    }

    auto theta_x() const -> const theta_x_type& { return m_theta_x; }

    auto theta_r() const -> const theta_r_type& { return m_theta_r; }

    auto tracking_error() const -> const state_type& { return m_tracking_error; }

    auto x_model() const -> const state_type& { return m_x_model; }

    void reset()
    {
        m_x_model = m_cfg.x_model_0;
        m_theta_x = m_cfg.theta_x_0;
        m_theta_r = m_cfg.theta_r_0;
        m_tracking_error = state_type::Zero();
    }

private:
    auto compute_error_norm(const state_type& e) const -> Scalar
    {
        return std::sqrt((e.transpose() * m_cfg.W * e)(0, 0));
    }

    config_type m_cfg;
    state_type m_x_model{};
    state_type m_tracking_error{};
    theta_x_type m_theta_x{};
    theta_r_type m_theta_r{};
};

}

#endif
