#ifndef HPP_GUARD_CTRLPP_CONTROL_MRAC_H
#define HPP_GUARD_CTRLPP_CONTROL_MRAC_H

/// @brief Stateful MRAC controller with Lyapunov-based adaptation and compile-time robustification.
///
/// @cite slotine1991 -- Slotine & Li, "Applied Nonlinear Control", 1991, Ch. 8

#include "ctrlpp/types.h"

#include "ctrlpp/control/mrac_config.h"
#include "ctrlpp/control/mrac_policies.h"

#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"

#include <cmath>
#include <cstddef>
#include <type_traits>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX = 1, std::size_t NU = 1,
          typename Robustification = no_robustification>
class mrac_controller
{
public:
    using config_type = mrac_config<Scalar, NX, NU, Robustification>;
    using state_type = Vector<Scalar, NX>;
    using input_type = Vector<Scalar, NU>;

    explicit mrac_controller(const config_type& cfg)
        : m_cfg{cfg}
        , m_x_model{cfg.x_model_0}
        , m_theta_x{cfg.theta_x_0}
        , m_theta_r{cfg.theta_r_0}
    {
    }

    auto evaluate(const state_type& x, const input_type& r) -> Scalar
    {
        m_x_model = propagate(m_cfg.reference_model, m_x_model, r);

        m_tracking_error = x - m_x_model;

        Scalar e = m_tracking_error[0];

        Scalar u = m_theta_x * x[0] + m_theta_r * r[0];

        if constexpr(std::is_same_v<Robustification, no_robustification>)
        {
            m_theta_x -= m_cfg.gamma * m_cfg.sign_b * e * x[0];
            m_theta_r -= m_cfg.gamma * m_cfg.sign_b * e * r[0];
        }
        else if constexpr(std::is_same_v<Robustification, dead_zone>)
        {
            if(std::abs(e) > m_cfg.robustification.threshold)
            {
                m_theta_x -= m_cfg.gamma * m_cfg.sign_b * e * x[0];
                m_theta_r -= m_cfg.gamma * m_cfg.sign_b * e * r[0];
            }
        }
        else if constexpr(std::is_same_v<Robustification, sigma_modification>)
        {
            m_theta_x -= m_cfg.gamma * (m_cfg.sign_b * e * x[0] + m_cfg.robustification.sigma * m_theta_x);
            m_theta_r -= m_cfg.gamma * (m_cfg.sign_b * e * r[0] + m_cfg.robustification.sigma * m_theta_r);
        }
        else if constexpr(std::is_same_v<Robustification, e_modification>)
        {
            auto abs_e = std::abs(e);
            m_theta_x -= m_cfg.gamma * (m_cfg.sign_b * e * x[0] + m_cfg.robustification.delta * abs_e * m_theta_x);
            m_theta_r -= m_cfg.gamma * (m_cfg.sign_b * e * r[0] + m_cfg.robustification.delta * abs_e * m_theta_r);
        }

        return u;
    }

    auto theta_x() const -> Scalar { return m_theta_x; }

    auto theta_r() const -> Scalar { return m_theta_r; }

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
    config_type m_cfg;
    state_type m_x_model{};
    state_type m_tracking_error{};
    Scalar m_theta_x{};
    Scalar m_theta_r{};
};

}

#endif
