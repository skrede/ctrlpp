#ifndef HPP_GUARD_CTRLPP_CONTROL_MRAC_H
#define HPP_GUARD_CTRLPP_CONTROL_MRAC_H

// @cite slotine1991 -- Slotine & Li, "Applied Nonlinear Control", 1991, Ch. 8

#include "ctrlpp/control/mrac_config.h"
#include "ctrlpp/control/mrac_policies.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"

#include <cmath>
#include <cstddef>
#include <type_traits>

namespace ctrlpp
{

// Listed in the order check_step tests them, which follows the adaptation's data flow, so the cause named is the
// most upstream one.
enum class mrac_step_error
{
    non_finite_reference_state,
    non_finite_parameters,
    non_finite_state,
    non_finite_reference,
    non_finite_result,
};

enum class mrac_health
{
    ok,
    non_finite_carried_state,
};

template <ctrlpp_floating_scalar Scalar, std::size_t NX = 1, std::size_t NU = 1, typename Robustification = no_robustification>
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
        , m_health{mrac_health::ok}
        , m_theta_x{cfg.theta_x_0}
        , m_theta_r{cfg.theta_r_0}
        , m_tracking_error{state_type::Zero()}
    {
    }

    // The cycle is classified before any member is written, so a rejected cycle leaves every member bitwise
    // unchanged. The adaptation only accumulates into the parameters, so one admitted non-finite sample would be
    // permanent. The gradient is B^T P e specialized to P = I, a valid Lyapunov weight for the discrete-time reference
    // model only when I - A_m^T A_m is positive definite; it carries no dt factor, so gamma_x and gamma_r absorb the
    // sample time.
    expected<input_type, mrac_step_error> evaluate(const state_type& x, const input_type& r)
    {
        if(const auto step = check_step(x, r); !step)
            return unexpected(latch_health(step.error()));

        auto next_x_model = propagate(m_cfg.reference_model, m_x_model, r);
        auto next_tracking_error = (x - next_x_model).eval();
        auto u = (m_theta_x * x + m_theta_r * r).eval();
        auto e_proj = (m_cfg.reference_model.B.transpose() * next_tracking_error).eval();
        auto next_theta_x = m_theta_x;
        auto next_theta_r = m_theta_r;
        if constexpr(std::is_same_v<Robustification, sigma_modification>)
        {
            next_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x + m_cfg.robustification.sigma * m_theta_x;
            next_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r + m_cfg.robustification.sigma * m_theta_r;
        }
        else if constexpr(std::is_same_v<Robustification, e_modification>)
        {
            auto e_norm = compute_error_norm(next_tracking_error);
            next_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x + m_cfg.robustification.delta * e_norm * m_theta_x;
            next_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r + m_cfg.robustification.delta * e_norm * m_theta_r;
        }
        else if(gradient_enabled(next_tracking_error))
        {
            next_theta_x.noalias() -= m_cfg.sign_b * e_proj * x.transpose() * m_cfg.gamma_x;
            next_theta_r.noalias() -= m_cfg.sign_b * e_proj * r.transpose() * m_cfg.gamma_r;
        }
        return try_commit(std::move(next_x_model), std::move(next_tracking_error), u, e_proj, std::move(next_theta_x), std::move(next_theta_r));
    }

    const theta_x_type& theta_x() const { return m_theta_x; }

    const theta_r_type& theta_r() const { return m_theta_r; }

    const state_type& tracking_error() const { return m_tracking_error; }

    const state_type& x_model() const { return m_x_model; }

    // Latches until reset; a rejected cycle mutates nothing and does not set it.
    mrac_health health() const { return m_health; }

    void reset()
    {
        m_x_model = m_cfg.x_model_0;
        m_theta_x = m_cfg.theta_x_0;
        m_theta_r = m_cfg.theta_r_0;
        m_tracking_error = state_type::Zero();
        m_health = mrac_health::ok;
    }

private:
    config_type m_cfg;
    state_type m_x_model;
    mrac_health m_health;
    theta_x_type m_theta_x;
    theta_r_type m_theta_r;
    state_type m_tracking_error;

    expected<void, mrac_step_error> check_step(const state_type& x, const input_type& r) const
    {
        if(!m_x_model.allFinite())
            return unexpected(mrac_step_error::non_finite_reference_state);
        if(!m_theta_x.allFinite() || !m_theta_r.allFinite())
            return unexpected(mrac_step_error::non_finite_parameters);
        if(!x.allFinite())
            return unexpected(mrac_step_error::non_finite_state);
        if(!r.allFinite())
            return unexpected(mrac_step_error::non_finite_reference);
        return {};
    }

    // Only the faults that describe the controller rather than the cycle's arguments latch the health status.
    mrac_step_error latch_health(mrac_step_error fault)
    {
        if(fault == mrac_step_error::non_finite_reference_state || fault == mrac_step_error::non_finite_parameters)
            m_health = mrac_health::non_finite_carried_state;
        return fault;
    }

    // A policy outside the four known ones adapts nothing.
    bool gradient_enabled(const state_type& e) const
    {
        if constexpr(std::is_same_v<Robustification, dead_zone>)
            return compute_error_norm(e) > m_cfg.robustification.threshold;
        else
            return std::is_same_v<Robustification, no_robustification>;
    }

    expected<input_type, mrac_step_error> try_commit(
            state_type&& next_x_model, state_type&& next_tracking_error, const input_type& u, const input_type& e_proj, theta_x_type&& next_theta_x, theta_r_type&& next_theta_r)
    {
        if(!next_x_model.allFinite() || !next_tracking_error.allFinite() || !u.allFinite() || !e_proj.allFinite() || !next_theta_x.allFinite() || !next_theta_r.allFinite())
            return unexpected(mrac_step_error::non_finite_result);

        m_x_model = std::move(next_x_model);
        m_tracking_error = std::move(next_tracking_error);
        m_theta_x = std::move(next_theta_x);
        m_theta_r = std::move(next_theta_r);
        return u;
    }

    Scalar compute_error_norm(const state_type& e) const { return std::sqrt((e.transpose() * m_cfg.W * e)(0, 0)); }
};

}

#endif
