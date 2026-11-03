#ifndef HPP_GUARD_CTRLPP_CONTROL_L1_H
#define HPP_GUARD_CTRLPP_CONTROL_L1_H

// @cite hovakimyan2010 -- Hovakimyan & Cao, "L1 Adaptive Control Theory", 2010, Ch. 2

#include "ctrlpp/control/l1_config.h"

#include "ctrlpp/types.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/dsp/vector_biquad.h"
#include "ctrlpp/dsp/discrete_filter.h"

#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"

#include <cstddef>
#include <utility>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NX = 1, std::size_t NU = 1, typename Filter = vector_biquad<Scalar, NU>>
    requires vector_discrete_filter<Filter, Vector<Scalar, NU>>
class l1_controller
{
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");
    static_assert(NX == NU,
            "L1 reference feedforward requires NX == NU (a square predictor). The "
            "general non-square feedforward gain k_g = -(C (I - A_m)^{-1} B)^{-1} is "
            "not yet implemented");

public:
    using config_type = l1_config<Scalar, NX, NU>;
    using state_type = Vector<Scalar, NX>;
    using input_type = Vector<Scalar, NU>;

    static expected<l1_controller, l1_error> create(const config_type& cfg, Scalar cutoff_hz, Scalar sample_hz)
    {
        auto filter = Filter::low_pass(cutoff_hz, sample_hz);
        if(!filter.has_value())
            return unexpected(l1_error::invalid_filter_config);
        return create(cfg, *std::move(filter));
    }

    static expected<l1_controller, l1_error> create(const config_type& cfg, Filter filter)
    {
        auto k_r = compute_k_r(cfg);
        if(!k_r.has_value())
            return unexpected(k_r.error());
        return l1_controller{cfg, std::move(filter), *std::move(k_r)};
    }

    // The cycle is classified before any member is written, so a rejected cycle leaves every member bitwise unchanged.
    // The raw adaptation update is kept so its finiteness can be observed: Eigen's maxi and mini return the left operand
    // when the comparison is false (Eigen 3.4.1, Core/MathFunctions.h), so a finite projection bound replaces an
    // infinity invisibly downstream, while a NaN passes through.
    expected<input_type, l1_step_error> evaluate(const state_type& x, const input_type& r)
    {
        if(const auto step = check_step(x, r); !step)
            return unexpected(latch_health(step.error()));

        auto next_x_hat = propagate(m_cfg.predictor_model, m_x_hat, (m_u_prev + m_sigma_hat).eval());
        if(!next_x_hat.allFinite())
            return unexpected(l1_step_error::non_finite_result);

        auto next_x_tilde = (next_x_hat - x).eval();
        input_type sigma_raw = m_sigma_hat;
        sigma_raw.noalias() -= m_cfg.gamma * (m_cfg.predictor_model.B.transpose() * next_x_tilde);
        auto next_sigma_hat = sigma_raw.cwiseMax(m_cfg.theta_min).cwiseMin(m_cfg.theta_max).eval();
        bool const projection_substituted_non_finite = !sigma_raw.allFinite();
        auto u_raw = (m_k_r * r - next_sigma_hat).eval();
        if(!next_x_tilde.allFinite() || !next_sigma_hat.allFinite() || !u_raw.allFinite())
            return unexpected(l1_step_error::non_finite_result);

        Filter next_filter = m_filter;
        auto u_filtered = next_filter.process(u_raw);
        if(!u_filtered.allFinite())
            return unexpected(l1_step_error::non_finite_result);

        commit(next_x_hat, next_x_tilde, next_sigma_hat, std::move(next_filter), u_filtered, projection_substituted_non_finite);
        return u_filtered;
    }

    const state_type& x_hat() const { return m_x_hat; }

    const input_type& sigma_hat() const { return m_sigma_hat; }

    const state_type& tracking_error() const { return m_x_tilde; }

    const input_type& theta() const { return m_sigma_hat; }

    // Latches until reset and never downgrades; a cycle rejected on its arguments mutates nothing and does not set it.
    l1_health health() const { return m_health; }

    void reset()
    {
        m_x_hat = m_cfg.x_hat_0;
        m_sigma_hat = m_cfg.sigma_hat_0;
        m_x_tilde = state_type::Zero();
        m_u_prev = input_type::Zero();
        m_filter.reset();
        m_health = l1_health::ok;
    }

private:
    Filter m_filter;
    config_type m_cfg;
    state_type m_x_hat;
    l1_health m_health;
    input_type m_u_prev;
    state_type m_x_tilde;
    input_type m_sigma_hat;
    Matrix<Scalar, NU, NU> m_k_r;

    l1_controller(const config_type& cfg, Filter filter, Matrix<Scalar, NU, NU> k_r)
        : m_filter{std::move(filter)}
        , m_cfg{cfg}
        , m_x_hat{cfg.x_hat_0}
        , m_health{l1_health::ok}
        , m_u_prev{input_type::Zero()}
        , m_x_tilde{state_type::Zero()}
        , m_sigma_hat{cfg.sigma_hat_0}
        , m_k_r{std::move(k_r)}
    {
    }

    expected<void, l1_step_error> check_step(const state_type& x, const input_type& r) const
    {
        if(!m_x_hat.allFinite() || !m_u_prev.allFinite())
            return unexpected(l1_step_error::non_finite_predictor_state);
        if(!m_sigma_hat.allFinite())
            return unexpected(l1_step_error::non_finite_adaptation);
        if(!x.allFinite())
            return unexpected(l1_step_error::non_finite_state);
        if(!r.allFinite())
            return unexpected(l1_step_error::non_finite_reference);
        return {};
    }

    // Only the faults that describe the controller rather than the cycle's arguments latch the health status.
    l1_step_error latch_health(l1_step_error fault)
    {
        if(fault == l1_step_error::non_finite_predictor_state || fault == l1_step_error::non_finite_adaptation)
            note_health(l1_health::non_finite_carried_state);
        return fault;
    }

    void note_health(l1_health status)
    {
        if(status > m_health)
            m_health = status;
    }

    void commit(const state_type& next_x_hat, const state_type& next_x_tilde, const input_type& next_sigma_hat, Filter&& next_filter, const input_type& u_filtered,
            bool projection_substituted_non_finite)
    {
        m_x_hat = next_x_hat;
        m_x_tilde = next_x_tilde;
        m_sigma_hat = next_sigma_hat;
        m_filter = std::move(next_filter);
        m_u_prev = u_filtered;
        if(projection_substituted_non_finite)
            note_health(l1_health::projection_clamped_non_finite);
    }

    // K_r inverts the predictor DC gain (I - A_m)^{-1} B, so the steady state tracks the reference.
    static expected<Matrix<Scalar, NU, NU>, l1_error> compute_k_r(const config_type& cfg)
    {
        auto i_minus_a = (Matrix<Scalar, NX, NX>::Identity() - cfg.predictor_model.A).eval();
        auto lu_ima = i_minus_a.fullPivLu();
        if(!lu_ima.isInvertible())
            return unexpected(l1_error::singular_predictor);
        Matrix<Scalar, NX, NU> dc_gain = lu_ima.solve(cfg.predictor_model.B);
        auto lu_dc = dc_gain.fullPivLu();
        if(!lu_dc.isInvertible() || !dc_gain.allFinite())
            return unexpected(l1_error::singular_dc_gain);
        Matrix<Scalar, NU, NU> k_r = lu_dc.solve(Matrix<Scalar, NU, NU>::Identity());
        if(!k_r.allFinite())
            return unexpected(l1_error::non_finite_gain);
        return k_r;
    }
};

}

#endif
