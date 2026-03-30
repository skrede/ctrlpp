#ifndef HPP_GUARD_CTRLPP_CONTROL_L1_H
#define HPP_GUARD_CTRLPP_CONTROL_L1_H

/// @brief L1 adaptive controller with state predictor, projection-based adaptation,
/// and low-pass filtered control output.
///
/// @cite hovakimyan2010 -- Hovakimyan & Cao, "L1 Adaptive Control Theory", 2010, Ch. 2

#include "ctrlpp/types.h"

#include "ctrlpp/control/l1_config.h"

#include "ctrlpp/dsp/vector_biquad.h"
#include "ctrlpp/dsp/discrete_filter.h"

#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <type_traits>

namespace ctrlpp
{

template <typename Scalar, std::size_t NX = 1, std::size_t NU = 1,
          typename Filter = vector_biquad<Scalar, NU>>
    requires vector_discrete_filter<Filter, Vector<Scalar, NU>>
class l1_controller
{
    static_assert(std::is_floating_point_v<Scalar>, "Scalar must be a floating-point type");
    static_assert(NX > 0, "State dimension NX must be positive");
    static_assert(NU > 0, "Input dimension NU must be positive");

public:
    using config_type = l1_config<Scalar, NX, NU>;
    using state_type = Vector<Scalar, NX>;
    using input_type = Vector<Scalar, NU>;

    l1_controller(const config_type& cfg, Scalar cutoff_hz, Scalar sample_hz)
        : m_cfg{cfg}
        , m_filter{Filter::low_pass(cutoff_hz, sample_hz)}
        , m_x_hat{cfg.x_hat_0}
        , m_sigma_hat{cfg.sigma_hat_0}
    {
        compute_k_r();
    }

    l1_controller(const config_type& cfg, Filter filter)
        : m_cfg{cfg}
        , m_filter{std::move(filter)}
        , m_x_hat{cfg.x_hat_0}
        , m_sigma_hat{cfg.sigma_hat_0}
    {
        compute_k_r();
    }

    auto evaluate(const state_type& x, const input_type& r) -> input_type
    {
        // 1. State predictor: x_hat = A_m * x_hat + B * (u_prev + sigma_hat)
        m_x_hat = propagate(m_cfg.predictor_model, m_x_hat,
                            (m_u_prev + m_sigma_hat).eval());

        // 2. Prediction error (Hovakimyan convention)
        m_x_tilde = m_x_hat - x;

        // 3. Adaptation with projection (elementwise clamp)
        m_sigma_hat.noalias() -= m_cfg.gamma
            * (m_cfg.predictor_model.B.transpose() * m_x_tilde);
        m_sigma_hat = m_sigma_hat.cwiseMax(m_cfg.theta_min).cwiseMin(m_cfg.theta_max);

        // 4. Raw control: reference feedforward minus uncertainty estimate
        auto u_raw = (m_k_r * r - m_sigma_hat).eval();

        // 5. Low-pass filter (L1 robustification mechanism)
        auto u_filtered = m_filter.process(u_raw);

        // 6. Store for next predictor step
        m_u_prev = u_filtered;

        return u_filtered;
    }

    auto x_hat() const -> const state_type& { return m_x_hat; }

    auto sigma_hat() const -> const input_type& { return m_sigma_hat; }

    auto tracking_error() const -> const state_type& { return m_x_tilde; }

    auto theta() const -> const input_type& { return m_sigma_hat; }

    void reset()
    {
        m_x_hat = m_cfg.x_hat_0;
        m_sigma_hat = m_cfg.sigma_hat_0;
        m_x_tilde = state_type::Zero();
        m_u_prev = input_type::Zero();
        m_filter.reset();
    }

private:
    void compute_k_r()
    {
        if constexpr(NX == NU)
        {
            // DC gain of predictor: G_dc = (I - A_m)^{-1} * B
            // K_r = G_dc^{-1} so that in steady state x_ss = r
            auto i_minus_a = (Matrix<Scalar, NX, NX>::Identity()
                - m_cfg.predictor_model.A).eval();
            auto lu_ima = i_minus_a.fullPivLu();
            if(!lu_ima.isInvertible())
                throw std::invalid_argument(
                    "L1 predictor model has unit eigenvalue: (I - A) is singular");
            Matrix<Scalar, NX, NU> dc_gain = lu_ima.solve(m_cfg.predictor_model.B);
            auto lu_dc = dc_gain.fullPivLu();
            if(!lu_dc.isInvertible() || !dc_gain.allFinite())
                throw std::invalid_argument(
                    "L1 predictor model has near-zero DC gain: B / (I - A) is singular");
            m_k_r = lu_dc.solve(Matrix<Scalar, NU, NU>::Identity());
            if(!m_k_r.allFinite())
                throw std::invalid_argument(
                    "L1 feedforward gain K_r is non-finite: ill-conditioned predictor model");
        }
    }

    config_type m_cfg;
    Filter m_filter;
    state_type m_x_hat = state_type::Zero();
    state_type m_x_tilde = state_type::Zero();
    input_type m_sigma_hat = input_type::Zero();
    input_type m_u_prev = input_type::Zero();
    Matrix<Scalar, NU, NU> m_k_r = Matrix<Scalar, NU, NU>::Identity();
};

}

#endif
