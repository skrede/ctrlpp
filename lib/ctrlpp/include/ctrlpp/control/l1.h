#ifndef HPP_GUARD_CTRLPP_CONTROL_L1_H
#define HPP_GUARD_CTRLPP_CONTROL_L1_H

/// @brief L1 adaptive controller with state predictor, projection-based adaptation,
/// and low-pass filtered control output.
///
/// @cite hovakimyan2010 -- Hovakimyan & Cao, "L1 Adaptive Control Theory", 2010, Ch. 2

#include "ctrlpp/types.h"
#include "ctrlpp/config.h"
#include "ctrlpp/expected.h"

#include "ctrlpp/util/concepts.h"

#include "ctrlpp/control/l1_config.h"

#include "ctrlpp/dsp/vector_biquad.h"
#include "ctrlpp/dsp/discrete_filter.h"

#include "ctrlpp/model/propagate.h"
#include "ctrlpp/model/state_space.h"

#include <cstddef>
#include <utility>

namespace ctrlpp
{

template <ctrlpp_floating_scalar Scalar, std::size_t NX = 1, std::size_t NU = 1,
          typename Filter = vector_biquad<Scalar, NU>>
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

    /// Builds the output low-pass filter via `Filter::low_pass` and validates
    /// the predictor model. Returns `l1_error::invalid_filter_config` when the
    /// filter factory rejects the (cutoff_hz, sample_hz) design, or the
    /// predictor validation error of the filter-taking overload below.
    [[nodiscard]] static auto try_create(const config_type& cfg, Scalar cutoff_hz, Scalar sample_hz)
        -> expected<l1_controller, l1_error>
    {
        auto filter = Filter::low_pass(cutoff_hz, sample_hz);
        if(!filter.has_value())
            return unexpected(l1_error::invalid_filter_config);
        return try_create(cfg, *std::move(filter));
    }

    /// Validates the predictor model and constructs the controller. Returns
    /// `l1_error::singular_predictor` when (I - A_m) is singular,
    /// `l1_error::singular_dc_gain` when the DC gain (I - A_m)^{-1} B is
    /// singular or non-finite, and `l1_error::non_finite_gain` when the
    /// feedforward gain K_r is non-finite.
    [[nodiscard]] static auto try_create(const config_type& cfg, Filter filter)
        -> expected<l1_controller, l1_error>
    {
        auto k_r = compute_k_r(cfg);
        if(!k_r.has_value())
            return unexpected(k_r.error());
        return l1_controller{cfg, std::move(filter), *std::move(k_r)};
    }

#if CTRLPP_HAS_EXCEPTIONS
    /// Throwing convenience wrapper around `try_create`; throws the
    /// `bad_expected_access` of the active `ctrlpp::expected` target when the
    /// filter design or the predictor model is rejected.
    l1_controller(const config_type& cfg, Scalar cutoff_hz, Scalar sample_hz)
        : l1_controller{try_create(cfg, cutoff_hz, sample_hz).value()}
    {
    }

    /// Throwing convenience wrapper around `try_create`; throws the
    /// `bad_expected_access` of the active `ctrlpp::expected` target when the
    /// predictor model is rejected.
    l1_controller(const config_type& cfg, Filter filter)
        : l1_controller{try_create(cfg, std::move(filter)).value()}
    {
    }
#endif

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
    l1_controller(const config_type& cfg, Filter filter, Matrix<Scalar, NU, NU> k_r)
        : m_cfg{cfg}
        , m_filter{std::move(filter)}
        , m_x_hat{cfg.x_hat_0}
        , m_sigma_hat{cfg.sigma_hat_0}
        , m_k_r{std::move(k_r)}
    {
    }

    static auto compute_k_r(const config_type& cfg)
        -> expected<Matrix<Scalar, NU, NU>, l1_error>
    {
        // DC gain of predictor: G_dc = (I - A_m)^{-1} * B
        // K_r = G_dc^{-1} so that in steady state x_ss = r
        auto i_minus_a = (Matrix<Scalar, NX, NX>::Identity()
            - cfg.predictor_model.A).eval();
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
