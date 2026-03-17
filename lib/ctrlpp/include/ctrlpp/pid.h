#ifndef HPP_GUARD_CTRLPP_PID_H
#define HPP_GUARD_CTRLPP_PID_H

// Element-wise loops use std::size_t but Eigen::operator[] takes signed
// Eigen::Index. The conversion is always safe (small non-negative values).
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-conversion"
#endif

#include "ctrlpp/pid_config.h"
#include "ctrlpp/pid_policies.h"
#include "ctrlpp/linalg_policy.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace ctrlpp {

namespace detail {

template<typename Scalar, std::size_t N, LinalgPolicy Policy>
constexpr auto make_zero_vector() -> typename Policy::template vector_type<Scalar, N>
{
    typename Policy::template vector_type<Scalar, N> v{};
    for (std::size_t i = 0; i < N; ++i)
        v[i] = Scalar{0};
    return v;
}

template<std::size_t N, typename Vec>
constexpr auto element_subtract(const Vec& a, const Vec& b) -> Vec
{
    Vec r{};
    for (std::size_t i = 0; i < N; ++i)
        r[i] = a[i] - b[i];
    return r;
}

template<std::size_t N, typename Vec>
constexpr auto element_add(const Vec& a, const Vec& b) -> Vec
{
    Vec r{};
    for (std::size_t i = 0; i < N; ++i)
        r[i] = a[i] + b[i];
    return r;
}

template<std::size_t N, typename Vec>
constexpr auto element_multiply(const Vec& a, const Vec& b) -> Vec
{
    Vec r{};
    for (std::size_t i = 0; i < N; ++i)
        r[i] = a[i] * b[i];
    return r;
}

template<std::size_t N, typename Vec, typename Scalar>
constexpr auto element_scale(const Vec& v, Scalar s) -> Vec
{
    Vec r{};
    for (std::size_t i = 0; i < N; ++i)
        r[i] = v[i] * s;
    return r;
}

template<std::size_t N, typename Vec>
constexpr auto element_clamp(const Vec& v, const Vec& lo, const Vec& hi) -> Vec
{
    Vec r{};
    for (std::size_t i = 0; i < N; ++i)
        r[i] = std::clamp(v[i], lo[i], hi[i]);
    return r;
}

template<std::size_t N, typename Vec, typename Scalar>
constexpr auto element_divide(const Vec& v, Scalar s) -> Vec
{
    Vec r{};
    for (std::size_t i = 0; i < N; ++i)
        r[i] = v[i] / s;
    return r;
}

template<std::size_t N, typename Vec>
constexpr bool any_not_equal(const Vec& a, const Vec& b)
{
    for (std::size_t i = 0; i < N; ++i)
        if (a[i] != b[i])
            return true;
    return false;
}

}

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY,
         LinalgPolicy Policy, typename... Policies>
class Pid {
public:
    using config_type = PidConfig<Scalar, NY, Policy, Policies...>;
    using vector_t = typename Policy::template vector_type<Scalar, NY>;

    explicit Pid(const config_type& cfg)
        : cfg_{cfg}
    {
        if constexpr (detail::contains_v<IsaForm, Policies...>) {
            for (std::size_t i = 0; i < NY; ++i) {
                if (cfg_.ki[i] != Scalar{0})
                    ki_[i] = cfg_.kp[i] / cfg_.ki[i];
                else
                    ki_[i] = Scalar{0};
                kd_[i] = cfg_.kp[i] * cfg_.kd[i];
                kp_[i] = cfg_.kp[i];
            }
        } else {
            kp_ = cfg_.kp;
            ki_ = cfg_.ki;
            kd_ = cfg_.kd;
        }

        if constexpr (detail::has_policy_v<AntiWindup, Policies...>) {
            using AW = detail::find_policy_t<AntiWindup, Policies...>;
            if constexpr (std::is_same_v<AW, AntiWindup<BackCalc>>) {
                const auto& aw_cfg = cfg_.template policy<AW>();
                for (std::size_t i = 0; i < NY; ++i) {
                    if (aw_cfg.kb[i] != Scalar{0}) {
                        kb_[i] = aw_cfg.kb[i];
                    } else {
                        if (kd_[i] != Scalar{0})
                            kb_[i] = std::sqrt(ki_[i] * kd_[i]);
                        else
                            kb_[i] = ki_[i];
                    }
                }
            }
        }
    }

    auto compute(const vector_t& sp, const vector_t& meas, Scalar dt) -> vector_t
    {
        if (dt <= Scalar{0})
            return prev_output_;

        // Input filtering: setpoint filter
        auto filtered_sp = sp;
        if constexpr (detail::contains_v<SetpointFilter, Policies...>) {
            const auto& tf_sp = cfg_.template policy<SetpointFilter>().tf;
            for (std::size_t i = 0; i < NY; ++i) {
                if (tf_sp[i] > Scalar{0}) {
                    auto alpha = tf_sp[i] / (tf_sp[i] + dt);
                    filtered_sp_[i] = alpha * filtered_sp_[i] + (Scalar{1} - alpha) * sp[i];
                } else {
                    filtered_sp_[i] = sp[i];
                }
            }
            filtered_sp = filtered_sp_;
        }

        // Input filtering: process variable filter
        auto filtered_meas = meas;
        if constexpr (detail::contains_v<PvFilter, Policies...>) {
            const auto& tf_pv = cfg_.template policy<PvFilter>().tf;
            for (std::size_t i = 0; i < NY; ++i) {
                if (tf_pv[i] > Scalar{0}) {
                    auto alpha = tf_pv[i] / (tf_pv[i] + dt);
                    filtered_meas_[i] = alpha * filtered_meas_[i] + (Scalar{1} - alpha) * meas[i];
                } else {
                    filtered_meas_[i] = meas[i];
                }
            }
            filtered_meas = filtered_meas_;
        }

        // Error (full, using filtered signals)
        auto e = detail::element_subtract<NY>(filtered_sp, filtered_meas);

        // Performance metric accumulation
        if constexpr (detail::has_policy_v<PerfAssessment, Policies...>) {
            using PA = detail::find_policy_t<PerfAssessment, Policies...>;
            accumulated_time_ += dt;

            if constexpr (detail::perf_has_metric_v<IAE, PA>) {
                for (std::size_t i = 0; i < NY; ++i) {
                    Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                    iae_[i] += abs_e * dt;
                }
            }
            if constexpr (detail::perf_has_metric_v<ISE, PA>) {
                for (std::size_t i = 0; i < NY; ++i)
                    ise_[i] += e[i] * e[i] * dt;
            }
            if constexpr (detail::perf_has_metric_v<ITAE, PA>) {
                for (std::size_t i = 0; i < NY; ++i) {
                    Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                    itae_[i] += accumulated_time_ * abs_e * dt;
                }
            }
            if constexpr (detail::perf_has_metric_v<OscillationDetect, PA>) {
                for (std::size_t i = 0; i < NY; ++i) {
                    Scalar sign_e = (e[i] > Scalar{0}) ? Scalar{1}
                                  : (e[i] < Scalar{0}) ? Scalar{-1}
                                  : Scalar{0};
                    if (!first_step_ && sign_e != Scalar{0} &&
                        prev_error_sign_[i] != Scalar{0} &&
                        sign_e != prev_error_sign_[i]) {
                        zero_crossings_[i] += Scalar{1};
                    }
                    if (sign_e != Scalar{0})
                        prev_error_sign_[i] = sign_e;
                }
            }
        }

        if constexpr (detail::contains_v<VelocityForm, Policies...>) {
            // Velocity (incremental) form: outputs delta_u per step
            // P contribution: Kp * (e(k) - e(k-1))
            auto dp = detail::element_multiply<NY>(
                kp_, detail::element_subtract<NY>(e, prev_error_));

            // I contribution: Ki * e(k) * dt (no accumulator — naturally anti-windup)
            auto di = detail::element_scale<NY>(
                detail::element_multiply<NY>(ki_, e), dt);

            // D contribution: Kd * (e(k) - 2*e(k-1) + e(k-2)) / dt
            auto d_num = detail::element_add<NY>(
                detail::element_subtract<NY>(e,
                    detail::element_scale<NY>(prev_error_, Scalar{2})),
                prev_prev_error_);
            auto dd = detail::element_multiply<NY>(
                kd_, detail::element_divide<NY>(d_num, dt));

            auto delta_u = detail::element_add<NY>(
                detail::element_add<NY>(dp, di), dd);

            // Feed-forward (adds to delta_u)
            if constexpr (detail::has_policy_v<FeedForward, Policies...>) {
                using ff_policy_t = detail::find_policy_t<FeedForward, Policies...>;
                if constexpr (!std::is_same_v<ff_policy_t, FeedForward<void>>) {
                    auto ff = cfg_.template policy<ff_policy_t>().ff_func(sp, dt);
                    delta_u = detail::element_add<NY>(delta_u, ff);
                }
            }

            // Output clamp on delta_u
            delta_u = detail::element_clamp<NY>(delta_u, cfg_.output_min, cfg_.output_max);

            // Update state
            prev_prev_error_ = prev_error_;
            prev_error_ = e;
            prev_meas_ = filtered_meas;
            prev_sp_ = filtered_sp;
            prev_output_ = delta_u;
            first_step_ = false;

            return delta_u;
        } else {
            // Position (absolute) form
            // P term with setpoint weighting: ep = b*sp - meas
            vector_t ep{};
            for (std::size_t i = 0; i < NY; ++i)
                ep[i] = cfg_.b[i] * filtered_sp[i] - filtered_meas[i];
            auto p = detail::element_multiply<NY>(kp_, ep);

            // I term -- always uses full error (sp - meas)
            auto integral_increment = detail::make_zero_vector<Scalar, NY, Policy>();
            if (!integral_frozen_) {
                if constexpr (detail::contains_v<ForwardEuler, Policies...>) {
                    auto ki_eprev = detail::element_multiply<NY>(ki_, prev_error_);
                    integral_increment = detail::element_scale<NY>(ki_eprev, dt);
                } else if constexpr (detail::contains_v<Tustin, Policies...>) {
                    auto avg = detail::element_scale<NY>(
                        detail::element_add<NY>(e, prev_error_),
                        Scalar{0.5});
                    auto ki_avg = detail::element_multiply<NY>(ki_, avg);
                    integral_increment = detail::element_scale<NY>(ki_avg, dt);
                } else {
                    auto ki_e = detail::element_multiply<NY>(ki_, e);
                    integral_increment = detail::element_scale<NY>(ki_e, dt);
                }
                integral_ = detail::element_add<NY>(integral_, integral_increment);
            }

            // D term with optional setpoint weighting (c parameter) and derivative filter
            vector_t d = detail::make_zero_vector<Scalar, NY, Policy>();
            if (!first_step_) {
                if (cfg_.derivative_on_error) {
                    vector_t ed_curr{};
                    vector_t ed_prev{};
                    for (std::size_t i = 0; i < NY; ++i) {
                        ed_curr[i] = cfg_.c[i] * filtered_sp[i] - filtered_meas[i];
                        ed_prev[i] = cfg_.c[i] * prev_sp_[i] - prev_meas_[i];
                    }
                    auto de = detail::element_subtract<NY>(ed_curr, ed_prev);
                    d = detail::element_multiply<NY>(kd_, detail::element_divide<NY>(de, dt));
                } else {
                    auto dm = detail::element_subtract<NY>(filtered_meas, prev_meas_);
                    auto dm_dt = detail::element_divide<NY>(dm, dt);
                    d = detail::element_subtract<NY>(
                        detail::make_zero_vector<Scalar, NY, Policy>(),
                        detail::element_multiply<NY>(kd_, dm_dt));
                }

                // Derivative filter (first-order low-pass)
                if constexpr (detail::contains_v<DerivFilter, Policies...>) {
                    const auto& df_cfg = cfg_.template policy<DerivFilter>();
                    for (std::size_t i = 0; i < NY; ++i) {
                        if (df_cfg.n[i] > Scalar{0} && kp_[i] != Scalar{0}) {
                            Scalar tf = kd_[i] / (kp_[i] * df_cfg.n[i]);
                            Scalar alpha = tf / (tf + dt);
                            d[i] = alpha * prev_deriv_filtered_[i] + (Scalar{1} - alpha) * d[i];
                        }
                    }
                    prev_deriv_filtered_ = d;
                }
            } else {
                if constexpr (detail::contains_v<DerivFilter, Policies...>) {
                    prev_deriv_filtered_ = d;
                }
            }

            // Sum
            auto u_raw = detail::element_add<NY>(
                detail::element_add<NY>(p, integral_), d);

            // Feed-forward
            if constexpr (detail::has_policy_v<FeedForward, Policies...>) {
                using ff_policy_t = detail::find_policy_t<FeedForward, Policies...>;
                if constexpr (!std::is_same_v<ff_policy_t, FeedForward<void>>) {
                    auto ff = cfg_.template policy<ff_policy_t>().ff_func(sp, dt);
                    u_raw = detail::element_add<NY>(u_raw, ff);
                }
            }

            // Rate limit (if RateLimit policy present)
            if constexpr (detail::contains_v<RateLimit, Policies...>) {
                const auto& rl_cfg = cfg_.template policy<RateLimit>();
                auto delta = detail::element_subtract<NY>(u_raw, prev_output_);
                vector_t max_delta{};
                vector_t neg_max_delta{};
                for (std::size_t i = 0; i < NY; ++i) {
                    if (rl_cfg.rate_max[i] > Scalar{0} &&
                        rl_cfg.rate_max[i] != std::numeric_limits<Scalar>::infinity()) {
                        max_delta[i] = rl_cfg.rate_max[i] * dt;
                        neg_max_delta[i] = -max_delta[i];
                    } else {
                        max_delta[i] = std::numeric_limits<Scalar>::max();
                        neg_max_delta[i] = std::numeric_limits<Scalar>::lowest();
                    }
                }
                auto clamped_delta = detail::element_clamp<NY>(delta, neg_max_delta, max_delta);
                u_raw = detail::element_add<NY>(prev_output_, clamped_delta);
            }

            // Output clamp
            auto u_sat = detail::element_clamp<NY>(u_raw, cfg_.output_min, cfg_.output_max);

            // Saturated flag
            saturated_ = detail::any_not_equal<NY>(u_sat, u_raw);

            // Anti-windup feedback (adjusts integral for next step)
            if constexpr (detail::has_policy_v<AntiWindup, Policies...>) {
                using AW = detail::find_policy_t<AntiWindup, Policies...>;
                if constexpr (std::is_same_v<AW, AntiWindup<BackCalc>>) {
                    auto sat_diff = detail::element_subtract<NY>(u_sat, u_raw);
                    auto feedback = detail::element_multiply<NY>(kb_, sat_diff);
                    integral_ = detail::element_add<NY>(
                        integral_, detail::element_scale<NY>(feedback, dt));
                } else if constexpr (std::is_same_v<AW, AntiWindup<Clamping>>) {
                    if (saturated_) {
                        for (std::size_t i = 0; i < NY; ++i) {
                            if ((e[i] > Scalar{0} && integral_[i] > Scalar{0}) ||
                                (e[i] < Scalar{0} && integral_[i] < Scalar{0})) {
                                integral_[i] -= integral_increment[i];
                            }
                        }
                    }
                } else if constexpr (std::is_same_v<AW, AntiWindup<ConditionalIntegration>>) {
                    const auto& ci_cfg = cfg_.template policy<AW>();
                    for (std::size_t i = 0; i < NY; ++i) {
                        Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                        if (abs_e > ci_cfg.error_threshold[i]) {
                            integral_[i] -= integral_increment[i];
                        }
                    }
                }
            }

            // Update state
            prev_error_ = e;
            prev_meas_ = filtered_meas;
            prev_sp_ = filtered_sp;
            prev_output_ = u_sat;
            first_step_ = false;

            return u_sat;
        }
    }

    auto compute(const vector_t& sp, const vector_t& meas, Scalar dt,
                 const vector_t& tracking_signal) -> vector_t
    {
        // Compute normal output first
        auto u = compute(sp, meas, dt);

        if (dt <= Scalar{0})
            return u;

        if constexpr (!detail::contains_v<VelocityForm, Policies...>) {
            // Tracking: set integral so that output matches tracking_signal
            // u = p + integral + d + ff
            // integral_target = tracking_signal - (u - integral)
            // i.e. integral = tracking_signal - p - d - ff
            auto non_integral = detail::element_subtract<NY>(u, integral_);
            integral_ = detail::element_subtract<NY>(tracking_signal, non_integral);
        }

        return u;
    }

    void set_params(const config_type& new_cfg)
    {
        // Bumpless integral rescaling: integral_new = integral_old * ki_old / ki_new
        vector_t ki_old = ki_;

        // Update config
        cfg_ = new_cfg;

        // Recompute internal gains
        if constexpr (detail::contains_v<IsaForm, Policies...>) {
            for (std::size_t i = 0; i < NY; ++i) {
                if (new_cfg.ki[i] != Scalar{0})
                    ki_[i] = new_cfg.kp[i] / new_cfg.ki[i];
                else
                    ki_[i] = Scalar{0};
                kd_[i] = new_cfg.kp[i] * new_cfg.kd[i];
                kp_[i] = new_cfg.kp[i];
            }
        } else {
            kp_ = new_cfg.kp;
            ki_ = new_cfg.ki;
            kd_ = new_cfg.kd;
        }

        // Rescale integral for bumpless gain change
        for (std::size_t i = 0; i < NY; ++i) {
            if (ki_[i] != Scalar{0} && ki_old[i] != Scalar{0})
                integral_[i] = integral_[i] * ki_old[i] / ki_[i];
            else if (ki_[i] == Scalar{0})
                integral_[i] = Scalar{0};
        }

        // Recompute BackCalc Kb if applicable
        if constexpr (detail::has_policy_v<AntiWindup, Policies...>) {
            using AW = detail::find_policy_t<AntiWindup, Policies...>;
            if constexpr (std::is_same_v<AW, AntiWindup<BackCalc>>) {
                const auto& aw_cfg = cfg_.template policy<AW>();
                for (std::size_t i = 0; i < NY; ++i) {
                    if (aw_cfg.kb[i] != Scalar{0}) {
                        kb_[i] = aw_cfg.kb[i];
                    } else {
                        if (kd_[i] != Scalar{0})
                            kb_[i] = std::sqrt(ki_[i] * kd_[i]);
                        else
                            kb_[i] = ki_[i];
                    }
                }
            }
        }
    }

    auto error() const -> const vector_t& { return prev_error_; }
    auto integral() const -> const vector_t& { return integral_; }
    auto params() const -> const config_type& { return cfg_; }
    auto saturated() const -> bool { return saturated_; }

    void reset()
    {
        integral_ = detail::make_zero_vector<Scalar, NY, Policy>();
        prev_error_ = detail::make_zero_vector<Scalar, NY, Policy>();
        prev_meas_ = detail::make_zero_vector<Scalar, NY, Policy>();
        prev_sp_ = detail::make_zero_vector<Scalar, NY, Policy>();
        prev_output_ = detail::make_zero_vector<Scalar, NY, Policy>();
        prev_prev_error_ = detail::make_zero_vector<Scalar, NY, Policy>();
        if constexpr (detail::contains_v<SetpointFilter, Policies...>) {
            filtered_sp_ = detail::make_zero_vector<Scalar, NY, Policy>();
        }
        if constexpr (detail::contains_v<PvFilter, Policies...>) {
            filtered_meas_ = detail::make_zero_vector<Scalar, NY, Policy>();
        }
        if constexpr (detail::contains_v<DerivFilter, Policies...>) {
            prev_deriv_filtered_ = detail::make_zero_vector<Scalar, NY, Policy>();
        }
        first_step_ = true;
        integral_frozen_ = false;
        saturated_ = false;
        if constexpr (detail::has_policy_v<PerfAssessment, Policies...>)
            reset_metrics();
    }

    void freeze_integral(bool freeze = true) { integral_frozen_ = freeze; }

    void set_integral(const vector_t& val) { integral_ = val; }

    template<typename Metric>
    auto metric() const -> const vector_t&
        requires detail::has_policy_v<PerfAssessment, Policies...>
    {
        using PA = detail::find_policy_t<PerfAssessment, Policies...>;
        static_assert(detail::perf_has_metric_v<Metric, PA>,
            "Metric type not in PerfAssessment pack");

        if constexpr (std::is_same_v<Metric, IAE>)
            return iae_;
        else if constexpr (std::is_same_v<Metric, ISE>)
            return ise_;
        else if constexpr (std::is_same_v<Metric, ITAE>)
            return itae_;
        else if constexpr (std::is_same_v<Metric, OscillationDetect>)
            return zero_crossings_;
    }

    auto oscillating() const -> bool
        requires detail::has_policy_v<PerfAssessment, Policies...>
    {
        using PA = detail::find_policy_t<PerfAssessment, Policies...>;
        static_assert(detail::perf_has_metric_v<OscillationDetect, PA>,
            "OscillationDetect not in PerfAssessment pack");

        if (accumulated_time_ <= Scalar{0})
            return false;

        for (std::size_t i = 0; i < NY; ++i) {
            if (zero_crossings_[i] / accumulated_time_ >
                static_cast<Scalar>(osc_threshold_))
                return true;
        }
        return false;
    }

    void reset_metrics()
        requires detail::has_policy_v<PerfAssessment, Policies...>
    {
        using PA = detail::find_policy_t<PerfAssessment, Policies...>;
        if constexpr (detail::perf_has_metric_v<IAE, PA>)
            iae_ = detail::make_zero_vector<Scalar, NY, Policy>();
        if constexpr (detail::perf_has_metric_v<ISE, PA>)
            ise_ = detail::make_zero_vector<Scalar, NY, Policy>();
        if constexpr (detail::perf_has_metric_v<ITAE, PA>)
            itae_ = detail::make_zero_vector<Scalar, NY, Policy>();
        if constexpr (detail::perf_has_metric_v<OscillationDetect, PA>) {
            zero_crossings_ = detail::make_zero_vector<Scalar, NY, Policy>();
            prev_error_sign_ = detail::make_zero_vector<Scalar, NY, Policy>();
        }
        accumulated_time_ = Scalar{0};
    }

private:
    config_type cfg_;
    vector_t kp_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t ki_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t kd_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t kb_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t integral_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t prev_error_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t prev_prev_error_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t prev_meas_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t prev_sp_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t prev_output_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t filtered_sp_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t filtered_meas_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t prev_deriv_filtered_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t iae_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t ise_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t itae_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t zero_crossings_ = detail::make_zero_vector<Scalar, NY, Policy>();
    vector_t prev_error_sign_ = detail::make_zero_vector<Scalar, NY, Policy>();
    Scalar accumulated_time_{0};
    Scalar osc_threshold_{5.0};
    bool first_step_{true};
    bool integral_frozen_{false};
    bool saturated_{false};
};

}

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

#endif
