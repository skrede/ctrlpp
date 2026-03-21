#ifndef HPP_GUARD_CTRLPP_PID_H
#define HPP_GUARD_CTRLPP_PID_H

#include "ctrlpp/pid_config.h"
#include "ctrlpp/pid_policies.h"
#include "ctrlpp/types.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>

namespace ctrlpp {

template<typename Scalar, std::size_t NX, std::size_t NU, std::size_t NY,
         typename... Policies>
class pid {
public:
    using config_type = pid_config<Scalar, NY, Policies...>;
    using vector_t = Vector<Scalar, NY>;

    explicit pid(const config_type& cfg)
        : cfg_{cfg}
    {
        if constexpr (detail::contains_v<isa_form, Policies...>) {
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
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

        if constexpr (detail::has_policy_v<anti_windup, Policies...>) {
            using AW = detail::find_policy_t<anti_windup, Policies...>;
            if constexpr (std::is_same_v<AW, anti_windup<back_calc>>) {
                const auto& aw_cfg = cfg_.template policy<AW>();
                for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                    if (aw_cfg.kb[static_cast<std::size_t>(i)] != Scalar{0}) {
                        kb_[i] = aw_cfg.kb[static_cast<std::size_t>(i)];
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
        if constexpr (detail::contains_v<setpoint_filter, Policies...>) {
            const auto& tf_sp = cfg_.template policy<setpoint_filter>().tf;
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                if (tf_sp[static_cast<std::size_t>(i)] > Scalar{0}) {
                    auto alpha = tf_sp[static_cast<std::size_t>(i)] / (tf_sp[static_cast<std::size_t>(i)] + dt);
                    filtered_sp_[i] = alpha * filtered_sp_[i] + (Scalar{1} - alpha) * sp[i];
                } else {
                    filtered_sp_[i] = sp[i];
                }
            }
            filtered_sp = filtered_sp_;
        }

        // Input filtering: process variable filter
        auto filtered_meas = meas;
        if constexpr (detail::contains_v<pv_filter, Policies...>) {
            const auto& tf_pv = cfg_.template policy<pv_filter>().tf;
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                if (tf_pv[static_cast<std::size_t>(i)] > Scalar{0}) {
                    auto alpha = tf_pv[static_cast<std::size_t>(i)] / (tf_pv[static_cast<std::size_t>(i)] + dt);
                    filtered_meas_[i] = alpha * filtered_meas_[i] + (Scalar{1} - alpha) * meas[i];
                } else {
                    filtered_meas_[i] = meas[i];
                }
            }
            filtered_meas = filtered_meas_;
        }

        // Error (full, using filtered signals)
        auto e = (filtered_sp - filtered_meas).eval();

        // Performance metric accumulation
        if constexpr (detail::has_policy_v<perf_assessment, Policies...>) {
            using PA = detail::find_policy_t<perf_assessment, Policies...>;
            accumulated_time_ += dt;

            if constexpr (detail::perf_has_metric_v<IAE, PA>) {
                for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                    Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                    iae_[i] += abs_e * dt;
                }
            }
            if constexpr (detail::perf_has_metric_v<ISE, PA>) {
                for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                    ise_[i] += e[i] * e[i] * dt;
            }
            if constexpr (detail::perf_has_metric_v<ITAE, PA>) {
                for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                    Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                    itae_[i] += accumulated_time_ * abs_e * dt;
                }
            }
            if constexpr (detail::perf_has_metric_v<oscillation_detect, PA>) {
                for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
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

        if constexpr (detail::contains_v<velocity_form, Policies...>) {
            // Velocity (incremental) form: outputs delta_u per step
            // P contribution: Kp * (e(k) - e(k-1))
            auto dp = kp_.cwiseProduct(e - prev_error_).eval();

            // I contribution: Ki * e(k) * dt (no accumulator — naturally anti-windup)
            auto di = (ki_.cwiseProduct(e) * dt).eval();

            // D contribution: Kd * (e(k) - 2*e(k-1) + e(k-2)) / dt
            auto d_num = (e - prev_error_ * Scalar{2} + prev_prev_error_).eval();
            auto dd = kd_.cwiseProduct(d_num / dt).eval();

            auto delta_u = (dp + di + dd).eval();

            // Feed-forward (adds to delta_u)
            if constexpr (detail::has_policy_v<feed_forward, Policies...>) {
                using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
                if constexpr (!std::is_same_v<ff_policy_t, feed_forward<void>>) {
                    auto ff = cfg_.template policy<ff_policy_t>().ff_func(sp, dt);
                    delta_u = (delta_u + ff).eval();
                }
            }

            // Output clamp on delta_u
            delta_u = delta_u.cwiseMax(cfg_.output_min).cwiseMin(cfg_.output_max).eval();

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
            vector_t ep;
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                ep[i] = cfg_.b[i] * filtered_sp[i] - filtered_meas[i];
            auto p = kp_.cwiseProduct(ep).eval();

            // I term -- always uses full error (sp - meas)
            vector_t integral_increment = vector_t::Zero();
            if (!integral_frozen_) {
                if constexpr (detail::contains_v<forward_euler, Policies...>) {
                    integral_increment = (ki_.cwiseProduct(prev_error_) * dt).eval();
                } else if constexpr (detail::contains_v<tustin, Policies...>) {
                    auto avg = ((e + prev_error_) * Scalar{0.5}).eval();
                    integral_increment = (ki_.cwiseProduct(avg) * dt).eval();
                } else {
                    integral_increment = (ki_.cwiseProduct(e) * dt).eval();
                }
                integral_ = (integral_ + integral_increment).eval();
            }

            // D term with optional setpoint weighting (c parameter) and derivative filter
            vector_t d = vector_t::Zero();
            if (!first_step_) {
                if (cfg_.derivative_on_error) {
                    vector_t ed_curr;
                    vector_t ed_prev;
                    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                        ed_curr[i] = cfg_.c[i] * filtered_sp[i] - filtered_meas[i];
                        ed_prev[i] = cfg_.c[i] * prev_sp_[i] - prev_meas_[i];
                    }
                    auto de = (ed_curr - ed_prev).eval();
                    d = kd_.cwiseProduct(de / dt).eval();
                } else {
                    auto dm = (filtered_meas - prev_meas_).eval();
                    auto dm_dt = (dm / dt).eval();
                    d = (-kd_.cwiseProduct(dm_dt)).eval();
                }

                // Derivative filter (first-order low-pass)
                if constexpr (detail::contains_v<deriv_filter, Policies...>) {
                    const auto& df_cfg = cfg_.template policy<deriv_filter>();
                    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                        if (df_cfg.n[static_cast<std::size_t>(i)] > Scalar{0} && kp_[i] != Scalar{0}) {
                            Scalar tf = kd_[i] / (kp_[i] * df_cfg.n[static_cast<std::size_t>(i)]);
                            Scalar alpha = tf / (tf + dt);
                            d[i] = alpha * prev_deriv_filtered_[i] + (Scalar{1} - alpha) * d[i];
                        }
                    }
                    prev_deriv_filtered_ = d;
                }
            } else {
                if constexpr (detail::contains_v<deriv_filter, Policies...>) {
                    prev_deriv_filtered_ = d;
                }
            }

            // Sum
            auto u_raw = (p + integral_ + d).eval();

            // Feed-forward
            if constexpr (detail::has_policy_v<feed_forward, Policies...>) {
                using ff_policy_t = detail::find_policy_t<feed_forward, Policies...>;
                if constexpr (!std::is_same_v<ff_policy_t, feed_forward<void>>) {
                    auto ff = cfg_.template policy<ff_policy_t>().ff_func(sp, dt);
                    u_raw = (u_raw + ff).eval();
                }
            }

            // Rate limit (if rate_limit policy present)
            if constexpr (detail::contains_v<rate_limit, Policies...>) {
                const auto& rl_cfg = cfg_.template policy<rate_limit>();
                auto delta = (u_raw - prev_output_).eval();
                vector_t max_delta;
                vector_t neg_max_delta;
                for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                    if (rl_cfg.rate_max[static_cast<std::size_t>(i)] > Scalar{0} &&
                        rl_cfg.rate_max[static_cast<std::size_t>(i)] != std::numeric_limits<Scalar>::infinity()) {
                        max_delta[i] = rl_cfg.rate_max[static_cast<std::size_t>(i)] * dt;
                        neg_max_delta[i] = -max_delta[i];
                    } else {
                        max_delta[i] = std::numeric_limits<Scalar>::max();
                        neg_max_delta[i] = std::numeric_limits<Scalar>::lowest();
                    }
                }
                auto clamped_delta = delta.cwiseMax(neg_max_delta).cwiseMin(max_delta).eval();
                u_raw = (prev_output_ + clamped_delta).eval();
            }

            // Output clamp
            auto u_sat = u_raw.cwiseMax(cfg_.output_min).cwiseMin(cfg_.output_max).eval();

            // Saturated flag
            saturated_ = (u_sat.array() != u_raw.array()).any();

            // Anti-windup feedback (adjusts integral for next step)
            if constexpr (detail::has_policy_v<anti_windup, Policies...>) {
                using AW = detail::find_policy_t<anti_windup, Policies...>;
                if constexpr (std::is_same_v<AW, anti_windup<back_calc>>) {
                    auto sat_diff = (u_sat - u_raw).eval();
                    auto feedback = kb_.cwiseProduct(sat_diff).eval();
                    integral_ = (integral_ + feedback * dt).eval();
                } else if constexpr (std::is_same_v<AW, anti_windup<clamping>>) {
                    if (saturated_) {
                        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                            if ((e[i] > Scalar{0} && integral_[i] > Scalar{0}) ||
                                (e[i] < Scalar{0} && integral_[i] < Scalar{0})) {
                                integral_[i] -= integral_increment[i];
                            }
                        }
                    }
                } else if constexpr (std::is_same_v<AW, anti_windup<conditional_integration>>) {
                    const auto& ci_cfg = cfg_.template policy<AW>();
                    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                        Scalar abs_e = e[i] < Scalar{0} ? -e[i] : e[i];
                        if (abs_e > ci_cfg.error_threshold[static_cast<std::size_t>(i)]) {
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

        if constexpr (!detail::contains_v<velocity_form, Policies...>) {
            // Tracking: set integral so that output matches tracking_signal
            auto non_integral = (u - integral_).eval();
            integral_ = (tracking_signal - non_integral).eval();
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
        if constexpr (detail::contains_v<isa_form, Policies...>) {
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
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
        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
            if (ki_[i] != Scalar{0} && ki_old[i] != Scalar{0})
                integral_[i] = integral_[i] * ki_old[i] / ki_[i];
            else if (ki_[i] == Scalar{0})
                integral_[i] = Scalar{0};
        }

        // Recompute back_calc Kb if applicable
        if constexpr (detail::has_policy_v<anti_windup, Policies...>) {
            using AW = detail::find_policy_t<anti_windup, Policies...>;
            if constexpr (std::is_same_v<AW, anti_windup<back_calc>>) {
                const auto& aw_cfg = cfg_.template policy<AW>();
                for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
                    if (aw_cfg.kb[static_cast<std::size_t>(i)] != Scalar{0}) {
                        kb_[i] = aw_cfg.kb[static_cast<std::size_t>(i)];
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
        integral_ = vector_t::Zero();
        prev_error_ = vector_t::Zero();
        prev_meas_ = vector_t::Zero();
        prev_sp_ = vector_t::Zero();
        prev_output_ = vector_t::Zero();
        prev_prev_error_ = vector_t::Zero();
        if constexpr (detail::contains_v<setpoint_filter, Policies...>) {
            filtered_sp_ = vector_t::Zero();
        }
        if constexpr (detail::contains_v<pv_filter, Policies...>) {
            filtered_meas_ = vector_t::Zero();
        }
        if constexpr (detail::contains_v<deriv_filter, Policies...>) {
            prev_deriv_filtered_ = vector_t::Zero();
        }
        first_step_ = true;
        integral_frozen_ = false;
        saturated_ = false;
        if constexpr (detail::has_policy_v<perf_assessment, Policies...>)
            reset_metrics();
    }

    void freeze_integral(bool freeze = true) { integral_frozen_ = freeze; }

    void set_integral(const vector_t& val) { integral_ = val; }

    template<typename Metric>
    auto metric() const -> const vector_t&
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        using PA = detail::find_policy_t<perf_assessment, Policies...>;
        static_assert(detail::perf_has_metric_v<Metric, PA>,
            "Metric type not in perf_assessment pack");

        if constexpr (std::is_same_v<Metric, IAE>)
            return iae_;
        else if constexpr (std::is_same_v<Metric, ISE>)
            return ise_;
        else if constexpr (std::is_same_v<Metric, ITAE>)
            return itae_;
        else if constexpr (std::is_same_v<Metric, oscillation_detect>)
            return zero_crossings_;
    }

    auto oscillating() const -> bool
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        using PA = detail::find_policy_t<perf_assessment, Policies...>;
        static_assert(detail::perf_has_metric_v<oscillation_detect, PA>,
            "oscillation_detect not in perf_assessment pack");

        if (accumulated_time_ <= Scalar{0})
            return false;

        for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i) {
            if (zero_crossings_[i] / accumulated_time_ >
                static_cast<Scalar>(osc_threshold_))
                return true;
        }
        return false;
    }

    void reset_metrics()
        requires detail::has_policy_v<perf_assessment, Policies...>
    {
        using PA = detail::find_policy_t<perf_assessment, Policies...>;
        if constexpr (detail::perf_has_metric_v<IAE, PA>)
            iae_ = vector_t::Zero();
        if constexpr (detail::perf_has_metric_v<ISE, PA>)
            ise_ = vector_t::Zero();
        if constexpr (detail::perf_has_metric_v<ITAE, PA>)
            itae_ = vector_t::Zero();
        if constexpr (detail::perf_has_metric_v<oscillation_detect, PA>) {
            zero_crossings_ = vector_t::Zero();
            prev_error_sign_ = vector_t::Zero();
        }
        accumulated_time_ = Scalar{0};
    }

private:
    config_type cfg_;
    vector_t kp_ = vector_t::Zero();
    vector_t ki_ = vector_t::Zero();
    vector_t kd_ = vector_t::Zero();
    vector_t kb_ = vector_t::Zero();
    vector_t integral_ = vector_t::Zero();
    vector_t prev_error_ = vector_t::Zero();
    vector_t prev_prev_error_ = vector_t::Zero();
    vector_t prev_meas_ = vector_t::Zero();
    vector_t prev_sp_ = vector_t::Zero();
    vector_t prev_output_ = vector_t::Zero();
    vector_t filtered_sp_ = vector_t::Zero();
    vector_t filtered_meas_ = vector_t::Zero();
    vector_t prev_deriv_filtered_ = vector_t::Zero();
    vector_t iae_ = vector_t::Zero();
    vector_t ise_ = vector_t::Zero();
    vector_t itae_ = vector_t::Zero();
    vector_t zero_crossings_ = vector_t::Zero();
    vector_t prev_error_sign_ = vector_t::Zero();
    Scalar accumulated_time_{0};
    Scalar osc_threshold_{5.0};
    bool first_step_{true};
    bool integral_frozen_{false};
    bool saturated_{false};
};

}

#endif
