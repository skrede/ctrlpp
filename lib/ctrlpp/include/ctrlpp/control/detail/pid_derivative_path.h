#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_DERIVATIVE_PATH_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_DERIVATIVE_PATH_H

#include "ctrlpp/control/detail/pid_gains.h"

#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"

#include "ctrlpp/types.h"

#include <cstddef>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NY, typename... Policies>
class pid_derivative_path
{
    using vector_t = Vector<Scalar, NY>;
    using config_type = pid_config<Scalar, NY, Policies...>;
    using gains_type = pid_gains<Scalar, NY>;

public:
    pid_derivative_path()
        : m_prev_derivative{vector_t::Zero()}
        , m_prev_deriv_filtered{vector_t::Zero()}
    {
    }

    // @cite astrom2006 -- Ch. 3.4 (derivative action, setpoint weighting c parameter)
    vector_t compute_derivative_term(
            const config_type& cfg, const gains_type& gains, const vector_t& sp, const vector_t& meas, const vector_t& prev_sp, const vector_t& prev_meas, bool first_step, Scalar dt)
    {
        vector_t d = vector_t::Zero();
        if(!first_step)
        {
            d = compute_raw_derivative(cfg, gains.kd, sp, meas, prev_sp, prev_meas, dt);
            d = apply_derivative_filter(gains.tf, d, dt);
        }
        else if constexpr(contains_v<deriv_filter, Policies...>)
            m_prev_deriv_filtered = d;
        return d;
    }

    vector_t increment(const vector_t& d)
    {
        auto delta = (d - m_prev_derivative).eval();
        m_prev_derivative = d;
        return delta;
    }

    bool finite() const
    {
        return (!contains_v<velocity_form, Policies...> || m_prev_derivative.allFinite()) && (!contains_v<deriv_filter, Policies...> || m_prev_deriv_filtered.allFinite());
    }

private:
    vector_t m_prev_derivative;
    vector_t m_prev_deriv_filtered;

    static vector_t compute_raw_derivative(
            const config_type& cfg, const vector_t& kd, const vector_t& sp, const vector_t& meas, const vector_t& prev_sp, const vector_t& prev_meas, Scalar dt)
    {
        if(cfg.derivative_on_error)
        {
            vector_t ed_curr, ed_prev;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                ed_curr[i] = cfg.c[i] * sp[i] - meas[i];
                ed_prev[i] = cfg.c[i] * prev_sp[i] - prev_meas[i];
            }
            return kd.cwiseProduct((ed_curr - ed_prev) / dt).eval();
        }
        else
        {
            auto dm_dt = ((meas - prev_meas) / dt).eval();
            return (-kd.cwiseProduct(dm_dt)).eval();
        }
    }

    vector_t apply_derivative_filter(const vector_t& tf, vector_t d, Scalar dt)
    {
        if constexpr(contains_v<deriv_filter, Policies...>)
        {
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
                if(tf[i] > Scalar{0})
                {
                    Scalar alpha = tf[i] / (tf[i] + dt);
                    d[i] = alpha * m_prev_deriv_filtered[i] + (Scalar{1} - alpha) * d[i];
                }
            m_prev_deriv_filtered = d;
        }
        return d;
    }
};

}

#endif
