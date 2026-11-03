#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_INPUT_FILTERS_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_INPUT_FILTERS_H

#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"

#include "ctrlpp/types.h"

#include <array>
#include <cstddef>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NY, typename... Policies>
class pid_input_filters
{
    using vector_t = Vector<Scalar, NY>;
    using config_type = pid_config<Scalar, NY, Policies...>;

public:
    pid_input_filters()
        : m_filtered_sp{vector_t::Zero()}
        , m_filtered_meas{vector_t::Zero()}
    {
    }

    vector_t apply_setpoint_filter(const config_type& cfg, const vector_t& sp, Scalar dt)
    {
        if constexpr(contains_v<setpoint_filter, Policies...>)
            return first_order_lag(m_filtered_sp, cfg.template policy<setpoint_filter>().tf, sp, dt);
        else
            return sp;
    }

    vector_t apply_pv_filter(const config_type& cfg, const vector_t& meas, Scalar dt)
    {
        if constexpr(contains_v<pv_filter, Policies...>)
            return first_order_lag(m_filtered_meas, cfg.template policy<pv_filter>().tf, meas, dt);
        else
            return meas;
    }

    bool finite() const { return (!contains_v<setpoint_filter, Policies...> || m_filtered_sp.allFinite()) && (!contains_v<pv_filter, Policies...> || m_filtered_meas.allFinite()); }

private:
    vector_t m_filtered_sp;
    vector_t m_filtered_meas;

    static vector_t first_order_lag(vector_t& state, const std::array<Scalar, NY>& tf, const vector_t& input, Scalar dt)
    {
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            if(tf[static_cast<std::size_t>(i)] > Scalar{0})
            {
                auto alpha = tf[static_cast<std::size_t>(i)] / (tf[static_cast<std::size_t>(i)] + dt);
                state[i] = alpha * state[i] + (Scalar{1} - alpha) * input[i];
            }
            else
                state[i] = input[i];
        }
        return state;
    }
};

}

#endif
