#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_OUTPUT_STAGE_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_OUTPUT_STAGE_H

#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"

#include "ctrlpp/types.h"

#include <cstddef>
#include <type_traits>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NY, typename... Policies>
class pid_output_stage
{
    using vector_t = Vector<Scalar, NY>;
    using mask_t = Eigen::Array<bool, static_cast<int>(NY), 1>;
    using config_type = pid_config<Scalar, NY, Policies...>;

    struct limited_command
    {
        vector_t value;
        mask_t clipped;
    };

public:
    pid_output_stage()
        : m_prev_ff{vector_t::Zero()}
        , m_prev_output{vector_t::Zero()}
        , m_accumulated_output{vector_t::Zero()}
    {
    }

    void record(const vector_t& output) { m_prev_output = output; }

    // The actuator received the tracked command, not the automatic one.
    void track(const vector_t& target)
    {
        if constexpr(contains_v<velocity_form, Policies...>)
            m_accumulated_output = target;
        else
            m_prev_output = target;
    }

    limited_command apply_rate_limit(const config_type& cfg, const vector_t& u_raw, Scalar dt) const
    {
        if constexpr(!contains_v<rate_limit, Policies...>)
            return {u_raw, mask_t::Constant(false)};
        auto limited = clip_increment(cfg, (u_raw - m_prev_output).eval(), dt);
        // Re-forming p + (u - p) rounds, so a channel within its bound keeps its input bit for bit.
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            limited.value[i] = limited.clipped[i] ? m_prev_output[i] + limited.value[i] : u_raw[i];
        return limited;
    }

    // A zero rate leaves its channel unlimited, and an infinite one gives an infinite bound that never binds.
    limited_command clip_increment(const config_type& cfg, const vector_t& delta, Scalar dt) const
    {
        limited_command limited{delta, mask_t::Constant(false)};
        if constexpr(contains_v<rate_limit, Policies...>)
        {
            const auto& rate_max = cfg.template policy<rate_limit>().rate_max;
            for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            {
                const Scalar rate = rate_max[static_cast<std::size_t>(i)];
                const Scalar step = rate * dt;
                if(rate > Scalar{0} && (delta[i] > step || delta[i] < -step))
                {
                    limited.value[i] = delta[i] > step ? step : -step;
                    limited.clipped[i] = true;
                }
            }
        }
        return limited;
    }

    // Clamping the accumulated output rather than the increment keeps it in range while still letting it
    // rise and fall against an asymmetric limit. The clamp returns its input when inside, so != is exact.
    vector_t accumulate_output(const config_type& cfg, const vector_t& delta_u, bool& saturated)
    {
        auto target = (m_accumulated_output + delta_u).eval();
        auto clamped = target.cwiseMax(cfg.output_min).cwiseMin(cfg.output_max).eval();
        saturated = (clamped.array() != target.array()).any();
        auto emitted = (clamped - m_accumulated_output).eval();
        m_accumulated_output = clamped;
        return emitted;
    }

    // The velocity form emits increments, so it injects the change in the feed-forward level, not the level:
    // a constant feed-forward then contributes nothing to the increment and the actuator does not drift.
    vector_t apply_feed_forward_velocity(config_type& cfg, vector_t delta_u, const vector_t& sp, Scalar dt)
    {
        if constexpr(has_policy_v<feed_forward, Policies...>)
        {
            using ff_policy_t = find_policy_t<feed_forward, Policies...>;
            if constexpr(!std::is_same_v<ff_policy_t, feed_forward<void>>)
            {
                auto ff = cfg.template policy<ff_policy_t>().ff_func(sp, dt).eval();
                auto delta_ff = (ff - m_prev_ff).eval();
                m_prev_ff = ff;
                return (delta_u + delta_ff).eval();
            }
        }
        return delta_u;
    }

    bool finite() const
    {
        return (!contains_v<velocity_form, Policies...> || m_accumulated_output.allFinite()) && (!contains_v<rate_limit, Policies...> || m_prev_output.allFinite()) &&
                (!has_policy_v<feed_forward, Policies...> || m_prev_ff.allFinite());
    }

private:
    vector_t m_prev_ff;
    vector_t m_prev_output;
    vector_t m_accumulated_output;
};

}

#endif
