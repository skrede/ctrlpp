#ifndef HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_INTEGRATOR_H
#define HPP_GUARD_CTRLPP_CONTROL_DETAIL_PID_INTEGRATOR_H

#include "ctrlpp/control/detail/pid_gains.h"

#include "ctrlpp/control/pid_config.h"
#include "ctrlpp/control/pid_policies.h"

#include "ctrlpp/types.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>

namespace ctrlpp::detail
{

template <typename Scalar, std::size_t NY, typename... Policies>
class pid_integrator
{
    using vector_t = Vector<Scalar, NY>;
    using mask_t = Eigen::Array<bool, static_cast<int>(NY), 1>;
    using direction_t = Eigen::Array<std::int8_t, static_cast<int>(NY), 1>;
    using config_type = pid_config<Scalar, NY, Policies...>;
    using gains_type = pid_gains<Scalar, NY>;

public:
    pid_integrator()
        : m_integral{vector_t::Zero()}
        , m_constraint_direction{direction_t::Zero()}
    {
    }

    const vector_t& integral() const { return m_integral; }

    void assign(const vector_t& value) { m_integral = value; }

    void track(vector_t integral)
    {
        m_integral = std::move(integral);
        m_constraint_direction = direction_t::Zero();
    }

    // @cite astrom2006 -- Ch. 3.3 (integral action, discretization methods)
    vector_t integral_increment(const gains_type& gains, bool frozen, const vector_t& prev_error, const vector_t& e, Scalar dt) const
    {
        if(frozen)
            return vector_t::Zero();
        if constexpr(contains_v<forward_euler, Policies...>)
            return (gains.ki.cwiseProduct(prev_error) * dt).eval();
        else if constexpr(contains_v<tustin, Policies...>)
            return (gains.ki.cwiseProduct(((e + prev_error) * Scalar{0.5}).eval()) * dt).eval();
        else
            return (gains.ki.cwiseProduct(e) * dt).eval();
    }

    vector_t compute_integral_term(const config_type& cfg, const gains_type& gains, bool frozen, const vector_t& prev_error, const vector_t& e, Scalar dt)
    {
        if(!frozen)
            m_integral = (m_integral + gate_increment(cfg, integral_increment(gains, frozen, prev_error, e, dt), e)).eval();
        return m_integral;
    }

    // A frozen integrator is held against back-calculation too, while the clamping verdict is recorded regardless,
    // since it describes the emitted command rather than the integrator.
    void apply_anti_windup(const gains_type& gains, bool frozen, const vector_t& u_sat, const vector_t& u_unconstrained, const mask_t& constrained, Scalar dt)
    {
        if constexpr(contains_v<anti_windup<clamping>, Policies...>)
            m_constraint_direction = constraint_direction(u_sat, u_unconstrained, constrained);
        else if constexpr(contains_v<anti_windup<back_calc>, Policies...>)
        {
            if(!frozen)
                apply_back_calculation(gains.kb, u_sat, u_unconstrained, dt);
        }
    }

    bool finite() const { return contains_v<velocity_form, Policies...> || m_integral.allFinite(); }

private:
    vector_t m_integral;
    direction_t m_constraint_direction;

    // Clamping and conditional integration decide before the output is formed, so no increment that formed a
    // command is removed afterward; clamping therefore reads the previous cycle's constraint.
    vector_t gate_increment(const config_type& cfg, vector_t increment, const vector_t& e) const
    {
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
        {
            if constexpr(contains_v<anti_windup<clamping>, Policies...>)
            {
                const std::int8_t direction = m_constraint_direction[i];
                if((direction > 0 && increment[i] > Scalar{0}) || (direction < 0 && increment[i] < Scalar{0}))
                    increment[i] = Scalar{0};
            }
            else if constexpr(contains_v<anti_windup<conditional_integration>, Policies...>)
            {
                const auto& threshold = cfg.template policy<anti_windup<conditional_integration>>().error_threshold;
                if(std::abs(e[i]) > threshold[static_cast<std::size_t>(i)])
                    increment[i] = Scalar{0};
            }
        }
        return increment;
    }

    void apply_back_calculation(const vector_t& kb, const vector_t& u_sat, const vector_t& u_unconstrained, Scalar dt)
    {
        auto feedback = kb.cwiseProduct((u_sat - u_unconstrained).eval()).eval();
        m_integral = (m_integral + feedback * dt).eval();
    }

    // @cite caparroz2026 -- clamping decided from the previous sample's saturation error, Algorithm 2. The
    // verdict is taken from the masks; the unconstrained command only names the side of the violated limit.
    static direction_t constraint_direction(const vector_t& u_sat, const vector_t& u_unconstrained, const mask_t& constrained)
    {
        direction_t direction = direction_t::Zero();
        for(Eigen::Index i = 0; i < static_cast<Eigen::Index>(NY); ++i)
            if(constrained[i] && u_unconstrained[i] != u_sat[i])
                direction[i] = u_unconstrained[i] > u_sat[i] ? std::int8_t{1} : std::int8_t{-1};
        return direction;
    }
};

}

#endif
