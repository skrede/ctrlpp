#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_STATE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_DOUBLE_S_STATE_H

#include "ctrlpp/trajectory/detail/boundary_velocity.h"
#include "ctrlpp/trajectory/detail/jerk_limited_ramp.h"

#include "ctrlpp/trajectory/trajectory_error.h"

#include "ctrlpp/expected.h"

#include <cmath>
#include <algorithm>

namespace ctrlpp::detail
{

template <typename Scalar>
class double_s_state
{
protected:
    Scalar m_q0;
    Scalar m_q1;
    Scalar m_sigma;
    Scalar m_v0;
    Scalar m_v1;
    Scalar m_pv0; // positive-displacement frame
    Scalar m_pv1; // positive-displacement frame
    Scalar m_v_lim;
    Scalar m_a_lim_a;
    Scalar m_a_lim_d;
    Scalar m_v_max;
    Scalar m_a_max;
    Scalar m_j_max;
    Scalar m_accel_jerk_duration;
    Scalar m_accel_duration;
    Scalar m_cruise_duration;
    Scalar m_decel_duration;
    Scalar m_decel_jerk_duration;
    Scalar m_duration;
    bool m_degenerate;

    template <typename Config>
    explicit double_s_state(Config const& cfg)
        : m_q0{cfg.q0}
        , m_q1{cfg.q1}
        , m_sigma{1}
        , m_v0{cfg.v0}
        , m_v1{cfg.v1}
        , m_pv0{}
        , m_pv1{}
        , m_v_lim{}
        , m_a_lim_a{}
        , m_a_lim_d{}
        , m_v_max{cfg.v_max}
        , m_a_max{cfg.a_max}
        , m_j_max{cfg.j_max}
        , m_accel_jerk_duration{}
        , m_accel_duration{}
        , m_cruise_duration{}
        , m_decel_duration{}
        , m_decel_jerk_duration{}
        , m_duration{}
        , m_degenerate{false}
    {
    }

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.4, p.79
    template <typename Config>
    static ctrlpp::expected<void, trajectory_error> check_limits(Config const& requested)
    {
        if(!std::isfinite(requested.q0) || !std::isfinite(requested.q1) || !std::isfinite(requested.v0) || !std::isfinite(requested.v1))
        {
            return ctrlpp::unexpected(trajectory_error::non_finite_input);
        }
        if(!std::isfinite(requested.v_max) || requested.v_max <= Scalar{0})
        {
            return ctrlpp::unexpected(trajectory_error::non_positive_velocity_limit);
        }
        if(!std::isfinite(requested.a_max) || requested.a_max <= Scalar{0})
        {
            return ctrlpp::unexpected(trajectory_error::non_positive_acceleration_limit);
        }
        if(!std::isfinite(requested.j_max) || requested.j_max <= Scalar{0})
        {
            return ctrlpp::unexpected(trajectory_error::non_positive_jerk_limit);
        }
        return {};
    }

    // No profile sweeps less than the fastest transition between the boundary velocities, and a zero
    // displacement is the standstill alone.
    template <typename Config>
    static bool is_reachable(Config const& cfg)
    {
        auto const h_signed = cfg.q1 - cfg.q0;
        if(h_signed == Scalar{0})
        {
            return cfg.v0 == Scalar{0} && cfg.v1 == Scalar{0};
        }
        auto const sigma = (h_signed > Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const pv0 = sigma * cfg.v0;
        auto const pv1 = sigma * cfg.v1;
        auto const h_min = no_cruise_displacement(Scalar{0}, std::max(pv0, pv1), std::min(pv0, pv1), cfg.a_max, cfg.j_max);
        return !(sigma * h_signed < h_min);
    }

    void set_ramps(velocity_ramp<Scalar> const& accel, velocity_ramp<Scalar> const& decel)
    {
        m_accel_jerk_duration = accel.T_j;
        m_accel_duration = accel.T;
        m_a_lim_a = accel.a_lim;
        m_decel_jerk_duration = decel.T_j;
        m_decel_duration = decel.T;
        m_a_lim_d = decel.a_lim;
    }
};

}

#endif
