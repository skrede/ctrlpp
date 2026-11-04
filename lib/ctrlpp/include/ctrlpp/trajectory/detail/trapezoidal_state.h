#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRAPEZOIDAL_STATE_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DETAIL_TRAPEZOIDAL_STATE_H

#include "ctrlpp/trajectory/detail/trapezoidal_profile.h"

#include "ctrlpp/trajectory/trajectory_error.h"
#include "ctrlpp/trajectory/trapezoidal_types.h"

#include "ctrlpp/expected.h"

#include <cmath>

namespace ctrlpp::detail
{

template <typename Scalar>
struct trapezoidal_rescaled_state
{
    trapezoidal_rescaled_state()
        : v_v{}
        , a_a{}
        , a_d{}
        , T_a{}
        , T_v{}
        , T_d{}
        , T{}
        , triangular{}
    {
    }

    Scalar v_v;
    Scalar a_a;
    Scalar a_d;
    Scalar T_a;
    Scalar T_v;
    Scalar T_d;
    Scalar T;
    bool triangular;
};

template <typename Scalar>
class trapezoidal_state
{
protected:
    Scalar m_q0;
    Scalar m_q1;
    Scalar m_sigma;
    Scalar m_v0;
    Scalar m_v1;
    Scalar m_v_v;
    Scalar m_v_max;
    Scalar m_a;
    Scalar m_a_a;
    Scalar m_a_d;
    Scalar m_T_a;
    Scalar m_T_v;
    Scalar m_T_d;
    Scalar m_T;
    bool m_triangular;
    trapezoidal_disposition<Scalar> m_disposition;

    // @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.13a)-(3.13c), p.71
    template <typename Config>
    explicit trapezoidal_state(Config const& cfg)
        : m_q0{cfg.q0}
        , m_q1{cfg.q1}
        , m_sigma{}
        , m_v0{}
        , m_v1{}
        , m_v_v{}
        , m_v_max{cfg.v_max}
        , m_a{solve_trapezoidal_acceleration<Scalar>(cfg)}
        , m_a_a{}
        , m_a_d{}
        , m_T_a{}
        , m_T_v{}
        , m_T_d{}
        , m_T{}
        , m_triangular{}
        , m_disposition{}
    {
        m_disposition.commanded_acceleration = cfg.a_max;
        m_disposition.realized_acceleration = m_a;
        auto const h = cfg.q1 - cfg.q0;
        m_sigma = positive_frame_sign(h);
        auto const sv0 = m_sigma * cfg.v0;
        auto const sv1 = m_sigma * cfg.v1;
        set_phases(std::abs(h), sv0, sv1, cfg.v_max);
        m_T = m_T_a + m_T_v + m_T_d;
        m_v0 = sv0;
        m_v1 = sv1;
    }

    void set_phases(Scalar abs_h, Scalar sv0, Scalar sv1, Scalar v)
    {
        auto const a = m_a;
        auto const v_tri = triangular_peak(abs_h, sv0, sv1, a);
        m_triangular = v_tri < v;
        m_v_v = m_triangular ? v_tri : v;
        m_a_a = a;
        m_a_d = a;
        if(m_triangular)
        {
            m_T_a = triangular_ramp_duration(abs_h, v_tri, sv0, sv1, a);
            m_T_d = triangular_ramp_duration(abs_h, v_tri, sv1, sv0, a);
            m_T_v = Scalar{0};
            return;
        }
        set_cruise_phases(abs_h, sv0, sv1);
    }

    // The residual is the slack of a h >= v^2 - (v0^2 + v1^2) / 2, which holds wherever the cruise shape is
    // chosen, so the clamp below is a rounding guard and not a stand-in for a rejected command.
    void set_cruise_phases(Scalar abs_h, Scalar sv0, Scalar sv1)
    {
        m_T_a = (m_v_v - sv0) / m_a_a;
        m_T_d = (m_v_v - sv1) / m_a_d;
        auto const d_a = sv0 * m_T_a + Scalar{0.5} * m_a_a * m_T_a * m_T_a;
        auto const d_d = sv1 * m_T_d + Scalar{0.5} * m_a_d * m_T_d * m_T_d;
        m_T_v = (abs_h - d_a - d_d) / m_v_v;
        if(m_T_v < Scalar{0})
        {
            m_T_v = Scalar{0};
        }
    }

    trapezoidal_rescaled_state<Scalar> current_state() const
    {
        trapezoidal_rescaled_state<Scalar> current;
        current.v_v = m_v_v;
        current.a_a = m_a_a;
        current.a_d = m_a_d;
        current.T_a = m_T_a;
        current.T_v = m_T_v;
        current.T_d = m_T_d;
        current.T = m_T;
        current.triangular = m_triangular;
        return current;
    }

    void apply(trapezoidal_rescaled_state<Scalar> const& s)
    {
        m_v_v = s.v_v;
        m_a_a = s.a_a;
        m_a_d = s.a_d;
        m_T_a = s.T_a;
        m_T_v = s.T_v;
        m_T_d = s.T_d;
        m_T = s.T;
        m_triangular = s.triangular;
    }

    // The ramps carry the sign of the velocity change they realize: a valley decelerates away from v0 and
    // accelerates back up to v1.
    ctrlpp::expected<trapezoidal_rescaled_state<Scalar>, trajectory_error> rescaled_state_at(Scalar h, Scalar a, Scalar v) const
    {
        auto const T_a = std::abs(v - m_v0) / a;
        auto const T_d = std::abs(v - m_v1) / a;
        auto const d_a = (m_v0 + v) * T_a / Scalar{2};
        auto const d_d = (m_v1 + v) * T_d / Scalar{2};
        auto const T_v = (h - d_a - d_d) / v;
        if(!(T_a >= Scalar{0}) || !(T_v >= Scalar{0}) || !(T_d >= Scalar{0}))
        {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        trapezoidal_rescaled_state<Scalar> solved;
        solved.v_v = v;
        solved.a_a = (v >= m_v0) ? a : -a;
        solved.a_d = (v >= m_v1) ? a : -a;
        solved.T_a = T_a;
        solved.T_v = T_v;
        solved.T_d = T_d;
        solved.T = T_a + T_v + T_d;
        solved.triangular = (T_v <= Scalar{0});
        return solved;
    }
};

}

#endif
