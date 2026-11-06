#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_MODIFIED_TRAP_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_MODIFIED_TRAP_TRAJECTORY_H

// @cite biagiotti2009 -- Sec. 3.7, eq. (3.49)-(3.51), p.119-122

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/util/concepts.h"

#include <cmath>
#include <numbers>
#include <algorithm>

namespace ctrlpp
{

// @cite biagiotti2009 -- Sec. 3.7, p.119-122
template <ctrlpp_floating_scalar Scalar>
class modified_trap_trajectory
{
public:
    struct config
    {
        Scalar q0;
        Scalar q1;
        Scalar T;
    };

    explicit modified_trap_trajectory(config const& cfg)
        : m_q0{cfg.q0}
        , m_h{std::abs(cfg.q1 - cfg.q0)}
        , m_sigma{cfg.q1 >= cfg.q0 ? Scalar{1} : Scalar{-1}}
        , m_T{cfg.T}
    {
    }

    // @cite biagiotti2009 -- Sec. 3.7, eq. (3.49)-(3.51), p.119-122
    trajectory_point<Scalar, 1> evaluate(Scalar t) const
    {
        auto const tc = std::clamp(t, Scalar{0}, m_T);

        Scalar q{};
        Scalar dq{};
        Scalar ddq{};

        if(tc <= m_T / Scalar{2})
        {
            eval_first_half(tc, q, dq, ddq);
        }
        else
        {
            // @cite biagiotti2009 -- Sec. 5.1
            auto const t_mirror = m_T - tc;
            eval_first_half(t_mirror, q, dq, ddq);
            q = m_h - q;
            ddq = -ddq;
        }

        return {
                .position = Vector<Scalar, 1>{m_q0 + m_sigma * q},
                .velocity = Vector<Scalar, 1>{m_sigma * dq},
                .acceleration = Vector<Scalar, 1>{m_sigma * ddq},
        };
    }

    Scalar duration() const { return m_T; }

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.7, p.123
    Scalar peak_velocity() const { return Scalar{2} * m_h / m_T; }

private:
    Scalar m_q0;
    Scalar m_h;
    Scalar m_sigma;
    Scalar m_T;

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2008, Sec. 3.7, eq. (3.49)-(3.51), p.120-122
    void eval_first_half(Scalar t, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        constexpr auto pi = std::numbers::pi_v<Scalar>;
        auto const hp = m_h / (Scalar{2} + pi);
        auto const T_sq = m_T * m_T;
        auto const T8 = m_T / Scalar{8};
        auto const a_max = hp * Scalar{8} * pi / T_sq;

        if(t < T8)
        {
            // @cite biagiotti2009 -- Sec. 3.7, eq. (3.49), p.120
            auto const w = Scalar{4} * pi * t / m_T;
            q = hp * (Scalar{2} * t / m_T - std::sin(w) / (Scalar{2} * pi));
            dq = hp / m_T * (Scalar{2} - Scalar{2} * std::cos(w));
            ddq = a_max * std::sin(w);
        }
        else if(t < Scalar{3} * T8)
        {
            // q_b and dq_b are the ramp-up's position and velocity at T/8, where sin(w) = 1 and cos(w) = 0.
            // @cite biagiotti2009 -- Sec. 3.7, eq. (3.50), p.121
            auto const dt = t - T8;
            auto const q_b = hp * (Scalar{1} / Scalar{4} - Scalar{1} / (Scalar{2} * pi));
            auto const dq_b = hp * Scalar{2} / m_T;

            ddq = a_max;
            dq = dq_b + a_max * dt;
            q = q_b + dq_b * dt + Scalar{0.5} * a_max * dt * dt;
        }
        else
        {
            // q_c and dq_c are the constant phase's position and velocity at 3T/8, after T/4 at a_max.
            // @cite biagiotti2009 -- Sec. 3.7, eq. (3.51), p.121
            auto const dt_bc = m_T / Scalar{4};
            auto const q_b = hp * (Scalar{1} / Scalar{4} - Scalar{1} / (Scalar{2} * pi));
            auto const dq_b = hp * Scalar{2} / m_T;
            auto const q_c = q_b + dq_b * dt_bc + Scalar{0.5} * a_max * dt_bc * dt_bc;
            auto const dq_c = dq_b + a_max * dt_bc;

            auto const dt = t - Scalar{3} * T8;
            auto const w = Scalar{4} * pi * dt / m_T;

            ddq = a_max * std::cos(w);
            dq = dq_c + a_max * m_T / (Scalar{4} * pi) * std::sin(w);

            auto const k = m_T / (Scalar{4} * pi);
            q = q_c + dq_c * dt + a_max * k * k * (Scalar{1} - std::cos(w));
        }
    }
};

static_assert(trajectory_segment<modified_trap_trajectory<double>, double, 1>);

}

#endif
