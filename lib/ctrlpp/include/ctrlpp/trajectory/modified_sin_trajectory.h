#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_MODIFIED_SIN_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_MODIFIED_SIN_TRAJECTORY_H

// @cite biagiotti2009 -- Sec. 3.8, p.124-127

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/util/concepts.h"

#include <cmath>
#include <numbers>
#include <algorithm>

namespace ctrlpp
{

// @cite biagiotti2009 -- Sec. 3.8, p.124-127
template <ctrlpp_floating_scalar Scalar>
class modified_sin_trajectory
{
public:
    struct config
    {
        Scalar q0;
        Scalar q1;
        Scalar T;
    };

    explicit modified_sin_trajectory(config const& cfg)
        : m_q0{cfg.q0}
        , m_h{std::abs(cfg.q1 - cfg.q0)}
        , m_sigma{cfg.q1 >= cfg.q0 ? Scalar{1} : Scalar{-1}}
        , m_T{cfg.T}
    {
    }

    // @cite biagiotti2009 -- Sec. 3.8, p.124-127
    trajectory_point<Scalar, 1> evaluate(Scalar t) const
    {
        auto const tc = std::clamp(t, Scalar{0}, m_T);
        auto const tau = tc / m_T;
        auto const inv_T = Scalar{1} / m_T;

        Scalar q_norm{}; // [0, 1]
        Scalar v_norm{}; // dq/dtau
        Scalar a_norm{}; // d2q/dtau2
        eval_normalized(tau, q_norm, v_norm, a_norm);

        return {
                .position = Vector<Scalar, 1>{m_q0 + m_sigma * m_h * q_norm},
                .velocity = Vector<Scalar, 1>{m_sigma * m_h * inv_T * v_norm},
                .acceleration = Vector<Scalar, 1>{m_sigma * m_h * inv_T * inv_T * a_norm},
        };
    }

    Scalar duration() const { return m_T; }

    // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2009, Sec. 3.8, p.126-127
    Scalar peak_velocity() const
    {
        constexpr auto pi = std::numbers::pi_v<Scalar>;
        return Scalar{4} * pi * m_h / ((pi + Scalar{4}) * m_T);
    }

private:
    Scalar m_q0;
    Scalar m_h;
    Scalar m_sigma;
    Scalar m_T;

    void eval_normalized(Scalar tau, Scalar& q_norm, Scalar& v_norm, Scalar& a_norm) const
    {
        constexpr auto pi = std::numbers::pi_v<Scalar>;
        // @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for Automatic Machines and Robots", 2009, Sec. 3.8, p.126
        auto const A = pi / (pi + Scalar{4});

        if(tau < Scalar{1} / Scalar{8})
        {
            eval_cycloidal_start(tau, A, q_norm, v_norm, a_norm);
        }
        else if(tau < Scalar{7} / Scalar{8})
        {
            eval_sinusoidal_middle(tau, A, q_norm, v_norm, a_norm);
        }
        else
        {
            eval_cycloidal_end(tau, A, q_norm, v_norm, a_norm);
        }
    }

    // @cite biagiotti2009 -- Sec. 3.8, p.126
    void eval_cycloidal_start(Scalar tau, Scalar A, Scalar& q_norm, Scalar& v_norm, Scalar& a_norm) const
    {
        constexpr auto pi = std::numbers::pi_v<Scalar>;
        auto const w = Scalar{4} * pi * tau;
        v_norm = A * (Scalar{1} - std::cos(w));
        q_norm = A * (tau - std::sin(w) / (Scalar{4} * pi));
        a_norm = A * Scalar{4} * pi * std::sin(w);
    }

    // @cite biagiotti2009 -- Sec. 3.8, p.126
    void eval_sinusoidal_middle(Scalar tau, Scalar A, Scalar& q_norm, Scalar& v_norm, Scalar& a_norm) const
    {
        constexpr auto pi = std::numbers::pi_v<Scalar>;
        auto const arg = Scalar{4} * pi * tau / Scalar{3} - pi / Scalar{6};
        auto const q_at_1_8 = A * (Scalar{1} / Scalar{8} - Scalar{1} / (Scalar{4} * pi));

        v_norm = A * (Scalar{1} + Scalar{3} * std::sin(arg));
        a_norm = A * Scalar{4} * pi * std::cos(arg);
        q_norm = q_at_1_8 + A * (tau - Scalar{1} / Scalar{8}) + Scalar{9} * A / (Scalar{4} * pi) * (Scalar{1} - std::cos(arg));
    }

    // The end region mirrors the start region: q(tau) = 1 - q_AB(1 - tau).
    // @cite biagiotti2009 -- Sec. 3.8, p.126
    void eval_cycloidal_end(Scalar tau, Scalar A, Scalar& q_norm, Scalar& v_norm, Scalar& a_norm) const
    {
        constexpr auto pi = std::numbers::pi_v<Scalar>;
        auto const tau_m = Scalar{1} - tau;
        auto const w = Scalar{4} * pi * tau_m;
        auto const v_m = A * (Scalar{1} - std::cos(w));
        auto const q_m = A * (tau_m - std::sin(w) / (Scalar{4} * pi));
        auto const a_m = A * Scalar{4} * pi * std::sin(w);

        q_norm = Scalar{1} - q_m;
        v_norm = v_m;
        a_norm = -a_m;
    }
};

static_assert(trajectory_segment<modified_sin_trajectory<double>, double, 1>);

}

#endif
