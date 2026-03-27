#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_MODIFIED_TRAP_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_MODIFIED_TRAP_TRAJECTORY_H

/// @brief Modified trapezoidal velocity profile with cycloidal acceleration phases.
///
/// Replaces the constant-acceleration ramps of a standard trapezoidal profile with
/// cycloidal (sinusoidal) transition phases, producing continuous acceleration and
/// reduced vibration excitation. The profile has 6 phases with fixed T/8 boundaries.
///
/// The first half [0, T/2] has three regions:
///   Phase A-B [0, T/8):     Cycloidal acceleration ramp-up
///   Phase B-C [T/8, 3T/8):  Constant acceleration (linear velocity)
///   Phase C-D [3T/8, T/2]:  Cycloidal acceleration ramp-down
///
/// The second half [T/2, T] is computed via symmetry: q(t) = h - q(T-t).
///
/// Key constants: max vel = 2h/T, max acc = 2h(2+pi)/(pi*T^2).
///
/// @cite biagiotti2009 -- Sec. 3.7, eq. (3.49)-(3.51), p.119-122

#include "ctrlpp/trajectory/trajectory_segment.h"
#include "ctrlpp/trajectory/trajectory_types.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace ctrlpp
{

/// @brief Modified trapezoidal trajectory with cycloidal acceleration phases.
///
/// @cite biagiotti2009 -- Sec. 3.7, p.119-122
template <typename Scalar>
class modified_trap_trajectory
{
public:
    struct config
    {
        Scalar q0; ///< Start position
        Scalar q1; ///< End position
        Scalar T;  ///< Total duration (user-specified)
    };

    explicit modified_trap_trajectory(config const& cfg)
        : q0_{cfg.q0}
        , h_{std::abs(cfg.q1 - cfg.q0)}
        , sigma_{cfg.q1 >= cfg.q0 ? Scalar{1} : Scalar{-1}}
        , T_{cfg.T}
    {
    }

    /// @brief Evaluate position, velocity, acceleration at time t.
    ///
    /// Uses six-phase structure with symmetry about T/2.
    ///
    /// @cite biagiotti2009 -- Sec. 3.7, eq. (3.49)-(3.51), p.119-122
    auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>
    {
        auto const tc = std::clamp(t, Scalar{0}, T_);

        Scalar q{};
        Scalar dq{};
        Scalar ddq{};

        if (tc <= T_ / Scalar{2}) {
            eval_first_half(tc, q, dq, ddq);
        }
        else {
            // Symmetry: q(t) = h - q(T-t), dq(t) = dq(T-t), ddq(t) = -ddq(T-t)
            // @cite biagiotti2009 -- Sec. 5.1
            auto const t_mirror = T_ - tc;
            eval_first_half(t_mirror, q, dq, ddq);
            q = h_ - q;
            // dq stays same (velocity is symmetric)
            ddq = -ddq;
        }

        return {
            .position = Vector<Scalar, 1>{q0_ + sigma_ * q},
            .velocity = Vector<Scalar, 1>{sigma_ * dq},
            .acceleration = Vector<Scalar, 1>{sigma_ * ddq},
        };
    }

    auto duration() const -> Scalar { return T_; }

    /// @brief Peak velocity magnitude = 2*h/T.
    /// @cite biagiotti2009 -- Sec. 3.7, p.121
    auto peak_velocity() const -> Scalar { return Scalar{2} * h_ / T_; }

private:
    Scalar q0_;
    Scalar h_;
    Scalar sigma_;
    Scalar T_;

    /// @brief Evaluate the first half [0, T/2] of the normalized profile.
    ///
    /// The displacement scale factor is hp = h / (2 + pi).
    /// Phase boundaries at T/8, 3T/8, T/2.
    ///
    /// @cite biagiotti2009 -- Sec. 3.7, eq. (3.49)-(3.51), p.120
    void eval_first_half(Scalar t, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        constexpr auto pi = std::numbers::pi_v<Scalar>;
        auto const hp = h_ / (Scalar{2} + pi);
        auto const T_sq = T_ * T_;
        auto const T8 = T_ / Scalar{8};

        // Max acceleration in the constant-accel region
        // a_max = hp * 8*pi / T^2 = 2*h*(2+pi) / (pi*T^2) ... wait, let me compute:
        // From Phase A-B: ddq = hp * 8*pi/T^2 * sin(4*pi*t/T)
        // Peak of sin is 1, so a_max = hp * 8*pi / T^2
        auto const a_max = hp * Scalar{8} * pi / T_sq;

        if (t < T8) {
            // Phase A-B [0, T/8): Cycloidal acceleration ramp-up
            // @cite biagiotti2009 -- Sec. 3.7, eq. (3.49), p.120
            auto const w = Scalar{4} * pi * t / T_;
            q = hp * (Scalar{2} * t / T_ - std::sin(w) / (Scalar{2} * pi));
            dq = hp / T_ * (Scalar{2} - Scalar{2} * std::cos(w));
            ddq = a_max * std::sin(w);
        }
        else if (t < Scalar{3} * T8) {
            // Phase B-C [T/8, 3T/8): Constant acceleration
            // @cite biagiotti2009 -- Sec. 3.7, eq. (3.50), p.121
            //
            // Values at t = T/8 (end of Phase A-B):
            //   w(T/8) = 4*pi*(T/8)/T = pi/2
            //   q(T/8) = hp*(1/4 - 1/(2*pi))   [sin(pi/2)=1]
            //   dq(T/8) = hp/T * (2 - 2*cos(pi/2)) = hp*2/T
            auto const dt = t - T8;
            auto const q_b = hp * (Scalar{1} / Scalar{4} - Scalar{1} / (Scalar{2} * pi));
            auto const dq_b = hp * Scalar{2} / T_;

            ddq = a_max;
            dq = dq_b + a_max * dt;
            q = q_b + dq_b * dt + Scalar{0.5} * a_max * dt * dt;
        }
        else {
            // Phase C-D [3T/8, T/2]: Cycloidal acceleration ramp-down from a_max to 0
            // @cite biagiotti2009 -- Sec. 3.7, eq. (3.51), p.121
            //
            // Values at t = 3T/8 (end of Phase B-C):
            //   dt_bc = 3T/8 - T/8 = T/4
            //   dq(3T/8) = dq_b + a_max * T/4
            //   q(3T/8) = q_b + dq_b*T/4 + 0.5*a_max*(T/4)^2
            auto const dt_bc = T_ / Scalar{4};
            auto const q_b = hp * (Scalar{1} / Scalar{4} - Scalar{1} / (Scalar{2} * pi));
            auto const dq_b = hp * Scalar{2} / T_;
            auto const q_c = q_b + dq_b * dt_bc + Scalar{0.5} * a_max * dt_bc * dt_bc;
            auto const dq_c = dq_b + a_max * dt_bc;

            auto const dt = t - Scalar{3} * T8;
            auto const w = Scalar{4} * pi * dt / T_;

            // Acceleration: a_max * cos(w), ramps from a_max to 0 over [0, T/8]
            // At dt=0: cos(0)=1 -> a_max. At dt=T/8: cos(pi/2)=0 -> 0.
            ddq = a_max * std::cos(w);

            // Velocity: integral of a_max*cos(4*pi*dt/T) = a_max*T/(4*pi)*sin(w)
            dq = dq_c + a_max * T_ / (Scalar{4} * pi) * std::sin(w);

            // Position: double integral
            //   integral of velocity above = dq_c*dt + a_max*T/(4*pi) * [-T/(4*pi)*cos(w) + T/(4*pi)]
            //   = dq_c*dt + a_max*(T/(4*pi))^2 * (1 - cos(w))
            auto const k = T_ / (Scalar{4} * pi);
            q = q_c + dq_c * dt + a_max * k * k * (Scalar{1} - std::cos(w));
        }
    }
};

static_assert(trajectory_segment<modified_trap_trajectory<double>, double, 1>);

} // namespace ctrlpp

#endif
