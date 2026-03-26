#ifndef HPP_GUARD_CTRLPP_TRAJ_MODIFIED_SIN_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJ_MODIFIED_SIN_TRAJECTORY_H

/// @brief Modified sinusoidal velocity profile with harmonic+cycloidal blend.
///
/// Uses a sinusoidal central region [T/8, 7T/8) with cycloidal transition
/// regions at start [0, T/8) and end [7T/8, T]. The profile is symmetric
/// about T/2 and produces very smooth acceleration with matched peak
/// accelerations across all regions.
///
/// Three regions spanning the full profile:
///   Region A-B [0, T/8):     Cycloidal start
///   Region B-C [T/8, 7T/8):  Sinusoidal middle
///   Region C-D [7T/8, T]:    Cycloidal end (symmetric to A-B)
///
/// Normalized velocity profile (tau = t/T):
///   A-B: v(tau) = A*(1 - cos(4*pi*tau))          where A = pi/(pi+4)
///   B-C: v(tau) = A*(1 + 3*sin(4*pi*tau/3 - pi/6))
///   C-D: v(tau) = A*(1 - cos(4*pi*(1-tau)))       (symmetry)
///
/// Key constants: max vel = 4*pi*h/((pi+4)*T), max acc = 4*pi^2*h/((pi+4)*T^2).
///
/// @cite biagiotti2009 -- Sec. 3.8, p.124-127

#include "ctrlpp/traj/trajectory_segment.h"
#include "ctrlpp/traj/trajectory_types.h"

#include <algorithm>
#include <cmath>
#include <numbers>

namespace ctrlpp
{

/// @brief Modified sinusoidal trajectory with harmonic+cycloidal blend.
///
/// @cite biagiotti2009 -- Sec. 3.8, p.124-127
template <typename Scalar>
class modified_sin_trajectory
{
public:
    struct config
    {
        Scalar q0; ///< Start position
        Scalar q1; ///< End position
        Scalar T;  ///< Total duration (user-specified)
    };

    explicit modified_sin_trajectory(config const& cfg)
        : q0_{cfg.q0}
        , h_{std::abs(cfg.q1 - cfg.q0)}
        , sigma_{cfg.q1 >= cfg.q0 ? Scalar{1} : Scalar{-1}}
        , T_{cfg.T}
    {
    }

    /// @brief Evaluate position, velocity, acceleration at time t.
    ///
    /// Three-region structure with cycloidal start/end and sinusoidal middle.
    /// Second half [T/2, T] uses symmetry: q(t) = h - q(T-t).
    ///
    /// @cite biagiotti2009 -- Sec. 3.8, p.124-127
    auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>
    {
        auto const tc = std::clamp(t, Scalar{0}, T_);
        constexpr auto pi = std::numbers::pi_v<Scalar>;

        // Base amplitude: A = pi/(pi+4) (normalized velocity coefficient)
        // @cite biagiotti2009 -- Sec. 3.8, p.125
        auto const A = pi / (pi + Scalar{4});
        auto const tau = tc / T_;
        auto const inv_T = Scalar{1} / T_;

        Scalar q_norm{};   // normalized position [0, 1]
        Scalar v_norm{};   // normalized velocity (dq/dtau)
        Scalar a_norm{};   // normalized acceleration (d2q/dtau2)

        if (tau < Scalar{1} / Scalar{8}) {
            // Region A-B [0, T/8): Cycloidal start
            // @cite biagiotti2009 -- Sec. 3.8, p.126
            //
            // v(tau) = A*(1 - cos(4*pi*tau))
            // q(tau) = A*(tau - sin(4*pi*tau) / (4*pi))
            // a(tau) = A*4*pi*sin(4*pi*tau)
            auto const w = Scalar{4} * pi * tau;
            v_norm = A * (Scalar{1} - std::cos(w));
            q_norm = A * (tau - std::sin(w) / (Scalar{4} * pi));
            a_norm = A * Scalar{4} * pi * std::sin(w);
        }
        else if (tau < Scalar{7} / Scalar{8}) {
            // Region B-C [T/8, 7T/8): Sinusoidal middle
            // @cite biagiotti2009 -- Sec. 3.8, p.126
            //
            // v(tau) = A*(1 + 3*sin(4*pi*tau/3 - pi/6))
            // a(tau) = A*4*pi*cos(4*pi*tau/3 - pi/6)
            // q(tau) = q(1/8) + A*(tau - 1/8) + 9*A/(4*pi)*(1 - cos(4*pi*tau/3 - pi/6))
            auto const arg = Scalar{4} * pi * tau / Scalar{3} - pi / Scalar{6};
            auto const q_at_1_8 =
                A * (Scalar{1} / Scalar{8} - Scalar{1} / (Scalar{4} * pi));

            v_norm = A * (Scalar{1} + Scalar{3} * std::sin(arg));
            a_norm = A * Scalar{4} * pi * std::cos(arg);
            q_norm = q_at_1_8 + A * (tau - Scalar{1} / Scalar{8})
                     + Scalar{9} * A / (Scalar{4} * pi)
                           * (Scalar{1} - std::cos(arg));
        }
        else {
            // Region C-D [7T/8, T]: Cycloidal end (symmetric to A-B)
            // @cite biagiotti2009 -- Sec. 3.8, p.126
            //
            // By symmetry: q(tau) = 1 - q_AB(1-tau), v(tau) = v_AB(1-tau), a(tau) = -a_AB(1-tau)
            auto const tau_m = Scalar{1} - tau;
            auto const w = Scalar{4} * pi * tau_m;
            auto const v_m = A * (Scalar{1} - std::cos(w));
            auto const q_m = A * (tau_m - std::sin(w) / (Scalar{4} * pi));
            auto const a_m = A * Scalar{4} * pi * std::sin(w);

            q_norm = Scalar{1} - q_m;
            v_norm = v_m;
            a_norm = -a_m;
        }

        // Convert from normalized to physical units:
        // q_phys = q0 + sigma * h * q_norm
        // v_phys = sigma * h / T * v_norm
        // a_phys = sigma * h / T^2 * a_norm
        return {
            .position = Vector<Scalar, 1>{q0_ + sigma_ * h_ * q_norm},
            .velocity = Vector<Scalar, 1>{sigma_ * h_ * inv_T * v_norm},
            .acceleration = Vector<Scalar, 1>{sigma_ * h_ * inv_T * inv_T * a_norm},
        };
    }

    auto duration() const -> Scalar { return T_; }

    /// @brief Peak velocity magnitude = 4*pi*h / ((pi+4)*T).
    /// @cite biagiotti2009 -- Sec. 3.8, p.125
    auto peak_velocity() const -> Scalar
    {
        constexpr auto pi = std::numbers::pi_v<Scalar>;
        return Scalar{4} * pi * h_ / ((pi + Scalar{4}) * T_);
    }

private:
    Scalar q0_;
    Scalar h_;
    Scalar sigma_;
    Scalar T_;
};

static_assert(trajectory_segment<modified_sin_trajectory<double>, double, 1>);

} // namespace ctrlpp

#endif
