#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_TRAPEZOIDAL_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_TRAPEZOIDAL_TRAJECTORY_H

/// @brief Trapezoidal velocity profile (LSPB) with degenerate case handling.
///
/// Computes a three-phase (acceleration, cruise, deceleration) trajectory from
/// kinematic limits v_max, a_max. Handles triangular degenerate case when
/// displacement is too short for full trapezoidal, and non-null initial/final
/// velocities with feasibility adjustment per B&M eq. (3.14)-(3.15).
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 3.2, eq. (3.9)-(3.16), p.65-73

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

/// @brief Trapezoidal (LSPB) velocity profile with introspection.
///
/// Construction solves phase durations from v_max, a_max constraints.
/// Degenerate triangular case handled silently.
///
/// @cite biagiotti2009 -- Sec. 3.2, eq. (3.9)-(3.16), p.65-73
template <ctrlpp_floating_scalar Scalar>
class trapezoidal_trajectory
{
  public:
    using scalar_type = Scalar;

    struct config
    {
        Scalar q0, q1;
        Scalar v_max, a_max;
        Scalar v0{}, v1{};
    };

    /// @brief Construct trapezoidal profile from kinematic limits.
    ///
    /// @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.13a)-(3.13c), p.71
    explicit trapezoidal_trajectory(config const& cfg)
        : q0_{cfg.q0}
        , q1_{cfg.q1}
    {
        auto const h = cfg.q1 - cfg.q0;
        sigma_ = (h >= Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const abs_h = std::abs(h);

        // Transform velocities into the positive-displacement frame
        auto const sv0 = sigma_ * cfg.v0;
        auto const sv1 = sigma_ * cfg.v1;
        auto a = cfg.a_max;

        // Feasibility check for non-null BCs -- B&M eq. (3.14)
        // @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.14), p.72
        auto const v_diff_sq = std::abs(sv0 * sv0 - sv1 * sv1) / Scalar{2};
        if (a * abs_h < v_diff_sq) {
            // Infeasible: compute a_lim per eq. (3.15)
            // @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.15), p.72
            a = v_diff_sq / abs_h + std::numeric_limits<Scalar>::epsilon();
        }

        auto v = cfg.v_max;

        // Compute cruise velocity and phase durations
        // @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.13a)-(3.13c), p.71
        // T_a = (v_v - v0) / a, T_d = (v_v - v1) / a; the cruise duration is the
        // residual displacement (after the accel and decel distances) divided by
        // the cruise velocity, which is exact for nonzero boundary velocities.

        // Check triangular degenerate case
        // When v_max cannot be reached: v_v = sqrt((2*a*h + v0^2 + v1^2) / 2)
        auto const v_tri_sq = (Scalar{2} * a * abs_h + sv0 * sv0 + sv1 * sv1) / Scalar{2};
        auto const v_tri = std::sqrt(v_tri_sq);

        if (v_tri < v) {
            // Triangular: cruise velocity limited by displacement
            v_v_ = v_tri;
            triangular_ = true;
        } else {
            v_v_ = v;
            triangular_ = false;
        }

        // Acceleration and deceleration rates (symmetric a for now)
        a_a_ = a;
        a_d_ = a;

        // Phase durations
        T_a_ = (v_v_ - sv0) / a_a_;
        T_d_ = (v_v_ - sv1) / a_d_;

        if (triangular_) {
            T_v_ = Scalar{0};
        } else {
            // Cruise duration from the residual displacement: the accel and decel
            // phases cover d_a = v0*T_a + a_a*T_a^2/2 and d_d = v1*T_d + a_d*T_d^2/2,
            // so the cruise phase covers (abs_h - d_a - d_d) at v_v. This keeps the
            // position continuous at the cruise-to-decel boundary for nonzero v0/v1.
            auto const d_a = sv0 * T_a_ + Scalar{0.5} * a_a_ * T_a_ * T_a_;
            auto const d_d = sv1 * T_d_ + Scalar{0.5} * a_d_ * T_d_ * T_d_;
            T_v_ = (abs_h - d_a - d_d) / v_v_;
            if (T_v_ < Scalar{0}) {
                T_v_ = Scalar{0};
            }
        }

        T_ = T_a_ + T_v_ + T_d_;
        v0_ = sv0;
        v1_ = sv1;
    }

    /// @brief Evaluate trajectory at time t, clamped to [0, T].
    ///
    /// Three-phase branching: acceleration, cruise, deceleration.
    /// Deceleration uses backward time (T - t) for numerical precision.
    ///
    /// @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.13a)-(3.13c), p.71
    auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>
    {
        auto const tc = std::clamp(t, Scalar{0}, T_);
        Scalar q{}, dq{}, ddq{};

        if (tc <= T_a_ && T_a_ > Scalar{0}) {
            // Acceleration phase
            // @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.13a), p.71
            ddq = a_a_;
            dq = v0_ + a_a_ * tc;
            q = v0_ * tc + Scalar{0.5} * a_a_ * tc * tc;
        } else if (tc < T_a_ + T_v_) {
            // Constant velocity (cruise) phase
            auto const dt = tc - T_a_;
            ddq = Scalar{0};
            dq = v_v_;
            // Position at end of accel phase + cruise distance
            q = v0_ * T_a_ + Scalar{0.5} * a_a_ * T_a_ * T_a_ + v_v_ * dt;
        } else {
            // Deceleration phase -- use backward time for precision (Pitfall 1)
            auto const dt_end = T_ - tc;
            ddq = -a_d_;
            dq = v1_ + a_d_ * dt_end;
            // Position from end: q1 - backward integration
            auto const abs_h = std::abs(q1_ - q0_);
            q = abs_h - v1_ * dt_end - Scalar{0.5} * a_d_ * dt_end * dt_end;
        }

        // Apply sigma transformation for sign
        return {
            .position = Vector<Scalar, 1>{q0_ + sigma_ * q},
            .velocity = Vector<Scalar, 1>{sigma_ * dq},
            .acceleration = Vector<Scalar, 1>{sigma_ * ddq},
        };
    }

    /// @brief Total trajectory duration [s].
    auto duration() const -> Scalar { return T_; }

    /// @brief True if the profile is triangular (cruise phase duration is zero).
    ///
    /// @cite biagiotti2009 -- Sec. 3.2.5, p.69
    auto is_triangular() const -> bool { return triangular_; }

    /// @brief Peak velocity (signed, in original frame).
    auto peak_velocity() const -> Scalar { return sigma_ * v_v_; }

    /// @brief Phase durations {T_accel, T_cruise, T_decel}.
    auto phase_durations() const -> std::array<Scalar, 3> { return {T_a_, T_v_, T_d_}; }

    /// @brief Rescale profile to a new (longer) duration for multi-axis synchronization.
    ///
    /// Recomputes cruise velocity and phase durations to achieve T_new while
    /// maintaining the same start/end positions and constraint satisfaction.
    /// Only slowing down is valid (T_new >= current duration).
    ///
    /// @cite biagiotti2009 -- Sec. 5.3
    void rescale_to(Scalar T_new)
    {
        if (T_new <= T_) {
            return; // Already at or faster than requested -- no-op
        }

        auto const abs_h = std::abs(q1_ - q0_);
        if (abs_h < std::numeric_limits<Scalar>::epsilon()) {
            T_ = T_new;
            return;
        }

        // Solve for the new cruise velocity v that realizes total time T_new with the
        // given boundary velocities. With T_a = (v - v0)/a, T_d = (v - v1)/a and the
        // residual-displacement cruise duration T_v = (abs_h - d_a - d_d)/v, the total
        // time T_new = T_a + T_v + T_d reduces to the quadratic
        //   v^2 - [(v0 + v1) + a*T_new] * v + [a*abs_h + (v0^2 + v1^2)/2] = 0.
        // The smaller root is the trapezoidal cruise velocity (the larger root exceeds
        // v_max). This collapses to v^2 - a*T_new*v + a*abs_h = 0 when v0 = v1 = 0.
        auto const a = a_a_; // Use existing acceleration limit
        auto const b_quad = (v0_ + v1_) + a * T_new;
        auto const c_quad = a * abs_h + (v0_ * v0_ + v1_ * v1_) / Scalar{2};
        auto const disc = b_quad * b_quad - Scalar{4} * c_quad;

        Scalar v_cruise{};
        if (disc < Scalar{0}) {
            // Should not happen if T_new >= T_, but handle gracefully with the
            // triangular (unreachable-cruise) velocity for nonzero boundary velocities.
            v_cruise = std::sqrt((Scalar{2} * a * abs_h + v0_ * v0_ + v1_ * v1_) / Scalar{2});
        } else {
            v_cruise = (b_quad - std::sqrt(disc)) / Scalar{2};
        }

        // Clamp v_cruise to avoid negative or zero
        v_cruise = std::max(v_cruise, std::numeric_limits<Scalar>::epsilon());

        // Recompute phase durations with the residual-displacement cruise duration.
        v_v_ = v_cruise;
        T_a_ = (v_cruise - v0_) / a;
        T_d_ = (v_cruise - v1_) / a;
        auto const d_a = v0_ * T_a_ + Scalar{0.5} * a * T_a_ * T_a_;
        auto const d_d = v1_ * T_d_ + Scalar{0.5} * a * T_d_ * T_d_;
        T_v_ = (abs_h - d_a - d_d) / v_cruise;

        if (T_v_ < Scalar{0}) {
            T_v_ = Scalar{0};
        }

        T_ = T_a_ + T_v_ + T_d_;
        triangular_ = (T_v_ <= Scalar{0});
    }

  private:
    Scalar q0_{};
    Scalar q1_{};
    Scalar sigma_{};
    Scalar v0_{};
    Scalar v1_{};
    Scalar v_v_{};
    Scalar a_a_{};
    Scalar a_d_{};
    Scalar T_a_{};
    Scalar T_v_{};
    Scalar T_d_{};
    Scalar T_{};
    bool triangular_{};
};

static_assert(trajectory_segment<trapezoidal_trajectory<double>, double, 1>);

}

#endif
