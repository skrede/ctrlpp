#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_2ND_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_2ND_H

/// @brief 2nd-order online trajectory planner (velocity and acceleration bounded).
///
/// Stateful filter that generates trapezoidal-like velocity profiles in real time.
/// On each update(target), the planner computes a time-optimal profile from the
/// current state (position, velocity) to the target position, respecting v_max and
/// a_max constraints. sample(t) evaluates the profile at arbitrary times for servo
/// control or lookahead.
///
/// Unlike pre-computed trajectory segments, online planners have no fixed duration
/// and do NOT satisfy trajectory_segment. They are stateful filters per D-12.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 4.6.2

#include "ctrlpp/trajectory/trajectory_types.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace ctrlpp
{

/// @brief 2nd-order online planner producing trapezoidal velocity profiles.
///
/// Generates bounded-velocity, bounded-acceleration trajectories that can be
/// replanned mid-motion when a new target arrives.
///
/// @cite biagiotti2009 -- Sec. 4.6.2
template <typename Scalar>
class online_planner_2nd
{
  public:
    struct config
    {
        Scalar v_max;
        Scalar a_max;
    };

    /// @brief Construct planner with kinematic limits. Initial state at rest at q=0.
    explicit online_planner_2nd(config const& cfg)
        : v_max_{cfg.v_max}
        , a_max_{cfg.a_max}
    {
    }

    /// @brief Set new target position. Replans from current state.
    ///
    /// Computes time-optimal trapezoidal profile from (q_, v_) to (target, 0)
    /// respecting v_max and a_max. Handles the case where current velocity
    /// requires overshoot recovery.
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.2
    void update(Scalar target)
    {
        target_ = target;
        t_ref_ = t_last_;

        // Snapshot current state
        q_ref_ = q_;
        v_ref_ = v_;

        compute_profile();
    }

    /// @brief Evaluate trajectory at time t.
    ///
    /// Returns trajectory_point<Scalar, 1> with position, velocity, acceleration.
    /// Updates internal state for future update() calls.
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.2
    auto sample(Scalar t) const -> trajectory_point<Scalar, 1>
    {
        auto const dt = std::max(t - t_ref_, Scalar{0});

        Scalar q{};
        Scalar v{};
        Scalar a{};

        evaluate_profile(dt, q, v, a);

        // Update mutable state for future update() calls
        q_ = q;
        v_ = v;
        t_last_ = t;

        // Check settled state
        auto constexpr settle_tol = static_cast<Scalar>(1e-9);
        settled_ = (std::abs(q - target_) < settle_tol)
                   && (std::abs(v) < settle_tol);

        return {
            .position = Vector<Scalar, 1>{q},
            .velocity = Vector<Scalar, 1>{v},
            .acceleration = Vector<Scalar, 1>{a},
        };
    }

    /// @brief True when at target with zero velocity.
    auto is_settled() const -> bool { return settled_; }

    /// @brief Reset state to position q0 with zero velocity.
    void reset(Scalar q0)
    {
        q_ = q0;
        v_ = Scalar{0};
        q_ref_ = q0;
        v_ref_ = Scalar{0};
        target_ = q0;
        t_ref_ = Scalar{0};
        t_last_ = Scalar{0};
        T_a_ = Scalar{0};
        T_v_ = Scalar{0};
        T_d_ = Scalar{0};
        T_ = Scalar{0};
        settled_ = true;
    }

  private:
    Scalar v_max_{};
    Scalar a_max_{};

    // Internal state (mutable for const sample)
    mutable Scalar q_{};
    mutable Scalar v_{};
    mutable Scalar t_last_{};
    mutable bool settled_{true};

    // Reference state at last update()
    Scalar q_ref_{};
    Scalar v_ref_{};
    Scalar target_{};
    Scalar t_ref_{};

    // Profile parameters
    Scalar sigma_{1};
    Scalar v_v_{};       ///< cruise velocity (positive frame)
    Scalar a_a_{};       ///< acceleration magnitude
    Scalar a_d_{};       ///< deceleration magnitude
    Scalar T_a_{};       ///< acceleration phase duration
    Scalar T_v_{};       ///< cruise phase duration
    Scalar T_d_{};       ///< deceleration phase duration
    Scalar T_{};         ///< total duration
    Scalar q_stop_{};    ///< position where braking from v_ref ends (for overshoot)
    bool needs_brake_{false};  ///< whether we need to brake first
    Scalar T_brake_{};   ///< braking duration
    Scalar q_after_brake_{};
    Scalar v_after_brake_{};

    /// @brief Compute trapezoidal profile from (q_ref_, v_ref_) to (target_, 0).
    ///
    /// Handles overshoot recovery when current velocity points away from target
    /// or is too large to decelerate in time.
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.2
    void compute_profile()
    {
        auto const h_signed = target_ - q_ref_;

        // Check if already at target with zero velocity
        auto constexpr eps = static_cast<Scalar>(1e-12);
        if (std::abs(h_signed) < eps && std::abs(v_ref_) < eps) {
            sigma_ = Scalar{1};
            v_v_ = Scalar{0};
            a_a_ = Scalar{0};
            a_d_ = Scalar{0};
            T_a_ = Scalar{0};
            T_v_ = Scalar{0};
            T_d_ = Scalar{0};
            T_ = Scalar{0};
            needs_brake_ = false;
            settled_ = true;
            return;
        }

        settled_ = false;

        // Determine direction
        // If we have velocity, we may need to brake first
        auto const v0 = v_ref_;

        // Stopping distance from current velocity
        auto const stop_dist = v0 * v0 / (Scalar{2} * a_max_);

        // Check if velocity is pointing the wrong way or overshooting
        bool const wrong_direction = (h_signed > Scalar{0} && v0 < Scalar{0})
                                     || (h_signed < Scalar{0} && v0 > Scalar{0})
                                     || (std::abs(h_signed) < eps && std::abs(v0) > eps);

        // Check if we would overshoot (can't decelerate in time)
        bool const overshoot = !wrong_direction
                               && (stop_dist > std::abs(h_signed) + eps)
                               && std::abs(v0) > eps;

        if (wrong_direction || overshoot) {
            // Need to brake first: decelerate to zero, then plan from stopped position
            needs_brake_ = true;
            T_brake_ = std::abs(v0) / a_max_;
            q_after_brake_ = q_ref_ + v0 * T_brake_ / Scalar{2};
            v_after_brake_ = Scalar{0};

            // Now plan from braked position to target
            compute_rest_to_rest(q_after_brake_, target_);
        } else {
            // Can plan directly incorporating current velocity
            needs_brake_ = false;
            T_brake_ = Scalar{0};
            compute_with_initial_velocity(q_ref_, v0, target_);
        }
    }

    /// @brief Plan a rest-to-rest trapezoidal profile.
    void compute_rest_to_rest(Scalar q0, Scalar target)
    {
        auto const h = target - q0;
        sigma_ = (h >= Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const abs_h = std::abs(h);

        if (abs_h < std::numeric_limits<Scalar>::epsilon()) {
            v_v_ = Scalar{0};
            a_a_ = Scalar{0};
            a_d_ = Scalar{0};
            T_a_ = Scalar{0};
            T_v_ = Scalar{0};
            T_d_ = Scalar{0};
            T_ = T_brake_;
            return;
        }

        auto const v = v_max_;
        auto const a = a_max_;

        // Check triangular degenerate: h < v^2/a
        auto const v_tri_sq = a * abs_h;
        auto const v_tri = std::sqrt(v_tri_sq);

        if (v_tri < v) {
            // Triangular
            v_v_ = v_tri;
        } else {
            v_v_ = v;
        }

        a_a_ = a;
        a_d_ = a;
        T_a_ = v_v_ / a_a_;
        T_d_ = v_v_ / a_d_;

        if (v_tri >= v) {
            T_v_ = abs_h / v_v_ - T_a_ / Scalar{2} - T_d_ / Scalar{2};
            if (T_v_ < Scalar{0}) {
                T_v_ = Scalar{0};
            }
        } else {
            T_v_ = Scalar{0};
        }

        T_ = T_brake_ + T_a_ + T_v_ + T_d_;
    }

    /// @brief Plan profile with non-zero initial velocity.
    ///
    /// Uses the same direction as h_signed and accounts for v0.
    void compute_with_initial_velocity(Scalar q0, Scalar v0, Scalar target)
    {
        auto const h = target - q0;
        sigma_ = (h >= Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const abs_h = std::abs(h);
        auto const sv0 = sigma_ * v0; // v0 in positive frame (should be >= 0 here)

        auto const a = a_max_;
        auto const v = v_max_;

        // Deceleration phase: from v_v_ to 0
        // Acceleration phase: from sv0 to v_v_

        // Check if cruise velocity can be reached
        // h = sv0*T_a + 0.5*a*T_a^2 + v_v*T_v + v_v*T_d - 0.5*a*T_d^2
        // With T_a = (v_v - sv0)/a, T_d = v_v/a
        // Triangular: v_v = sqrt((2*a*abs_h + sv0^2)/2)

        auto const v_tri_sq = (Scalar{2} * a * abs_h + sv0 * sv0) / Scalar{2};
        auto const v_tri = std::sqrt(v_tri_sq);

        if (v_tri < v) {
            v_v_ = v_tri;
        } else {
            v_v_ = v;
        }

        // Ensure v_v >= sv0 (otherwise we'd need to decelerate first, handled above)
        if (v_v_ < sv0) {
            v_v_ = sv0;
        }

        a_a_ = a;
        a_d_ = a;
        T_a_ = (v_v_ - sv0) / a_a_;
        T_d_ = v_v_ / a_d_;

        if (v_tri >= v) {
            // Cruise phase: abs_h = area under velocity curve
            // area = sv0*T_a + 0.5*a*T_a^2 + v_v*T_v + v_v*T_d - 0.5*a*T_d^2
            auto const area_accel = sv0 * T_a_ + Scalar{0.5} * a * T_a_ * T_a_;
            auto const area_decel = v_v_ * T_d_ - Scalar{0.5} * a * T_d_ * T_d_;
            T_v_ = (abs_h - area_accel - area_decel) / v_v_;
            if (T_v_ < Scalar{0}) {
                T_v_ = Scalar{0};
            }
        } else {
            T_v_ = Scalar{0};
        }

        T_ = T_a_ + T_v_ + T_d_;
    }

    /// @brief Evaluate the stored profile at relative time dt.
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.2
    void evaluate_profile(Scalar dt, Scalar& q, Scalar& v, Scalar& a) const
    {
        if (T_ <= Scalar{0} || dt >= T_) {
            // At or past the end: settled
            q = target_;
            v = Scalar{0};
            a = Scalar{0};
            return;
        }

        if (needs_brake_) {
            if (dt < T_brake_) {
                // Braking phase: decelerate from v_ref_ to 0
                auto const brake_sign = (v_ref_ >= Scalar{0}) ? Scalar{-1} : Scalar{1};
                a = brake_sign * a_max_;
                v = v_ref_ + a * dt;
                q = q_ref_ + v_ref_ * dt + Scalar{0.5} * a * dt * dt;
                return;
            }
            // Past braking: evaluate rest-to-rest from q_after_brake_
            dt -= T_brake_;
            evaluate_rest_to_rest(dt, q_after_brake_, q, v, a);
            return;
        }

        // Profile with initial velocity
        evaluate_with_velocity(dt, q, v, a);
    }

    /// @brief Evaluate rest-to-rest trapezoidal profile at time dt.
    void evaluate_rest_to_rest(Scalar dt, Scalar q0, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const T_total = T_a_ + T_v_ + T_d_;
        if (dt >= T_total) {
            q = target_;
            v = Scalar{0};
            a = Scalar{0};
            return;
        }

        auto const s = sigma_;
        auto const aa = a_a_;

        if (dt < T_a_) {
            // Acceleration phase
            a = s * aa;
            v = s * aa * dt;
            q = q0 + s * Scalar{0.5} * aa * dt * dt;
        } else if (dt < T_a_ + T_v_) {
            // Cruise phase
            auto const dt_c = dt - T_a_;
            auto const q_end_accel = q0 + s * Scalar{0.5} * aa * T_a_ * T_a_;
            a = Scalar{0};
            v = s * v_v_;
            q = q_end_accel + s * v_v_ * dt_c;
        } else {
            // Deceleration phase -- backward time for precision
            auto const dt_end = T_total - dt;
            a = -s * a_d_;
            v = s * a_d_ * dt_end;
            q = target_ - s * Scalar{0.5} * a_d_ * dt_end * dt_end;
        }
    }

    /// @brief Evaluate profile with non-zero initial velocity.
    void evaluate_with_velocity(Scalar dt, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const s = sigma_;
        auto const aa = a_a_;
        auto const T_total = T_a_ + T_v_ + T_d_;

        if (dt >= T_total) {
            q = target_;
            v = Scalar{0};
            a = Scalar{0};
            return;
        }

        if (dt < T_a_) {
            // Acceleration phase
            a = s * aa;
            v = v_ref_ + s * aa * dt;
            q = q_ref_ + v_ref_ * dt + s * Scalar{0.5} * aa * dt * dt;
        } else if (dt < T_a_ + T_v_) {
            // Cruise phase
            auto const dt_c = dt - T_a_;
            auto const q_end_accel = q_ref_ + v_ref_ * T_a_
                                     + s * Scalar{0.5} * aa * T_a_ * T_a_;
            a = Scalar{0};
            v = s * v_v_;
            q = q_end_accel + s * v_v_ * dt_c;
        } else {
            // Deceleration phase -- backward time for precision (Pitfall 1)
            auto const dt_end = T_total - dt;
            a = -s * a_d_;
            v = s * a_d_ * dt_end;
            q = target_ - s * Scalar{0.5} * a_d_ * dt_end * dt_end;
        }
    }
};

} // namespace ctrlpp

#endif
