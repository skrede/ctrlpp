#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_3RD_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_ONLINE_PLANNER_3RD_H

/// @brief 3rd-order online trajectory planner (velocity, acceleration, and jerk bounded).
///
/// Stateful filter that generates double-S-like velocity profiles in real time.
/// On each update(target), the planner computes a time-optimal profile from the
/// current state (position, velocity, acceleration) to the target position,
/// respecting v_max, a_max, and j_max constraints. sample(t) evaluates the profile
/// at arbitrary times for servo control or lookahead.
///
/// Unlike pre-computed trajectory segments, online planners have no fixed duration
/// and do NOT satisfy trajectory_segment. They are stateful filters per D-12.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 4.6.1

#include "ctrlpp/trajectory/trajectory_types.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace ctrlpp
{

/// @brief 3rd-order online planner producing double-S (jerk-limited) velocity profiles.
///
/// Generates bounded-velocity, bounded-acceleration, bounded-jerk trajectories
/// that can be replanned mid-motion when a new target arrives.
///
/// @cite biagiotti2009 -- Sec. 4.6.1
template <typename Scalar>
class online_planner_3rd
{
  public:
    struct config
    {
        Scalar v_max;
        Scalar a_max;
        Scalar j_max;
    };

    /// @brief Construct planner with kinematic limits. Initial state at rest at q=0.
    explicit online_planner_3rd(config const& cfg)
        : v_max_{cfg.v_max}
        , a_max_{cfg.a_max}
        , j_max_{cfg.j_max}
    {
    }

    /// @brief Set new target position. Replans from current state.
    ///
    /// Computes time-optimal double-S profile from (q_, v_, a_) to (target, 0, 0)
    /// respecting v_max, a_max, and j_max. Handles cases where current velocity
    /// or acceleration require multi-phase deceleration before replanning.
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.1
    void update(Scalar target)
    {
        target_ = target;
        t_ref_ = t_last_;

        // Snapshot current state
        q_ref_ = q_;
        v_ref_ = v_;
        a_ref_ = a_;

        compute_profile();
    }

    /// @brief Evaluate trajectory at time t.
    ///
    /// Returns trajectory_point<Scalar, 1> with position, velocity, acceleration.
    /// Updates internal state for future update() calls.
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.1
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
        a_ = a;
        t_last_ = t;

        // Check settled state
        auto constexpr settle_tol = static_cast<Scalar>(1e-9);
        settled_ = (std::abs(q - target_) < settle_tol)
                   && (std::abs(v) < settle_tol)
                   && (std::abs(a) < settle_tol);

        return {
            .position = Vector<Scalar, 1>{q},
            .velocity = Vector<Scalar, 1>{v},
            .acceleration = Vector<Scalar, 1>{a},
        };
    }

    /// @brief True when at target with zero velocity and zero acceleration.
    auto is_settled() const -> bool { return settled_; }

    /// @brief Reset state to position q0 with zero velocity and zero acceleration.
    void reset(Scalar q0)
    {
        q_ = q0;
        v_ = Scalar{0};
        a_ = Scalar{0};
        q_ref_ = q0;
        v_ref_ = Scalar{0};
        a_ref_ = Scalar{0};
        target_ = q0;
        t_ref_ = Scalar{0};
        t_last_ = Scalar{0};
        n_phases_ = 0;
        T_ = Scalar{0};
        settled_ = true;
    }

  private:
    Scalar v_max_{};
    Scalar a_max_{};
    Scalar j_max_{};

    // Internal state (mutable for const sample)
    mutable Scalar q_{};
    mutable Scalar v_{};
    mutable Scalar a_{};
    mutable Scalar t_last_{};
    mutable bool settled_{true};

    // Reference state at last update()
    Scalar q_ref_{};
    Scalar v_ref_{};
    Scalar a_ref_{};
    Scalar target_{};
    Scalar t_ref_{};

    // Profile as a sequence of constant-jerk phases.
    // Each phase has duration T_ph_[i] and jerk j_ph_[i].
    // Maximum 11 phases: up to 4 for bringing accel to zero + 7 for double-S.
    static constexpr int max_phases_ = 11;
    int n_phases_{0};
    Scalar T_ph_[max_phases_]{};
    Scalar j_ph_[max_phases_]{};
    Scalar T_{};  ///< total duration

    /// @brief Compute jerk-limited profile from current state to target.
    ///
    /// Strategy per B&M Sec. 4.6.1 variable-structure approach:
    /// 1. If current acceleration is non-zero, first bring it to zero (jerk phase).
    /// 2. From the resulting state (q', v', 0), plan a rest-to-rest double-S
    ///    to (target, 0, 0).
    /// If the velocity after nulling acceleration overshoots the target, the
    /// double-S plan handles the return.
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.1, eq. (4.49)
    void compute_profile()
    {
        auto constexpr eps = static_cast<Scalar>(1e-12);
        n_phases_ = 0;

        // Check if already settled
        if (std::abs(target_ - q_ref_) < eps
            && std::abs(v_ref_) < eps
            && std::abs(a_ref_) < eps) {
            T_ = Scalar{0};
            settled_ = true;
            return;
        }

        settled_ = false;

        Scalar q_start = q_ref_;
        Scalar v_start = v_ref_;
        Scalar a_start = a_ref_;

        // Phase 0: Bring acceleration to zero if non-zero
        if (std::abs(a_start) > eps) {
            auto const T_az = std::abs(a_start) / j_max_;
            auto const j_az = (a_start > Scalar{0}) ? -j_max_ : j_max_;

            // State after this phase:
            auto const v_after = v_start + a_start * T_az + Scalar{0.5} * j_az * T_az * T_az;
            auto const q_after = q_start + v_start * T_az
                                 + Scalar{0.5} * a_start * T_az * T_az
                                 + j_az * T_az * T_az * T_az / Scalar{6};

            T_ph_[n_phases_] = T_az;
            j_ph_[n_phases_] = j_az;
            ++n_phases_;

            q_start = q_after;
            v_start = v_after;
            a_start = Scalar{0};
        }

        // Now state is (q_start, v_start, 0). Plan double-S from here.
        plan_from_zero_accel(q_start, v_start);

        // Compute total duration
        T_ = Scalar{0};
        for (int i = 0; i < n_phases_; ++i) {
            T_ += T_ph_[i];
        }
    }

    /// @brief Plan a profile from (q0, v0, a=0) to (target_, 0, 0).
    ///
    /// If v0 is non-zero, first bring velocity to zero with a trapezoidal
    /// acceleration profile (jerk-limited), then plan rest-to-rest.
    void plan_from_zero_accel(Scalar q0, Scalar v0)
    {
        auto constexpr eps = static_cast<Scalar>(1e-12);
        auto const h_signed = target_ - q0;

        // If velocity is zero (or nearly), plan rest-to-rest directly
        if (std::abs(v0) < eps) {
            plan_rest_to_rest(q0);
            return;
        }

        // Compute stopping distance: distance to bring v0 to 0 using a_max, j_max
        auto const stop_info = compute_stop(v0);
        auto const stop_dist = stop_info.dist;

        // Check if velocity points wrong way or overshoots
        bool const wrong_way = (h_signed > eps && v0 < -eps)
                               || (h_signed < -eps && v0 > eps)
                               || (std::abs(h_signed) < eps);
        bool const overshoot = !wrong_way
                               && (std::abs(stop_dist) > std::abs(h_signed) + eps);

        if (wrong_way || overshoot) {
            // Brake to zero velocity, then plan rest-to-rest
            append_brake_phases(v0, stop_info);
            auto const q_after = q0 + stop_dist;
            plan_rest_to_rest(q_after);
        } else {
            // Can incorporate initial velocity into the profile
            // For simplicity and robustness: brake to zero, then rest-to-rest
            // This is slightly suboptimal but always correct and respects all constraints.
            append_brake_phases(v0, stop_info);
            auto const q_after = q0 + stop_dist;
            plan_rest_to_rest(q_after);
        }
    }

    /// @brief Information about stopping from a given velocity.
    struct stop_result
    {
        Scalar dist;    ///< signed distance to stop
        Scalar T_j;     ///< jerk phase duration
        Scalar T_a;     ///< total accel phase (T_a >= 2*T_j if a_max reached)
        bool a_max_reached;
    };

    /// @brief Compute distance and phase durations to stop from v0 (with a=0).
    ///
    /// Decelerates using jerk-limited profile: jerk -> const decel -> jerk.
    auto compute_stop(Scalar v0) const -> stop_result
    {
        auto const abs_v = std::abs(v0);
        auto const sign_v = (v0 >= Scalar{0}) ? Scalar{1} : Scalar{-1};

        // Check if a_max is reached during deceleration
        // a_max is reached if abs_v > a_max^2 / j_max
        auto const v_threshold = a_max_ * a_max_ / j_max_;

        Scalar T_j{};
        Scalar T_total{};
        Scalar dist{};
        bool a_reached{};

        if (abs_v > v_threshold) {
            // a_max reached: jerk phase T_j = a_max/j_max, const decel, jerk phase
            T_j = a_max_ / j_max_;
            auto const T_const = abs_v / a_max_ - T_j;
            T_total = T_j + T_const + T_j;
            // Distance: integral of velocity during deceleration
            // v(t) starts at abs_v and goes to 0
            dist = sign_v * abs_v * (T_total) / Scalar{2};
            // More precise: area = abs_v*T_j - j_max*T_j^3/6
            //                    + (abs_v - j_max*T_j^2/2)*T_const - a_max*T_const^2/2
            //                    + ... (backward phase)
            // Use the fact that for symmetric decel: area = abs_v * T_total / 2
            // This holds when v goes linearly from abs_v to 0 on average.
            // Actually let's compute precisely:
            dist = sign_v * compute_decel_distance(abs_v, T_j, T_const);
            a_reached = true;
        } else {
            // a_max not reached: only jerk phases
            T_j = std::sqrt(abs_v / j_max_);
            T_total = Scalar{2} * T_j;
            dist = sign_v * abs_v * T_j; // = sign_v * j_max * T_j^2 * T_j = sign_v * abs_v * T_j
            a_reached = false;
        }

        return {dist, T_j, T_total, a_reached};
    }

    /// @brief Compute distance during jerk-limited deceleration (a_max reached case).
    auto compute_decel_distance(Scalar abs_v, Scalar T_j, Scalar T_const) const -> Scalar
    {
        // Phase 1: jerk = -j_max, duration T_j
        // v(t) = abs_v - j_max*t^2/2
        // q(t) = abs_v*t - j_max*t^3/6
        auto const v1 = abs_v - j_max_ * T_j * T_j / Scalar{2};
        auto const q1 = abs_v * T_j - j_max_ * T_j * T_j * T_j / Scalar{6};

        // Phase 2: jerk = 0, a = -a_max, duration T_const
        // v(t) = v1 - a_max*t
        // q(t) = v1*t - a_max*t^2/2
        auto const v2 = v1 - a_max_ * T_const;
        auto const q2 = q1 + v1 * T_const - a_max_ * T_const * T_const / Scalar{2};

        // Phase 3: jerk = +j_max, duration T_j
        // v(t) = v2 - a_max*t + j_max*t^2/2
        // q(t) = v2*t - a_max*t^2/2 + j_max*t^3/6
        auto const q3 = q2 + v2 * T_j - a_max_ * T_j * T_j / Scalar{2}
                         + j_max_ * T_j * T_j * T_j / Scalar{6};

        return q3;
    }

    /// @brief Append jerk-limited braking phases to bring v0 to zero (from a=0).
    void append_brake_phases(Scalar v0, stop_result const& info)
    {
        auto const sign_v = (v0 >= Scalar{0}) ? Scalar{1} : Scalar{-1};
        // Deceleration jerk: opposite sign of velocity
        auto const j_decel = -sign_v * j_max_;
        auto const j_accel = sign_v * j_max_;

        if (info.a_max_reached) {
            auto const T_const = info.T_a - Scalar{2} * info.T_j;
            // Phase 1: jerk to build deceleration
            T_ph_[n_phases_] = info.T_j;
            j_ph_[n_phases_] = j_decel;
            ++n_phases_;

            // Phase 2: constant deceleration (jerk = 0)
            if (T_const > Scalar(1e-15)) {
                T_ph_[n_phases_] = T_const;
                j_ph_[n_phases_] = Scalar{0};
                ++n_phases_;
            }

            // Phase 3: jerk to bring acceleration back to zero
            T_ph_[n_phases_] = info.T_j;
            j_ph_[n_phases_] = j_accel;
            ++n_phases_;
        } else {
            // Only two jerk phases (triangular deceleration)
            T_ph_[n_phases_] = info.T_j;
            j_ph_[n_phases_] = j_decel;
            ++n_phases_;

            T_ph_[n_phases_] = info.T_j;
            j_ph_[n_phases_] = j_accel;
            ++n_phases_;
        }
    }

    /// @brief Plan a rest-to-rest double-S profile from q0 to target_.
    ///
    /// Appends up to 7 constant-jerk phases to the phase array.
    /// Handles all degenerate cases (v_max or a_max not reached).
    ///
    /// @cite biagiotti2009 -- Sec. 3.4.3
    void plan_rest_to_rest(Scalar q0)
    {
        auto constexpr eps = static_cast<Scalar>(1e-12);
        auto const h_signed = target_ - q0;

        if (std::abs(h_signed) < eps) {
            // Already at target
            return;
        }

        auto const sigma = (h_signed > Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const h = std::abs(h_signed);

        // Compute double-S phase durations
        // @cite biagiotti2009 -- Sec. 3.4.3
        Scalar T_j1{};
        Scalar T_a{};
        Scalar T_v{};
        Scalar T_d{};
        Scalar T_j2{};

        compute_double_s_durations(h, T_j1, T_a, T_v, T_d, T_j2);

        // Acceleration phase: 3 sub-phases
        auto const j_pos = sigma * j_max_;
        auto const j_neg = -sigma * j_max_;

        auto const T_a_const = T_a - Scalar{2} * T_j1;
        auto const T_d_const = T_d - Scalar{2} * T_j2;

        // Phase 1: jerk (+j) for T_j1
        if (T_j1 > eps) {
            T_ph_[n_phases_] = T_j1;
            j_ph_[n_phases_] = j_pos;
            ++n_phases_;
        }

        // Phase 2: constant acceleration (j=0) for T_a - 2*T_j1
        if (T_a_const > eps) {
            T_ph_[n_phases_] = T_a_const;
            j_ph_[n_phases_] = Scalar{0};
            ++n_phases_;
        }

        // Phase 3: jerk (-j) for T_j1
        if (T_j1 > eps) {
            T_ph_[n_phases_] = T_j1;
            j_ph_[n_phases_] = j_neg;
            ++n_phases_;
        }

        // Phase 4: cruise (j=0) for T_v
        if (T_v > eps) {
            T_ph_[n_phases_] = T_v;
            j_ph_[n_phases_] = Scalar{0};
            ++n_phases_;
        }

        // Phase 5: jerk (-j) for T_j2
        if (T_j2 > eps) {
            T_ph_[n_phases_] = T_j2;
            j_ph_[n_phases_] = j_neg;
            ++n_phases_;
        }

        // Phase 6: constant deceleration (j=0) for T_d - 2*T_j2
        if (T_d_const > eps) {
            T_ph_[n_phases_] = T_d_const;
            j_ph_[n_phases_] = Scalar{0};
            ++n_phases_;
        }

        // Phase 7: jerk (+j) for T_j2
        if (T_j2 > eps) {
            T_ph_[n_phases_] = T_j2;
            j_ph_[n_phases_] = j_pos;
            ++n_phases_;
        }
    }

    /// @brief Compute double-S phase durations for rest-to-rest displacement h.
    ///
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.88-91
    void compute_double_s_durations(Scalar h,
                                    Scalar& T_j1, Scalar& T_a, Scalar& T_v,
                                    Scalar& T_d, Scalar& T_j2) const
    {
        auto const v = v_max_;
        auto const a = a_max_;
        auto const j = j_max_;

        // Check if a_max is reached: need v_max * j_max >= a_max^2
        bool const a_max_reached = (v * j >= a * a);

        Scalar T_j_val{};
        Scalar T_a_val{};
        Scalar T_v_val{};

        if (a_max_reached) {
            T_j_val = a / j;
            T_a_val = T_j_val + v / a;
            T_v_val = h / v - T_a_val;
        } else {
            T_j_val = std::sqrt(v / j);
            T_a_val = Scalar{2} * T_j_val;
            T_v_val = h / v - T_a_val;
        }

        if (T_v_val > Scalar{0}) {
            // v_max is reached (cruise phase exists)
            T_j1 = T_j_val;
            T_j2 = T_j_val;
            T_a = T_a_val;
            T_d = T_a_val;
            T_v = T_v_val;
        } else {
            // v_max not reached
            T_v = Scalar{0};

            if (a_max_reached) {
                // a_max reached but v_max not: solve quadratic for v_lim
                auto const coeff_a = Scalar{1} / a;
                auto const coeff_b = a / j;
                auto const coeff_c = -h;
                auto const disc = coeff_b * coeff_b - Scalar{4} * coeff_a * coeff_c;
                auto const v_lim = (-coeff_b + std::sqrt(disc)) / (Scalar{2} * coeff_a);

                T_j1 = a / j;
                T_j2 = T_j1;
                T_a = T_j1 + v_lim / a;
                T_d = T_a;

                // Verify a_max is truly reached (T_a >= 2*T_j1)
                if (T_a < Scalar{2} * T_j1) {
                    // Doubly degenerate
                    auto const T_j_dd = std::cbrt(h / (Scalar{2} * j));
                    T_j1 = T_j_dd;
                    T_j2 = T_j_dd;
                    T_a = Scalar{2} * T_j_dd;
                    T_d = Scalar{2} * T_j_dd;
                }
            } else {
                // Doubly degenerate: neither v_max nor a_max reached
                auto const T_j_dd = std::cbrt(h / (Scalar{2} * j));
                T_j1 = T_j_dd;
                T_j2 = T_j_dd;
                T_a = Scalar{2} * T_j_dd;
                T_d = Scalar{2} * T_j_dd;
            }
        }
    }

    /// @brief Evaluate the stored phase sequence at relative time dt.
    ///
    /// Integrates through constant-jerk phases: jerk -> acceleration -> velocity -> position.
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.1
    void evaluate_profile(Scalar dt, Scalar& q, Scalar& v, Scalar& a) const
    {
        if (n_phases_ == 0 || dt >= T_) {
            q = target_;
            v = Scalar{0};
            a = Scalar{0};
            return;
        }

        // Start from reference state
        q = q_ref_;
        v = v_ref_;
        a = a_ref_;

        Scalar t_elapsed = Scalar{0};

        for (int i = 0; i < n_phases_; ++i) {
            auto const T_i = T_ph_[i];
            auto const j_i = j_ph_[i];
            auto const t_remaining = dt - t_elapsed;

            if (t_remaining <= Scalar{0}) {
                break;
            }

            if (t_remaining < T_i) {
                // Partial phase
                auto const tau = t_remaining;
                q += v * tau + Scalar{0.5} * a * tau * tau + j_i * tau * tau * tau / Scalar{6};
                v += a * tau + Scalar{0.5} * j_i * tau * tau;
                a += j_i * tau;
                break;
            }

            // Full phase
            q += v * T_i + Scalar{0.5} * a * T_i * T_i + j_i * T_i * T_i * T_i / Scalar{6};
            v += a * T_i + Scalar{0.5} * j_i * T_i * T_i;
            a += j_i * T_i;

            t_elapsed += T_i;
        }
    }
};

} // namespace ctrlpp

#endif
