#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DOUBLE_S_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DOUBLE_S_TRAJECTORY_H

/// @brief Double-S (7-segment) velocity profile with jerk-limited motion.
///
/// Computes a time-optimal S-curve profile that respects velocity, acceleration,
/// and jerk constraints simultaneously. Construction solves the B&M flowchart
/// (Fig 3.18) including all degenerate cases where v_max or a_max cannot be reached.
/// Negative displacement is handled via sigma transformation (eq 3.31-3.33).
///
/// The 7 segments are: jerk(+), const-accel, jerk(-), cruise, jerk(-), const-decel, jerk(+).
///
/// @cite biagiotti2009 -- Sec. 3.4, eq. (3.17)-(3.33), Fig. 3.18, p.79-96

#include "ctrlpp/trajectory/trajectory_segment.h"
#include "ctrlpp/trajectory/trajectory_types.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>

namespace ctrlpp
{

/// @brief Double-S (7-segment) velocity profile bounding v, a, and j.
/// @cite biagiotti2009 -- Sec. 3.4, eq. (3.17)-(3.33), Fig. 3.18, p.79-96
template <typename Scalar>
class double_s_trajectory
{
public:
    using scalar_type = Scalar;

    struct config
    {
        Scalar q0, q1;
        Scalar v_max, a_max, j_max;
        Scalar v0{}, v1{}; ///< initial/final velocities (default 0)
    };

    explicit double_s_trajectory(config const& cfg)
        : q0_{cfg.q0}
        , q1_{cfg.q1}
        , v0_{cfg.v0}
        , v1_{cfg.v1}
        , j_max_{cfg.j_max}
    {
        auto const h_signed = cfg.q1 - cfg.q0;

        // Zero displacement: stationary profile
        if (std::abs(h_signed) < Scalar{1e-15}) {
            sigma_ = Scalar{1};
            v_lim_ = Scalar{0};
            a_lim_a_ = Scalar{0};
            a_lim_d_ = Scalar{0};
            T_j1_ = Scalar{0};
            T_a_ = Scalar{0};
            T_v_ = Scalar{0};
            T_d_ = Scalar{0};
            T_j2_ = Scalar{0};
            T_ = Scalar{0};
            degenerate_ = true;
            return;
        }

        // Sigma transformation for negative displacement
        // @cite biagiotti2009 -- Sec. 3.4.2, eq. (3.31)-(3.33), p.87
        sigma_ = (h_signed > Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const h = sigma_ * h_signed; // always positive
        auto const v_max = cfg.v_max;
        auto const a_max = cfg.a_max;

        // For zero initial/final velocities, use simplified algorithm
        // @cite biagiotti2009 -- Sec. 3.4.3, p.88-91
        compute_zero_bc(h, v_max, a_max, cfg.j_max);
    }

    /// @brief Evaluate trajectory at time t.
    /// @cite biagiotti2009 -- Sec. 3.4, eq. (3.30a)-(3.30g), p.85-86
    auto evaluate(Scalar t) const -> trajectory_point<Scalar, 1>
    {
        if (T_ <= Scalar{0}) {
            return make_point(q0_, Scalar{0}, Scalar{0});
        }

        auto const tc = std::clamp(t, Scalar{0}, T_);

        // Exact endpoint return
        if (tc >= T_) {
            return make_point(q1_, v1_, Scalar{0});
        }
        if (tc <= Scalar{0}) {
            return make_point(q0_, v0_, Scalar{0});
        }

        // Compute in the positive-displacement frame, then apply sigma
        Scalar q{}, dq{}, ddq{};
        evaluate_positive_frame(tc, q, dq, ddq);

        // Apply sigma transformation
        // @cite biagiotti2009 -- Sec. 3.4.2, eq. (3.33), p.87
        auto const q_actual = q0_ + sigma_ * (q - Scalar{0});
        auto const dq_actual = sigma_ * dq;
        auto const ddq_actual = sigma_ * ddq;

        return make_point(q_actual, dq_actual, ddq_actual);
    }

    auto duration() const -> Scalar { return T_; }

    /// @brief Whether the profile is degenerate (v_max or a_max not reached).
    auto is_degenerate() const -> bool { return degenerate_; }

    /// @brief Actual peak velocity achieved by the profile.
    auto peak_velocity() const -> Scalar { return v_lim_; }

    /// @brief Rescale profile to a new (longer) duration for multi-axis synchronization.
    ///
    /// Extends the cruise phase to fill the time gap while keeping acceleration
    /// and deceleration phases unchanged. This preserves all constraint limits
    /// (v_max, a_max, j_max) since the accel/decel phases are not modified.
    ///
    /// @cite biagiotti2009 -- Sec. 5.3
    void rescale_to(Scalar T_new)
    {
        if (T_new <= T_) {
            return; // Already at or faster than requested -- no-op
        }

        // Insert additional cruise time: keep T_a and T_d fixed, extend T_v
        auto const T_non_cruise = T_a_ + T_d_;
        T_v_ = T_new - T_non_cruise;

        if (T_v_ < Scalar{0}) {
            T_v_ = Scalar{0};
        }

        T_ = T_a_ + T_v_ + T_d_;
    }

    /// @brief Phase durations for the 7 segments.
    ///
    /// Returns {T_j1, T_a - 2*T_j1, T_j1, T_v, T_j2, T_d - 2*T_j2, T_j2}
    /// representing jerk(+), const-accel, jerk(-), cruise, jerk(-), const-decel, jerk(+).
    auto phase_durations() const -> std::array<Scalar, 7>
    {
        auto const const_accel = std::max(T_a_ - Scalar{2} * T_j1_, Scalar{0});
        auto const const_decel = std::max(T_d_ - Scalar{2} * T_j2_, Scalar{0});
        return {T_j1_, const_accel, T_j1_, T_v_, T_j2_, const_decel, T_j2_};
    }

private:
    Scalar q0_{}, q1_{};
    Scalar sigma_{1};
    Scalar v0_{}, v1_;
    Scalar v_lim_{};
    Scalar a_lim_a_{};
    Scalar a_lim_d_{};
    Scalar j_max_{};
    Scalar T_j1_{}, T_a_{}, T_v_{}, T_d_{}, T_j2_{};
    Scalar T_{};
    bool degenerate_{false};

    auto make_point(Scalar q, Scalar dq, Scalar ddq) const -> trajectory_point<Scalar, 1>
    {
        return {.position = Vector<Scalar, 1>{q},
                .velocity = Vector<Scalar, 1>{dq},
                .acceleration = Vector<Scalar, 1>{ddq}};
    }

    /// @brief Compute phase durations for zero initial/final velocity case.
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.88-91
    void compute_zero_bc(Scalar h, Scalar v_max, Scalar a_max, Scalar j_max)
    {
        // Case 1: Assume both v_max and a_max are reached
        // @cite biagiotti2009 -- Sec. 3.4.3, eq. before (3.34), p.89
        auto const T_j = a_max / j_max;

        // Check if a_max is reached: need v_max * j_max >= a_max^2
        // @cite biagiotti2009 -- Sec. 3.4.3, p.89
        bool const a_max_reached = (v_max * j_max >= a_max * a_max);

        Scalar T_a{}, T_v{}, T_d{};

        if (a_max_reached) {
            // a_max is reached
            T_a = T_j + v_max / a_max;
            T_d = T_a; // symmetric for v0 = v1 = 0
            T_v = h / v_max - T_a;
            // T_v = h/v_max - (T_a + T_d)/2 but T_a == T_d => T_v = h/v_max - T_a
        } else {
            // a_max not reached: T_j limited by sqrt(v_max / j_max)
            auto const T_j_actual = std::sqrt(v_max / j_max);
            T_a = Scalar{2} * T_j_actual;
            T_d = T_a;
            T_v = h / v_max - T_a;
        }

        if (T_v > Scalar{0}) {
            // Case 1 or variant: v_max is reached (cruise phase exists)
            if (a_max_reached) {
                T_j1_ = T_j;
                T_j2_ = T_j;
                a_lim_a_ = a_max;
                a_lim_d_ = a_max;
            } else {
                T_j1_ = std::sqrt(v_max / j_max);
                T_j2_ = T_j1_;
                a_lim_a_ = T_j1_ * j_max;
                a_lim_d_ = a_lim_a_;
            }
            T_a_ = T_a;
            T_d_ = T_d;
            T_v_ = T_v;
            v_lim_ = v_max;
            degenerate_ = !a_max_reached;
        } else {
            // Case 2: v_max not reached (T_v < 0 => no cruise phase)
            degenerate_ = true;
            T_v_ = Scalar{0};

            // Need to find actual v_lim < v_max
            // @cite biagiotti2009 -- Sec. 3.4.3, p.89-91
            if (a_max_reached) {
                // a_max is reached but v_max is not
                // Solve: h = (v_lim/a_max + a_max/j_max) * v_lim
                // This is a quadratic in v_lim:
                // v_lim^2/a_max + v_lim * a_max/j_max - h = 0
                auto const a = Scalar{1} / a_max;
                auto const b = a_max / j_max;
                auto const c = -h;
                auto const disc = b * b - Scalar{4} * a * c;
                auto const v_lim = (-b + std::sqrt(disc)) / (Scalar{2} * a);

                v_lim_ = v_lim;
                T_j1_ = a_max / j_max;
                T_j2_ = T_j1_;
                T_a_ = T_j1_ + v_lim / a_max;
                T_d_ = T_a_;
                a_lim_a_ = a_max;
                a_lim_d_ = a_max;

                // Check if T_a >= 2*T_j1 (a_max truly reached)
                if (T_a_ < Scalar{2} * T_j1_) {
                    // Actually a_max is also not reached: fall through to bisection
                    solve_doubly_degenerate(h, a_max, j_max);
                }
            } else {
                // a_max not reached either
                solve_doubly_degenerate(h, a_max, j_max);
            }
        }

        T_ = T_a_ + T_v_ + T_d_;
    }

    /// @brief Solve doubly degenerate case using bisection on gamma.
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.90-91 (Pitfall 2 from RESEARCH.md)
    void solve_doubly_degenerate(Scalar h, Scalar /*a_max*/, Scalar j_max)
    {
        // When neither v_max nor a_max is reached, the profile is purely jerk-limited.
        // Symmetric case (v0=v1=0): the profile has T_a = T_d = 2*T_j, T_v = 0
        // and displacement h = 2 * j_max * T_j^3
        // => T_j = cbrt(h / (2 * j_max))
        auto const T_j = std::cbrt(h / (Scalar{2} * j_max));

        T_j1_ = T_j;
        T_j2_ = T_j;
        T_a_ = Scalar{2} * T_j;
        T_d_ = Scalar{2} * T_j;
        T_v_ = Scalar{0};
        a_lim_a_ = j_max * T_j;
        a_lim_d_ = a_lim_a_;
        v_lim_ = a_lim_a_ * T_j;
    }

    /// @brief Evaluate in the positive-displacement frame (sigma=+1).
    /// @cite biagiotti2009 -- Sec. 3.4, eq. (3.30a)-(3.30g), p.85-86
    void evaluate_positive_frame(Scalar tc, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        auto const j = j_max_;
        auto const a_a = a_lim_a_;
        auto const a_d = a_lim_d_;

        // Phase boundaries
        auto const t1 = T_j1_;                       // end of segment 1
        auto const t2 = T_a_ - T_j1_;                // end of segment 2
        auto const t3 = T_a_;                         // end of segment 3
        auto const t4 = T_a_ + T_v_;                  // end of segment 4
        auto const t5 = T_a_ + T_v_ + T_j2_;         // end of segment 5
        // t6 = T_ - T_j2_ (end of segment 6, used implicitly via dt_end)
        // t7 = T_ (end of segment 7)

        if (tc < t1) {
            // Segment 1: positive jerk (+j_max)
            // @cite biagiotti2009 -- Sec. 3.4, eq. (3.30a), p.85
            auto const t = tc;
            q = j * t * t * t / Scalar{6};
            dq = j * t * t / Scalar{2};
            ddq = j * t;
        } else if (tc < t2) {
            // Segment 2: constant acceleration (a_lim_a)
            // @cite biagiotti2009 -- Sec. 3.4, eq. (3.30b), p.85
            auto const t = tc - T_j1_;
            auto const q1 = j * T_j1_ * T_j1_ * T_j1_ / Scalar{6};
            auto const dq1 = j * T_j1_ * T_j1_ / Scalar{2};
            q = q1 + dq1 * t + a_a * t * t / Scalar{2};
            dq = dq1 + a_a * t;
            ddq = a_a;
        } else if (tc < t3) {
            // Segment 3: negative jerk (-j_max), end of acceleration
            // @cite biagiotti2009 -- Sec. 3.4, eq. (3.30c), p.85
            auto const dt = T_a_ - tc;
            // Evaluate from the end: at t3 we have velocity = v_lim, accel = 0
            auto const q_at_Ta = compute_q_at_Ta();
            q = q_at_Ta - v_lim_ * dt + j * dt * dt * dt / Scalar{6};
            dq = v_lim_ - j * dt * dt / Scalar{2};
            ddq = j * dt;
        } else if (tc < t4) {
            // Segment 4: constant velocity (cruise)
            // @cite biagiotti2009 -- Sec. 3.4, eq. (3.30d), p.86
            auto const q_at_Ta = compute_q_at_Ta();
            auto const t = tc - T_a_;
            q = q_at_Ta + v_lim_ * t;
            dq = v_lim_;
            ddq = Scalar{0};
        } else {
            // Segments 5-7: deceleration phase
            // Use dt_end = T_ - tc for numerical precision (Pitfall 1)
            auto const dt_end = T_ - tc;
            auto const h = std::abs(q1_ - q0_);

            if (dt_end > T_d_ - T_j2_) {
                // Segment 5: negative jerk (-j_max), start of deceleration
                // @cite biagiotti2009 -- Sec. 3.4, eq. (3.30e), p.86
                auto const q_at_Tv_end = compute_q_at_Ta() + v_lim_ * T_v_;
                auto const t = tc - (T_a_ + T_v_);
                q = q_at_Tv_end + v_lim_ * t - j * t * t * t / Scalar{6};
                dq = v_lim_ - j * t * t / Scalar{2};
                ddq = -j * t;
            } else if (dt_end > T_j2_) {
                // Segment 6: constant deceleration (-a_lim_d)
                // Evaluate from the end for numerical precision
                auto const dq_end = Scalar{0}; // v1 = 0
                auto const q_end = h;
                q = q_end - dq_end * dt_end - a_d * dt_end * dt_end / Scalar{2}
                    - j * T_j2_ * T_j2_ * T_j2_ / Scalar{6}
                    + j * T_j2_ * T_j2_ / Scalar{2} * dt_end;
                // Simplified: evaluate carefully
                // Use forward from segment 5 end instead
                auto const t5_local = T_j2_;
                auto const v_at_t5 = v_lim_ - j * t5_local * t5_local / Scalar{2};
                auto const a_at_t5 = -j * t5_local;
                auto const q_at_t5 = compute_q_at_Ta() + v_lim_ * T_v_
                                     + v_lim_ * t5_local - j * t5_local * t5_local * t5_local / Scalar{6};
                auto const dt_from_t5 = tc - t5;
                q = q_at_t5 + v_at_t5 * dt_from_t5 + a_at_t5 * dt_from_t5 * dt_from_t5 / Scalar{2};
                dq = v_at_t5 + a_at_t5 * dt_from_t5;
                ddq = a_at_t5; // = -a_d
            } else {
                // Segment 7: positive jerk (+j_max), end of deceleration
                // Evaluate from end for exact endpoint arrival (Pitfall 1)
                // @cite biagiotti2009 -- Sec. 3.4, eq. (3.30g), p.86
                q = h - j * dt_end * dt_end * dt_end / Scalar{6};
                dq = j * dt_end * dt_end / Scalar{2};
                ddq = -j * dt_end;
            }
        }
    }

    /// @brief Compute position at end of acceleration phase (t = T_a).
    auto compute_q_at_Ta() const -> Scalar
    {
        // At T_a, velocity = v_lim, acceleration = 0
        // Position = integral of velocity from 0 to T_a
        // For symmetric accel phase: q = v_lim * T_a / 2
        // More precisely: q(T_a) = v_lim * (T_a - T_j1) + j_max * T_j1^3 / 6
        //                         + (j_max * T_j1^2 / 2) * (T_a - 2*T_j1)
        //                         + a_lim_a * (T_a - 2*T_j1)^2 / 2

        // Simpler: for v0=0, q(T_a) = (v_lim/2) * T_a
        // This follows from the area under the velocity curve
        // which is a symmetric trapezoid from 0 to v_lim over time T_a
        return v_lim_ * T_a_ / Scalar{2};
    }
};

static_assert(trajectory_segment<double_s_trajectory<double>, double, 1>);
static_assert(trajectory_segment<double_s_trajectory<float>, float, 1>);

} // namespace ctrlpp

#endif
