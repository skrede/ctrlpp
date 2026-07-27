#ifndef HPP_GUARD_CTRLPP_TRAJECTORY_DOUBLE_S_TRAJECTORY_H
#define HPP_GUARD_CTRLPP_TRAJECTORY_DOUBLE_S_TRAJECTORY_H

/// @brief Double-S (7-segment) velocity profile with jerk-limited motion.
///
/// Computes a time-optimal S-curve profile that respects velocity, acceleration,
/// and jerk constraints simultaneously. Construction solves the B&M flowchart
/// (Fig 3.18) including all degenerate cases where v_max or a_max cannot be reached.
/// Nonzero initial and final velocities are supported via the general Sec 3.4.1
/// formulation; the zero-velocity case reduces to the Sec 3.4.3 special case.
/// Negative displacement is handled via sigma transformation (eq 3.31-3.33).
///
/// The 7 segments are: jerk(+), const-accel, jerk(-), cruise, jerk(-), const-decel, jerk(+).
///
/// @cite biagiotti2009 -- Sec. 3.4, eq. (3.17)-(3.33), Fig. 3.18, p.79-96

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

/// @brief Double-S (7-segment) velocity profile bounding v, a, and j.
/// @cite biagiotti2009 -- Sec. 3.4, eq. (3.17)-(3.33), Fig. 3.18, p.79-96
template <ctrlpp_floating_scalar Scalar>
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
        , v_max_{cfg.v_max}
        , a_max_{cfg.a_max}
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

        // Boundary velocities in the positive-displacement frame: a move in the
        // -q direction with v0 < 0 maps to +v0 here, so the same profile math
        // serves both directions and evaluate() folds the sign back with sigma_.
        pv0_ = sigma_ * cfg.v0;
        pv1_ = sigma_ * cfg.v1;

        if (cfg.v0 == Scalar{0} && cfg.v1 == Scalar{0}) {
            // Zero boundary velocities: the symmetric special case.
            // @cite biagiotti2009 -- Sec. 3.4.3, p.88-91
            compute_zero_bc(h, v_max, a_max, cfg.j_max);
        } else {
            // General nonzero boundary velocities.
            // @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.27), p.79-85
            compute_general_bc(h, v_max, a_max, cfg.j_max);
        }
    }

    /// @brief Construct a profile, reporting the commands this family cannot realize.
    ///
    /// A seven-segment profile cannot sweep less ground than the fastest
    /// admissible transition from the larger of the two positive-frame boundary
    /// velocities to the smaller one: covering less would require overshooting
    /// the target and coming back, which is a different velocity profile shape.
    /// Rejections:
    ///  * commanded displacement below that minimum ->
    ///    trajectory_error::unreachable_boundary_velocity
    ///
    /// The plain constructor stays available for the callers that have already
    /// established their command is realizable. On a command that is not, it
    /// yields a stationary zero-duration profile rather than a finite profile
    /// that does not traverse its own displacement.
    ///
    /// @cite biagiotti2009 -- Sec. 3.4.1, p.79-85
    [[nodiscard]] static auto try_create(config const& cfg)
        -> ctrlpp::expected<double_s_trajectory, trajectory_error>
    {
        double_s_trajectory profile{cfg};
        if (!profile.realizable_) {
            return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
        }
        return profile;
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
    Scalar v0_{}, v1_{};
    Scalar pv0_{}, pv1_{}; ///< boundary velocities in the positive-displacement frame
    Scalar v_lim_{};
    Scalar a_lim_a_{};
    Scalar a_lim_d_{};
    Scalar v_max_{};
    Scalar a_max_{};
    Scalar j_max_{};
    Scalar T_j1_{}, T_a_{}, T_v_{}, T_d_{}, T_j2_{};
    Scalar T_{};
    bool degenerate_{false};
    bool realizable_{true};

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

    /// @brief Compute phase durations for the general nonzero-boundary-velocity case.
    ///
    /// Works entirely in the positive-displacement frame with the sigma-normalized
    /// boundary velocities pv0_, pv1_. Reduces exactly to the zero-boundary special
    /// case when pv0_ = pv1_ = 0.
    ///
    /// @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.27), Fig. 3.18, p.79-85
    void compute_general_bc(Scalar h, Scalar v_max, Scalar a_max, Scalar j_max)
    {
        auto const pv_hi = std::max(pv0_, pv1_);
        auto const pv_lo = std::min(pv0_, pv1_);

        // Steps 1-3: try the shape that reaches the velocity limit. A cruise
        // segment at v_max only exists when the limit is at least as large as
        // both boundary velocities; when it is not, the peak must lie above both
        // of them and only the cruise-free construction below can express that.
        // @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.21), p.80-81
        if (v_max > pv_hi) {
            auto const accel = make_velocity_ramp(v_max - pv0_, a_max, j_max);
            auto const decel = make_velocity_ramp(v_max - pv1_, a_max, j_max);
            auto const swept = (pv0_ + v_max) * accel.T / Scalar{2}
                             + (v_max + pv1_) * decel.T / Scalar{2};
            auto const T_v = (h - swept) / v_max;

            if (T_v > Scalar{0}) {
                T_j1_ = accel.T_j;
                T_a_ = accel.T;
                a_lim_a_ = accel.a_lim;
                T_j2_ = decel.T_j;
                T_d_ = decel.T;
                a_lim_d_ = decel.a_lim;
                T_v_ = T_v;
                v_lim_ = v_max;
                degenerate_ = (accel.a_lim < a_max) || (decel.a_lim < a_max);
                T_ = T_a_ + T_v_ + T_d_;
                return;
            }
        }

        // Step 4: v_max is not reached, so there is no cruise segment and the
        // peak velocity is the only remaining unknown. Both ramps run toward the
        // peak, so it cannot lie below either boundary velocity.
        //
        // Where the peak is nonnegative the swept distance grows with it, so its
        // smallest value is reached when the ramp attached to the larger boundary
        // velocity vanishes -- the profile is then the fastest admissible
        // transition from the larger boundary velocity to the smaller one. A
        // command below that distance is not realizable by any profile of this
        // shape and is reported as such rather than clamped into a profile that
        // does not traverse its own command. When the larger boundary velocity is
        // itself negative that vanishing-ramp distance is not positive, so a
        // positive command always clears it and nothing is ever rejected there.
        //
        // The unknown carried through this step is the peak's RISE above the
        // larger boundary velocity, never the peak itself. The rise is what every
        // ramp duration depends on, and it is routinely orders of magnitude
        // smaller than the boundary velocity it sits on, so storing the peak and
        // subtracting the boundary velocity back out of it would throw away most
        // of the rise's significant digits before the ramps ever see it.
        auto const h_min = no_cruise_displacement(Scalar{0}, pv_hi, pv_lo, a_max, j_max);
        if (h < h_min) {
            set_unrealizable();
            return;
        }

        // Step 5: both ramps reach a_max once the rise is at least a_max^2 / j_max,
        // and there the swept distance is a quadratic in the rise x:
        //   x^2 + (2 pv_hi + a_max^2 / j_max) x
        //       + (pv_hi^2 - pv_lo^2) / 2
        //       + (a_max^2 / j_max) (3 pv_hi + pv_lo) / 2 - a_max h = 0
        // whose larger root is the admissible one. Derived in place by
        // substituting the trapezoidal ramp durations of make_velocity_ramp into
        // the swept distance and rewriting the result about pv_hi; it is the
        // cruise-free case of the general nonzero-boundary construction, and its
        // discriminant reduces to the same expression that case is usually stated
        // with. The root is taken in the form that avoids subtracting two nearly
        // equal terms.
        // @cite biagiotti2009 -- Sec. 3.4.1, p.83-84
        auto const a_sq_over_j = a_max * a_max / j_max;

        Scalar rise{};
        if (h >= no_cruise_displacement(a_sq_over_j, pv_hi, pv_lo, a_max, j_max)) {
            auto const linear = Scalar{2} * pv_hi + a_sq_over_j;
            auto const constant = (pv_hi - pv_lo) * (pv_hi + pv_lo) / Scalar{2}
                                + a_sq_over_j * (Scalar{3} * pv_hi + pv_lo) / Scalar{2}
                                - a_max * h;
            auto const root = std::sqrt(linear * linear - Scalar{4} * constant);
            rise = (linear >= Scalar{0}) ? (Scalar{-2} * constant / (linear + root))
                                         : ((root - linear) / Scalar{2});
        } else {
            rise = solve_no_cruise_rise(h, pv_hi, pv_lo, a_max, j_max);
        }

        auto const v_peak = pv_hi + rise;
        auto const near_ramp = make_velocity_ramp(rise, a_max, j_max);
        auto const far_ramp = make_velocity_ramp(rise + (pv_hi - pv_lo), a_max, j_max);
        auto const& accel = (pv0_ >= pv1_) ? near_ramp : far_ramp;
        auto const& decel = (pv0_ >= pv1_) ? far_ramp : near_ramp;

        T_j1_ = accel.T_j;
        T_a_ = accel.T;
        a_lim_a_ = accel.a_lim;
        T_j2_ = decel.T_j;
        T_d_ = decel.T;
        a_lim_d_ = decel.a_lim;
        T_v_ = Scalar{0};
        v_lim_ = v_peak;
        degenerate_ = true;
        T_ = T_a_ + T_v_ + T_d_;
    }

    /// @brief Shape of a jerk-limited ramp between two velocities whose endpoint
    /// accelerations are both zero.
    struct velocity_ramp
    {
        Scalar T{};     ///< total ramp duration
        Scalar T_j{};   ///< duration of each constant-jerk sub-segment
        Scalar a_lim{}; ///< acceleration magnitude actually reached
    };

    /// @brief Fastest admissible ramp spanning a velocity change.
    ///
    /// The ramp is jerk, constant acceleration, jerk. When the velocity change is
    /// too small to build the acceleration up to a_max the constant-acceleration
    /// sub-segment vanishes, the ramp is triangular in acceleration, and it peaks
    /// at j_max * T_j instead. Both branches are exact closed forms; there is no
    /// sub-case left over for which the acceleration bound has to be lowered and
    /// the solve retried.
    ///
    /// @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.20), p.80
    static auto make_velocity_ramp(Scalar dv, Scalar a_max, Scalar j_max) -> velocity_ramp
    {
        auto const d = std::abs(dv);
        if (d * j_max < a_max * a_max) {
            auto const T_j = std::sqrt(d / j_max);
            return {Scalar{2} * T_j, T_j, j_max * T_j};
        }
        auto const T_j = a_max / j_max;
        return {d / a_max + T_j, T_j, a_max};
    }

    /// @brief Distance swept by the two ramps of a cruise-free profile whose peak
    /// velocity rises by `rise` above the larger boundary velocity.
    ///
    /// Each ramp's acceleration is point-symmetric about the ramp's own midpoint,
    /// so the velocity it sweeps has mean exactly the average of its two
    /// endpoints and the swept distance is that mean times the ramp duration.
    /// Derived in place from that symmetry, which is the same argument
    /// compute_q_at_Ta() rests on. Both mean velocities are formed from the rise
    /// rather than from the peak so that a rise far below the boundary velocities
    /// keeps every digit it has.
    static auto no_cruise_displacement(Scalar rise, Scalar pv_hi, Scalar pv_lo, Scalar a_max, Scalar j_max) -> Scalar
    {
        auto const far_rise = rise + (pv_hi - pv_lo);
        auto const near_ramp = make_velocity_ramp(rise, a_max, j_max);
        auto const far_ramp = make_velocity_ramp(far_rise, a_max, j_max);
        return (Scalar{2} * pv_hi + rise) * near_ramp.T / Scalar{2}
             + (Scalar{2} * pv_lo + far_rise) * far_ramp.T / Scalar{2};
    }

    /// @brief Rise of the cruise-free profile that sweeps exactly h.
    ///
    /// Reached only below the region where both ramps attain a_max, which the
    /// closed form above covers. There the ramp attached to the larger boundary
    /// velocity is triangular in acceleration, so the two ramps carry different
    /// square roots of the rise, the swept distance is not a polynomial in it,
    /// and no elementary closed form for it exists.
    ///
    /// The bracket spans the rises for which that ramp exists without reaching
    /// a_max, and the swept distance grows across it. Termination is algebraic,
    /// not a trip count: each step either halves the rise bracket until its
    /// midpoint falls on an endpoint, or fails to shrink the swept-distance
    /// bracket, which cannot recur once that bracket has reached the resolution
    /// of the scalar type.
    static auto solve_no_cruise_rise(Scalar h, Scalar pv_hi, Scalar pv_lo, Scalar a_max, Scalar j_max) -> Scalar
    {
        auto const swept_at = [&](Scalar rise) {
            return no_cruise_displacement(rise, pv_hi, pv_lo, a_max, j_max);
        };

        // The bracket starts where the peak velocity reaches zero, not where the
        // near ramp vanishes. Differentiating the swept distance shows it grows
        // with the rise for every nonnegative peak, in both ramp shapes, while
        // below a zero peak the profile moves backwards and the distance can fall
        // as the rise grows. Starting at a zero peak therefore makes the bracket
        // monotone, and it costs no root: a peak below zero leaves both ramp mean
        // velocities negative, so it sweeps no positive distance at all.
        auto lo = std::max(Scalar{0}, -pv_hi);
        auto hi = a_max * a_max / j_max;
        auto swept_lo = swept_at(lo);
        auto swept_hi = swept_at(hi);
        auto width = swept_hi - swept_lo;

        for (;;) {
            auto const mid = lo + (hi - lo) / Scalar{2};
            if (!(mid > lo) || !(mid < hi)) {
                break;
            }
            auto const swept_mid = swept_at(mid);
            if (swept_mid < h) {
                lo = mid;
                swept_lo = swept_mid;
            } else {
                hi = mid;
                swept_hi = swept_mid;
            }
            auto const next_width = swept_hi - swept_lo;
            if (!(next_width < width)) {
                break;
            }
            width = next_width;
        }

        return (std::abs(swept_lo - h) <= std::abs(swept_hi - h)) ? lo : hi;
    }

    /// @brief Record a command this profile family cannot realize.
    ///
    /// The profile becomes stationary, so evaluate() holds the start position for
    /// a zero duration rather than reporting a traversal that never happens.
    /// try_create() turns this state into trajectory_error.
    void set_unrealizable()
    {
        realizable_ = false;
        degenerate_ = true;
        v_lim_ = Scalar{0};
        a_lim_a_ = Scalar{0};
        a_lim_d_ = Scalar{0};
        T_j1_ = Scalar{0};
        T_a_ = Scalar{0};
        T_v_ = Scalar{0};
        T_d_ = Scalar{0};
        T_j2_ = Scalar{0};
        T_ = Scalar{0};
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

    /// @brief Evaluate acceleration-phase segments 1-3.
    /// @cite biagiotti2009 -- Sec. 3.4, eq. (3.30a)-(3.30c), p.85
    void eval_accel_phase(Scalar tc, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        auto const j = j_max_;
        auto const t1 = T_j1_;
        auto const t2 = T_a_ - T_j1_;

        if (tc < t1) {
            // Segment 1: jerk(+), velocity starts at pv0_.
            q = pv0_ * tc + j * tc * tc * tc / Scalar{6};
            dq = pv0_ + j * tc * tc / Scalar{2};
            ddq = j * tc;
        } else if (tc < t2) {
            // Segment 2: constant accel, carrying the pv0_ offset forward.
            auto const t = tc - T_j1_;
            auto const q1 = pv0_ * T_j1_ + j * T_j1_ * T_j1_ * T_j1_ / Scalar{6};
            auto const dq1 = pv0_ + j * T_j1_ * T_j1_ / Scalar{2};
            q = q1 + dq1 * t + a_lim_a_ * t * t / Scalar{2};
            dq = dq1 + a_lim_a_ * t;
            ddq = a_lim_a_;
        } else {
            auto const dt = T_a_ - tc;
            auto const q_at_Ta = compute_q_at_Ta();
            q = q_at_Ta - v_lim_ * dt + j * dt * dt * dt / Scalar{6};
            dq = v_lim_ - j * dt * dt / Scalar{2};
            ddq = j * dt;
        }
    }

    /// @brief Evaluate deceleration-phase segments 5-7.
    /// @cite biagiotti2009 -- Sec. 3.4, eq. (3.30e)-(3.30g), p.86
    void eval_decel_phase(Scalar tc, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        auto const j = j_max_;
        auto const dt_end = T_ - tc;
        auto const h = std::abs(q1_ - q0_);
        auto const t5 = T_a_ + T_v_ + T_j2_;

        if (dt_end > T_d_ - T_j2_) {
            auto const q_at_Tv_end = compute_q_at_Ta() + v_lim_ * T_v_;
            auto const t = tc - (T_a_ + T_v_);
            q = q_at_Tv_end + v_lim_ * t - j * t * t * t / Scalar{6};
            dq = v_lim_ - j * t * t / Scalar{2};
            ddq = -j * t;
        } else if (dt_end > T_j2_) {
            auto const t5_local = T_j2_;
            auto const v_at_t5 = v_lim_ - j * t5_local * t5_local / Scalar{2};
            auto const a_at_t5 = -j * t5_local;
            auto const q_at_t5 = compute_q_at_Ta() + v_lim_ * T_v_
                                 + v_lim_ * t5_local - j * t5_local * t5_local * t5_local / Scalar{6};
            auto const dt_from_t5 = tc - t5;
            q = q_at_t5 + v_at_t5 * dt_from_t5 + a_at_t5 * dt_from_t5 * dt_from_t5 / Scalar{2};
            dq = v_at_t5 + a_at_t5 * dt_from_t5;
            ddq = a_at_t5;
        } else {
            // Segment 7: final jerk(+), velocity ends at pv1_ (not zero).
            q = h - pv1_ * dt_end - j * dt_end * dt_end * dt_end / Scalar{6};
            dq = pv1_ + j * dt_end * dt_end / Scalar{2};
            ddq = -j * dt_end;
        }
    }

    /// @brief Evaluate in the positive-displacement frame (sigma=+1).
    /// @cite biagiotti2009 -- Sec. 3.4, eq. (3.30a)-(3.30g), p.85-86
    void evaluate_positive_frame(Scalar tc, Scalar& q, Scalar& dq, Scalar& ddq) const
    {
        auto const t3 = T_a_;
        auto const t4 = T_a_ + T_v_;

        if (tc < t3) {
            eval_accel_phase(tc, q, dq, ddq);
        } else if (tc < t4) {
            auto const q_at_Ta = compute_q_at_Ta();
            auto const t = tc - T_a_;
            q = q_at_Ta + v_lim_ * t;
            dq = v_lim_;
            ddq = Scalar{0};
        } else {
            eval_decel_phase(tc, q, dq, ddq);
        }
    }

    /// @brief Compute position at end of acceleration phase (t = T_a).
    auto compute_q_at_Ta() const -> Scalar
    {
        // At T_a, velocity = v_lim, acceleration = 0. The acceleration profile of
        // the phase is symmetric about its midpoint, so the velocity ramp is
        // point-symmetric and its mean value is (pv0_ + v_lim) / 2. The swept
        // displacement is therefore the trapezoidal area (pv0_ + v_lim)/2 * T_a,
        // which is exact for nonzero initial velocity and reduces to v_lim*T_a/2
        // when pv0_ = 0.
        return (pv0_ + v_lim_) * T_a_ / Scalar{2};
    }
};

static_assert(trajectory_segment<double_s_trajectory<double>, double, 1>);
static_assert(trajectory_segment<double_s_trajectory<float>, float, 1>);

}

#endif
