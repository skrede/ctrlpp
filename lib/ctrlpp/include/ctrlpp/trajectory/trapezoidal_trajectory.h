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

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/trajectory_segment.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cmath>
#include <limits>
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
        , v_max_{cfg.v_max}
        , a_{cfg.a_max}
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
            if (!std::isfinite(a)) {
                // The raise is a division by the commanded displacement, so it
                // leaves the representable range exactly when that displacement
                // is too small to reconcile the two boundary velocities at any
                // acceleration the scalar type can hold -- at a zero
                // displacement, at any acceleration whatsoever. Every quantity
                // downstream is formed from this value multiplied by a phase
                // duration derived from it, so an unbounded one does not
                // propagate as an unbounded duration: it propagates as NaN, in
                // the cruise duration first and in the total from there.
                set_unrealizable();
                return;
            }
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
        a_ = a;
    }

    /// @brief Construct a profile, reporting the commands this family cannot realize.
    ///
    /// Both ramps of a three-phase profile run toward one cruise velocity that
    /// lies at or above each boundary velocity, so the ground the profile sweeps
    /// is at least the ground the transition between those two velocities
    /// already sweeps. A command below that is not feasible at the commanded
    /// acceleration, and B&M's remedy is to raise the acceleration to the
    /// smallest value which makes the two feasible together, eq. (3.15). At a
    /// zero commanded displacement no acceleration is large enough, and just
    /// above zero the value the remedy asks for is no longer representable.
    /// Rejections:
    ///  * a commanded displacement the two boundary velocities cannot be
    ///    reconciled with at a representable acceleration ->
    ///    trajectory_error::unreachable_boundary_velocity
    ///
    /// The plain constructor stays available for the callers that have already
    /// established their command is realizable. On a command that is not, it
    /// yields a stationary zero-duration profile rather than the NaN duration and
    /// NaN evaluation the unbounded acceleration would otherwise propagate.
    ///
    /// @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.14)-(3.15), p.72
    [[nodiscard]] static auto try_create(config const& cfg)
        -> ctrlpp::expected<trapezoidal_trajectory, trajectory_error>
    {
        trapezoidal_trajectory profile{cfg};
        if (!profile.realizable_) {
            return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
        }
        return profile;
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

    /// @brief Rescale the profile to a longer duration for multi-axis synchronization.
    ///
    /// The profile is rebuilt at a lower cruise velocity, never patched. The
    /// commanded displacement, both boundary velocities, and the acceleration
    /// magnitude are held fixed and the cruise velocity that realizes T_new is
    /// solved in closed form, so the traversed displacement and the terminal
    /// velocity hold by construction. The stored duration stays the sum of the
    /// three realized phase durations and is never assigned the requested value;
    /// it lands within a few units in the last place of it.
    ///
    /// Rejections, checked in order:
    ///  * a duration below the current one -> trajectory_error::duration_shorter_than_current
    ///  * NaN, infinite, or non-positive T_new -> trajectory_error::non_positive_duration
    ///  * a duration the displacement, limits, and boundary velocities cannot
    ///    realize together -> trajectory_error::unreachable_duration
    ///
    /// A request equal to the current duration succeeds and changes nothing, which
    /// is the path the slowest axis of a synchronized set always takes.
    ///
    /// @cite biagiotti2009 -- Sec. 5.3, eq. (5.13)-(5.14) -- time scaling for synchronization
    [[nodiscard]] auto rescale_to(Scalar T_new) -> ctrlpp::expected<void, trajectory_error>
    {
        auto const solved = solve_rescale(T_new);
        if (!solved.has_value()) {
            return ctrlpp::unexpected(solved.error());
        }

        auto const& s = solved.value();
        v_v_ = s.v_v;
        a_a_ = s.a_a;
        a_d_ = s.a_d;
        T_a_ = s.T_a;
        T_v_ = s.T_v;
        T_d_ = s.T_d;
        T_ = s.T;
        triangular_ = s.triangular;
        return {};
    }

    /// @brief Report whether rescale_to(T_new) would succeed, without mutating.
    ///
    /// Runs the identical solve rescale_to() runs and discards the result, so the
    /// two cannot disagree: identical inputs traverse identical code with no
    /// intervening state. That is what lets a multi-axis synchronization check
    /// every axis before it commits any of them.
    [[nodiscard]] auto can_rescale_to(Scalar T_new) const -> ctrlpp::expected<void, trajectory_error>
    {
        auto const solved = solve_rescale(T_new);
        if (!solved.has_value()) {
            return ctrlpp::unexpected(solved.error());
        }
        return {};
    }

  private:
    /// @brief Record a command this profile family cannot realize.
    ///
    /// The profile becomes stationary, so evaluate() holds one position for a
    /// zero duration rather than reporting a traversal that never happens. The
    /// two ramp rates and the two boundary velocities are cleared with it: they
    /// are the shape of a traversal there is none of, and leaving them set would
    /// have evaluate() report a moving axis at a standstill.
    /// try_create() turns this state into trajectory_error.
    void set_unrealizable()
    {
        realizable_ = false;
        triangular_ = true;
        v0_ = Scalar{0};
        v1_ = Scalar{0};
        v_v_ = Scalar{0};
        a_a_ = Scalar{0};
        a_d_ = Scalar{0};
        T_a_ = Scalar{0};
        T_v_ = Scalar{0};
        T_d_ = Scalar{0};
        T_ = Scalar{0};
    }

    /// @brief Complete profile state produced by a rescaling solve.
    struct rescaled_state
    {
        Scalar v_v{};
        Scalar a_a{};
        Scalar a_d{};
        Scalar T_a{};
        Scalar T_v{};
        Scalar T_d{};
        Scalar T{};
        bool triangular{};
    };

    /// @brief Distance the two ramps sweep together in the ramp-through shape.
    ///
    /// Both ramps run in the same direction there, so one covers the exact
    /// trapezoid area between the larger boundary velocity and the cruise velocity
    /// and the other covers the area between the cruise velocity and the smaller
    /// one; the cruise velocity cancels and what remains is independent of it.
    /// The expression is SIGNED and is negative when the smaller boundary velocity
    /// is the larger of the two in magnitude, which happens whenever the axis
    /// starts or ends moving away from its target faster than it moves toward it.
    /// The |v0^2 - v1^2| spelling loses exactly that case.
    [[nodiscard]] auto ramp_through_distance() const -> Scalar
    {
        auto const v_lo = std::min(v0_, v1_);
        auto const v_hi = std::max(v0_, v1_);
        return (v_hi - v_lo) * (v_hi + v_lo) / (Scalar{2} * a_);
    }

    /// @brief Solve the rescaled profile, or report why the request is not realizable.
    ///
    /// Everything below works in the positive-displacement frame the constructor
    /// establishes, with the transformed boundary velocities v0_ and v1_, the
    /// commanded displacement h = |q1 - q0|, and the acceleration magnitude a_.
    /// A ramp's swept distance is taken as the exact trapezoid area, the mean of
    /// its two velocities times its duration; the |v^2 - v0^2| / (2a) shortcut is
    /// equal to that only when the two velocities sum to a positive value, and a
    /// transformed boundary velocity may be negative when the axis starts out
    /// moving away from its target.
    ///
    /// The total duration
    ///     T(v) = |v - v0| / a + (h - d_a - d_d) / v + |v - v1| / a
    /// is continuous and strictly decreasing in the cruise velocity v over
    /// (0, v_tri], with v_tri = sqrt(a h + (v0^2 + v1^2) / 2) the triangular
    /// (time-optimal) cruise velocity, so at most one of the three shapes below can
    /// contain the root and the shape is identified by evaluating T at the two
    /// shape boundaries rather than by trying each in turn. Each expression is
    /// derived here by substituting the ramp durations and trapezoid areas into
    /// T(v) and clearing the division by v; only the parametrization is taken from
    /// the reference.
    ///
    ///  * plateau (v >= max(v0, v1)): both ramps run toward the cruise velocity, so
    ///    T(v) = (v - v0 - v1) / a + D_A / v with D_A = h + (v0^2 + v1^2) / (2a),
    ///    giving v^2 - [(v0 + v1) + a T] v + [a h + (v0^2 + v1^2) / 2] = 0. The
    ///    smaller root is the admissible one; the larger lies above v_tri.
    ///  * ramp-through (min(v0, v1) < v < max(v0, v1)): both ramps run in the same
    ///    direction, so their total duration |v1 - v0| / a and their total swept
    ///    distance (max(v0,v1)^2 - min(v0,v1)^2) / (2a), which is SIGNED, are both
    ///    independent of v. What remains is linear in 1 / v with exactly one root,
    ///    no discriminant and no selection.
    ///  * valley (v <= min(v0, v1)): both ramps run away from the cruise velocity,
    ///    so T(v) = (v0 + v1 - v) / a + D_C / v with D_C = h - (v0^2 + v1^2) / (2a),
    ///    giving v^2 + [a T - (v0 + v1)] v + [(v0^2 + v1^2) / 2 - a h] = 0. The
    ///    larger root is the admissible one. This shape is inside the
    ///    parametrization and is emitted, not rejected.
    ///
    /// Reachability is decided ahead of the solve and is free of chosen constants.
    /// The cruise duration (h - d_a - d_d) / v is what runs out: the shape whose
    /// validity interval reaches down toward a vanishing cruise velocity fixes the
    /// supremum of the reachable durations, and whether that supremum is finite is
    /// the sign of that shape's own residual displacement term.
    ///
    /// @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.13a)-(3.13c), p.71 -- the
    ///   three-phase parametrization these expressions solve for the cruise velocity
    [[nodiscard]] auto solve_rescale(Scalar T_new) const
        -> ctrlpp::expected<rescaled_state, trajectory_error>
    {
        auto const current = rescaled_state{
            v_v_, a_a_, a_d_, T_a_, T_v_, T_d_, T_, triangular_};

        // Exact equality first: a synchronized set passes the slowest axis a
        // bit-exact copy of its own duration, so that axis lands here rather than
        // in the shortening rejection below. There is no float-equality fragility
        // in that -- the value compared is the same object's own stored duration.
        if (T_new == T_) {
            return current;
        }
        if (!std::isfinite(T_new) || T_new <= Scalar{0}) {
            return ctrlpp::unexpected(trajectory_error::non_positive_duration);
        }
        if (T_new < T_) {
            return ctrlpp::unexpected(trajectory_error::duration_shorter_than_current);
        }
        if (!realizable_) {
            // A command this family cannot realize has no traversal to slow
            // down, and the stationary stand-in left in its place traverses
            // nothing at any duration.
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        auto const h = std::abs(q1_ - q0_);
        auto const a = a_;
        auto const v_lo = std::min(v0_, v1_);
        auto const v_hi = std::max(v0_, v1_);
        auto const v_sum_sq = v0_ * v0_ + v1_ * v1_;
        auto const v_tri = std::sqrt(a * h + v_sum_sq / Scalar{2});
        auto const v_top = std::min(v_tri, v_max_);

        // A cruise velocity below the larger boundary velocity cannot be reached by
        // ramping toward it, so no shape of this family serves the configuration.
        if (v_top < v_hi) {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        // The ramp-through shape's residual displacement is a difference of two
        // nearly equal quantities whenever the constructor had to raise the
        // acceleration to make the boundary velocities feasible, because that
        // raised value is defined by making the ramp distance equal the commanded
        // displacement. Five chained roundings form the ramp distance and a sixth
        // forms the subtraction below. Each is worth one unit in the last place at
        // the scale of the largest operand that entered, which is NOT the ramp
        // distance itself: that distance is a difference of the two squared
        // boundary velocities and inherits their scale, so the floor is written
        // against (v0^2 + v1^2) / (2a) and the commanded displacement. Below it the
        // residual is zero to the precision available, and a cruise velocity
        // divided out of it would be noise.
        constexpr int residual_rounding_ops = 6;
        auto const residual_scale = std::max(h, v_sum_sq / (Scalar{2} * a));
        auto const residual_floor = Scalar{residual_rounding_ops}
                                    * std::numeric_limits<Scalar>::epsilon() * residual_scale;
        auto const ramp_residual = h - ramp_through_distance();
        auto const ramp_resolved = ramp_residual > residual_floor;

        // Reachability. The shape whose validity interval reaches down toward a
        // vanishing cruise velocity fixes the supremum of the reachable durations.
        // The duration diverges as that velocity vanishes exactly when that shape's
        // residual displacement is positive; otherwise the cruise duration reaches
        // zero first and pins a finite supremum.
        bool unbounded{};
        Scalar T_sup{};
        if (v_hi <= Scalar{0}) {
            // Both boundary velocities point away from the target, so the plateau
            // shape covers the whole positive range and its residual displacement
            // h + (v0^2 + v1^2) / (2a) is a sum of nonnegative terms.
            unbounded = (h > Scalar{0});
            T_sup = -(v0_ + v1_) / a;
        } else if (v_lo <= Scalar{0}) {
            // The ramp-through shape covers the range down to a vanishing cruise
            // velocity, where the distance its two ramps sweep is all there is.
            unbounded = ramp_resolved;
            T_sup = std::abs(v1_ - v0_) / a;
        } else {
            // The valley shape reaches the vanishing-cruise limit at
            // v_min = sqrt((v0^2 + v1^2) / 2 - a h), where T = (v0 + v1 - 2 v_min) / a.
            auto const residual = h - v_sum_sq / (Scalar{2} * a);
            unbounded = (residual > Scalar{0});
            auto const v_min = std::sqrt(std::max(Scalar{0}, -a * residual));
            T_sup = (v0_ + v1_ - Scalar{2} * v_min) / a;
        }
        if (!unbounded && T_new > T_sup) {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        auto const duration_at = [&](Scalar v) -> Scalar {
            auto const T_a = std::abs(v - v0_) / a;
            auto const T_d = std::abs(v - v1_) / a;
            auto const d_a = (v0_ + v) * T_a / Scalar{2};
            auto const d_d = (v1_ + v) * T_d / Scalar{2};
            return T_a + (h - d_a - d_d) / v + T_d;
        };

        // Shape selection by monotonicity: T decreases as the cruise velocity
        // grows, so the shape whose duration interval brackets T_new is the one
        // that contains the root. A boundary at or below zero means that shape's
        // interval already reaches the vanishing-cruise limit.
        enum class shape
        {
            plateau,
            ramp_through,
            valley
        };
        auto const selected = (v_hi <= Scalar{0} || T_new <= duration_at(v_hi))
                                  ? shape::plateau
                              : (v_lo <= Scalar{0} || T_new <= duration_at(v_lo))
                                  ? shape::ramp_through
                                  : shape::valley;

        Scalar v{};
        if (selected == shape::plateau) {
            auto const b = (v0_ + v1_) + a * T_new;
            auto const c = a * h + v_sum_sq / Scalar{2};
            auto const disc = b * b - Scalar{4} * c;
            if (!(disc >= Scalar{0})) {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            auto const root = std::sqrt(disc);
            // Smaller root, in the form that avoids subtracting two nearly equal terms.
            v = (b > Scalar{0}) ? (Scalar{2} * c / (b + root)) : ((b - root) / Scalar{2});
        } else if (selected == shape::ramp_through) {
            auto const denom = T_new - std::abs(v1_ - v0_) / a;
            if (!ramp_resolved || !(denom > Scalar{0})) {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            v = ramp_residual / denom;
        } else {
            auto const b = a * T_new - (v0_ + v1_);
            auto const c = v_sum_sq / Scalar{2} - a * h;
            auto const disc = b * b - Scalar{4} * c;
            if (!(disc >= Scalar{0})) {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            auto const root = std::sqrt(disc);
            // Larger root, in the form that avoids subtracting two nearly equal terms.
            v = (b > Scalar{0}) ? (Scalar{-2} * c / (b + root)) : ((root - b) / Scalar{2});
        }

        // A root is accepted only inside its own shape's validity interval, with
        // every phase duration nonnegative and the cruise velocity within the
        // velocity limit. A root failing any of these is a rejection, never a
        // clamped value.
        auto const in_shape = (selected == shape::plateau)        ? (v >= v_hi)
                              : (selected == shape::ramp_through) ? (v >= v_lo && v <= v_hi)
                                                                  : (v <= v_lo);
        if (!std::isfinite(v) || !(v > Scalar{0}) || v > v_max_ || !in_shape) {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        auto const T_a = std::abs(v - v0_) / a;
        auto const T_d = std::abs(v - v1_) / a;
        auto const d_a = (v0_ + v) * T_a / Scalar{2};
        auto const d_d = (v1_ + v) * T_d / Scalar{2};
        auto const T_v = (h - d_a - d_d) / v;
        if (!(T_a >= Scalar{0}) || !(T_v >= Scalar{0}) || !(T_d >= Scalar{0})) {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        // The ramps carry the sign of the velocity change they realize: a valley
        // shape decelerates away from v0 and accelerates back up to v1, which the
        // evaluator renders from these signed rates.
        return rescaled_state{
            .v_v = v,
            .a_a = (v >= v0_) ? a : -a,
            .a_d = (v >= v1_) ? a : -a,
            .T_a = T_a,
            .T_v = T_v,
            .T_d = T_d,
            .T = T_a + T_v + T_d,
            .triangular = (T_v <= Scalar{0}),
        };
    }

    Scalar q0_{};
    Scalar q1_{};
    Scalar sigma_{};
    Scalar v0_{};
    Scalar v1_{};
    Scalar v_v_{};
    Scalar v_max_{};
    Scalar a_{};
    Scalar a_a_{};
    Scalar a_d_{};
    Scalar T_a_{};
    Scalar T_v_{};
    Scalar T_d_{};
    Scalar T_{};
    bool triangular_{};
    bool realizable_{true};
};

static_assert(trajectory_segment<trapezoidal_trajectory<double>, double, 1>);

}

#endif
