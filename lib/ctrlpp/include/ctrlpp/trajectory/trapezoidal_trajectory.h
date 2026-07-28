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

/// @brief Disposition of a constructed trapezoidal profile: the acceleration
/// limit the caller commanded, against the one the profile realizes.
///
/// This is a report on a SUCCESSFUL construction, not a failure. Where the two
/// boundary velocities cannot be reconciled over the commanded displacement at
/// the commanded limit, the construction raises the limit to the smallest value
/// that makes them feasible together (B&M eq. (3.15)) and builds a profile that
/// is correct and limit-respecting under the RAISED value. A caller whose
/// acceleration limit is physical rather than advisory reads these two fields
/// and decides what to do; a caller for whom it was advisory ignores them.
///
/// The two accelerations are carried rather than a flag, because a supervisory
/// layer that learns only that something was substituted cannot decide
/// anything with it. `realized_acceleration` equals `commanded_acceleration`
/// exactly when nothing was raised.
///
/// @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.14)-(3.15), p.72
template <typename Scalar>
struct [[nodiscard]] trapezoidal_disposition
{
    Scalar commanded_acceleration{};  ///< the acceleration limit passed to `create`
    Scalar realized_acceleration{};   ///< the magnitude every ramp of the profile actually uses
};

/// @brief Trapezoidal (LSPB) velocity profile with introspection.
///
/// Construction goes through `create`, which returns
/// `ctrlpp::expected<trapezoidal_trajectory, trajectory_error>` and is the only
/// way to obtain one. There is no non-validating public constructor: a profile
/// object either satisfies its contract -- a finite, nonnegative duration in
/// every phase and a peak within the commanded velocity limit -- or it was
/// never built. Every rejection therefore reaches the caller as a value it has
/// to inspect, and `evaluate` can rely on the invariants that keep its own
/// arithmetic defined.
///
/// The degenerate triangular case is not a rejection: it is one of the shapes
/// this family covers, and it is reported by `is_triangular()`.
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

    /// @brief Construct a profile, or report why the command has none.
    ///
    /// The only construction path. Rejections, checked in order:
    ///  * a NaN or infinite position or boundary velocity ->
    ///    trajectory_error::non_finite_input
    ///  * a NaN, infinite, or non-positive velocity limit ->
    ///    trajectory_error::non_positive_velocity_limit
    ///  * a NaN, infinite, or non-positive acceleration limit ->
    ///    trajectory_error::non_positive_acceleration_limit
    ///  * a boundary velocity above the velocity limit in magnitude ->
    ///    trajectory_error::boundary_velocity_exceeds_limit
    ///  * a commanded displacement the two boundary velocities cannot be
    ///    reconciled with at a representable acceleration ->
    ///    trajectory_error::unreachable_boundary_velocity
    ///  * a negative phase duration or total duration ->
    ///    trajectory_error::unreachable_boundary_velocity
    ///  * a duration outside the representable range, or a total duration that
    ///    underflowed to zero on a nonzero commanded displacement ->
    ///    trajectory_error::unrepresentable_duration
    ///
    /// The velocity limit is a PRECONDITION on the boundary velocities, not a
    /// bound the construction raises to accommodate them. Raising it would
    /// return a profile that violates a limit the caller stated, which is the
    /// same silent-success defect a rejection exists to prevent, and honoring
    /// the limit while accepting the command would need a ramp that runs
    /// backwards in time: the acceleration phase spans (v_v - v0) / a, and a
    /// cruise velocity held under the limit while v0 sits above it makes that
    /// span negative.
    ///
    /// One outcome deliberately does NOT appear in that list. Where the two
    /// boundary velocities cannot be reconciled over the commanded displacement
    /// at the commanded acceleration limit, the construction raises the limit to
    /// the smallest value that makes them feasible together (eq. (3.15)) and
    /// succeeds. That is a SUCCESS whose realized limit differs from the
    /// commanded one, not a failure: the profile that comes back is correct and
    /// respects the raised limit in every phase, and it is the profile the
    /// command asks for at the only acceleration that can deliver it. It is
    /// reported on the disposition channel -- `disposition()` carries the
    /// commanded and the realized acceleration -- so a caller whose limit is
    /// physical rather than advisory can compare them and decide, and no caller
    /// is made to handle a non-failure through the failure path.
    ///
    /// The last two checks are one guard against two spellings of the same
    /// failure. Both ramps of a three-phase profile run toward one cruise
    /// velocity lying at or above each boundary velocity, so the profile sweeps
    /// at least the ground the transition between those two velocities already
    /// sweeps. B&M's remedy for a shorter command is to raise the acceleration
    /// to the smallest value that makes the two feasible together, eq. (3.15);
    /// at a zero commanded displacement no acceleration is large enough, and
    /// just above zero the value the remedy asks for is no longer
    /// representable. Below the square root of the smallest normal value the
    /// remedy cannot even see the case: a boundary velocity squares to zero
    /// there, so eq. (3.14)'s feasibility test reads as satisfied, the
    /// triangular peak underflows below the boundary velocity it is
    /// analytically bounded by, and the ramp duration formed from the
    /// difference comes out negative. The realized durations are therefore
    /// checked directly rather than inferred from the test that was supposed to
    /// guarantee them.
    ///
    /// @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.14)-(3.15), p.72
    static auto create(config const& cfg)
        -> ctrlpp::expected<trapezoidal_trajectory, trajectory_error>
    {
        if (!std::isfinite(cfg.q0) || !std::isfinite(cfg.q1) || !std::isfinite(cfg.v0)
            || !std::isfinite(cfg.v1)) {
            return ctrlpp::unexpected(trajectory_error::non_finite_input);
        }
        if (!std::isfinite(cfg.v_max) || cfg.v_max <= Scalar{0}) {
            return ctrlpp::unexpected(trajectory_error::non_positive_velocity_limit);
        }
        if (!std::isfinite(cfg.a_max) || cfg.a_max <= Scalar{0}) {
            return ctrlpp::unexpected(trajectory_error::non_positive_acceleration_limit);
        }
        if (std::abs(cfg.v0) > cfg.v_max || std::abs(cfg.v1) > cfg.v_max) {
            return ctrlpp::unexpected(trajectory_error::boundary_velocity_exceeds_limit);
        }
        if (!std::isfinite(solve_acceleration(cfg))) {
            return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
        }

        trapezoidal_trajectory profile{unchecked_t{}, cfg};

        for (auto const phase : profile.phase_durations()) {
            if (phase < Scalar{0}) {
                return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
            }
            if (!std::isfinite(phase)) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
        }
        if (profile.T_ < Scalar{0}) {
            return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
        }
        if (!std::isfinite(profile.T_)) {
            return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
        }
        // A nonzero commanded displacement is not traversed in zero time at any
        // finite velocity, so a total duration that underflowed to zero under
        // one describes no motion the command asked for. The predicate is the
        // physical statement itself and carries no threshold.
        if (cfg.q1 != cfg.q0 && !(profile.T_ > Scalar{0})) {
            return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
        }
        return profile;
    }

    /// @brief Evaluate trajectory at time t, clamped to [0, T].
    ///
    /// Three-phase branching: acceleration, cruise, deceleration.
    /// Deceleration uses backward time (T - t) for numerical precision.
    ///
    /// The clamp below is well defined because no object of this type exists
    /// with a negative duration: `create` is the only construction path and it
    /// rejects one. `std::clamp` with a lower bound above its upper bound is
    /// undefined behavior, not a clamp that returns something unhelpful, so
    /// that construction-time guarantee is what keeps this line safe rather
    /// than a test performed here on every evaluation.
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

    /// @brief The commanded acceleration limit against the one this profile
    /// realizes.
    ///
    /// Fixed when the profile was built and never recomputed: the raise is
    /// decided once, by `solve_acceleration`, and rescaling holds the
    /// acceleration magnitude fixed by construction.
    auto disposition() const -> trapezoidal_disposition<Scalar> const& { return disposition_; }

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
    ///  * a request the expressions that would locate its cruise velocity cannot
    ///    resolve or represent -> trajectory_error::unrepresentable_duration
    ///  * a duration the displacement, limits, and boundary velocities cannot
    ///    realize together -> trajectory_error::unreachable_duration
    ///
    /// The last two are kept apart on purpose. The first says a shape may well
    /// exist and the arithmetic cannot locate it, so the caller should change the
    /// request; the second says no admissible shape exists, so the caller should
    /// change the limits or the command.
    ///
    /// A request equal to the current duration succeeds and changes nothing, which
    /// is the path the slowest axis of a synchronized set always takes.
    ///
    /// @cite biagiotti2009 -- Sec. 5.3, eq. (5.13)-(5.14) -- time scaling for synchronization
    auto rescale_to(Scalar T_new) -> ctrlpp::expected<void, trajectory_error>
    {
        auto const solved = solve_rescale(T_new);
        if (!solved.has_value()) {
            return ctrlpp::unexpected(solved.error());
        }

        auto const& s = *solved;
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
    auto can_rescale_to(Scalar T_new) const -> ctrlpp::expected<void, trajectory_error>
    {
        auto const solved = solve_rescale(T_new);
        if (!solved.has_value()) {
            return ctrlpp::unexpected(solved.error());
        }
        return {};
    }

  private:
    /// @brief Tag selecting the non-validating constructor reserved for `create`.
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    /// @brief Acceleration magnitude the command is realized at.
    ///
    /// The commanded value, except where the two boundary velocities are not
    /// feasible over the commanded displacement at it: B&M eq. (3.14) is the
    /// test and eq. (3.15) is the remedy, which raises the acceleration to the
    /// smallest value that makes the two ramps cover the displacement exactly.
    /// The added unit in the last place keeps the raised value on the feasible
    /// side of the test it was derived from after rounding. The remedy divides
    /// by the commanded displacement, so it leaves the representable range
    /// exactly when that displacement is too small to reconcile the two
    /// boundary velocities at any acceleration the scalar type can hold, and at
    /// a zero displacement at any acceleration whatsoever. `create` tests that
    /// before it builds anything.
    ///
    /// The sign frame does not enter: the test and the remedy are both built
    /// from the squares of the boundary velocities, which the sigma transform
    /// leaves alone.
    ///
    /// @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.14)-(3.15), p.72
    static auto solve_acceleration(config const& cfg) -> Scalar
    {
        auto const abs_h = std::abs(cfg.q1 - cfg.q0);
        auto const v_diff_sq = std::abs(cfg.v0 * cfg.v0 - cfg.v1 * cfg.v1) / Scalar{2};
        if (cfg.a_max * abs_h >= v_diff_sq) {
            return cfg.a_max;
        }
        return v_diff_sq / abs_h + std::numeric_limits<Scalar>::epsilon();
    }

    /// @brief Solve the profile from a configuration `create` has checked.
    ///
    /// Non-validating by construction: it computes the three phase durations and
    /// nothing else, and it is private so the only way to reach it is through
    /// `create`, which decides both before and after whether what came out is a
    /// profile.
    ///
    /// @cite biagiotti2009 -- Sec. 3.2.7, eq. (3.13a)-(3.13c), p.71
    explicit trapezoidal_trajectory(unchecked_t, config const& cfg)
        : q0_{cfg.q0}
        , q1_{cfg.q1}
        , v_max_{cfg.v_max}
        , a_{solve_acceleration(cfg)}
        , disposition_{cfg.a_max, a_}
    {
        auto const h = cfg.q1 - cfg.q0;
        sigma_ = (h >= Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const abs_h = std::abs(h);

        // Transform velocities into the positive-displacement frame
        auto const sv0 = sigma_ * cfg.v0;
        auto const sv1 = sigma_ * cfg.v1;
        auto const a = a_;
        auto const v = cfg.v_max;

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
            //
            // The residual is nonnegative wherever this branch runs, and the
            // guard below is a rounding guard rather than a stand-in for a
            // rejected command: the branch is entered only when the triangular
            // peak reaches the velocity limit, which is a h >= v^2 - (v0^2 +
            // v1^2) / 2, and the two ramp distances substituted from the phase
            // durations sum to exactly (v^2 - (v0^2 + v1^2) / 2) / a. The
            // residual is that inequality's slack, zero at the shape boundary
            // and positive inside.
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
    auto ramp_through_distance() const -> Scalar
    {
        auto const v_lo = std::min(v0_, v1_);
        auto const v_hi = std::max(v0_, v1_);
        return (v_hi - v_lo) * (v_hi + v_lo) / (Scalar{2} * a_);
    }

    /// @brief A total duration paired with the scale its own rounding is measured
    /// against.
    struct duration_with_scale
    {
        Scalar T{};
        Scalar scale{};
    };

    /// @brief Total duration at a cruise velocity, with the largest operand that
    /// entered it, expressed in the units of the duration.
    ///
    /// The scale is NOT the duration alone. The cruise term's numerator is the
    /// commanded displacement less the two swept ramp distances, and that
    /// difference cancels to nothing whenever the ramps sweep almost all of the
    /// displacement -- the same configuration that makes the ramp-through
    /// residual below vanish. A resolution floor written against the duration
    /// alone would be blind to exactly that case, so the largest operand of the
    /// cancelling difference is carried out alongside the result and divided by
    /// the cruise velocity to put it in units of time.
    auto duration_and_scale_at(Scalar v) const -> duration_with_scale
    {
        auto const h = std::abs(q1_ - q0_);
        auto const T_a = std::abs(v - v0_) / a_;
        auto const T_d = std::abs(v - v1_) / a_;
        auto const d_a = (v0_ + v) * T_a / Scalar{2};
        auto const d_d = (v1_ + v) * T_d / Scalar{2};
        auto const T = T_a + (h - d_a - d_d) / v + T_d;
        auto const cruise_scale = std::max({h, std::abs(d_a), std::abs(d_d)}) / std::abs(v);
        return {T, std::max({T_a, T_d, cruise_scale, std::abs(T)})};
    }

    /// @brief Total duration at a cruise velocity.
    ///
    /// Forwards to the sibling above rather than restating the expression, so the
    /// duration and the scale its floor is written against cannot drift apart.
    auto duration_at(Scalar v) const -> Scalar { return duration_and_scale_at(v).T; }

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
    auto solve_rescale(Scalar T_new) const
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
        //
        // Only the two shapes whose boundary velocities straddle zero are decided
        // here. When both are positive the valley is the shape that reaches the
        // vanishing-cruise limit, and its supremum is decided by the sign of its
        // own linear coefficient at the point of solving, below. That test is
        // exact where a precomputed supremum is not: the velocity at the limit is
        // sqrt((v0^2 + v1^2) / 2 - a h), a square root of a difference of two
        // nearly equal quantities, and it therefore carries no significant digits
        // in the very configuration the valley solve has to get right.
        bool unbounded{true};
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
        }
        if (!unbounded && T_new > T_sup) {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        // The two shape boundaries, each with the scale its own rounding is
        // measured against. A boundary at or below zero is not an admissible
        // cruise velocity and its duration is never formed.
        auto const at_hi
            = (v_hi > Scalar{0}) ? duration_and_scale_at(v_hi) : duration_with_scale{};
        auto const at_lo
            = (v_lo > Scalar{0}) ? duration_and_scale_at(v_lo) : duration_with_scale{};

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
        auto const selected = (v_hi <= Scalar{0} || T_new <= at_hi.T)   ? shape::plateau
                              : (v_lo <= Scalar{0} || T_new <= at_lo.T) ? shape::ramp_through
                                                                        : shape::valley;

        // Both solved shapes measure the request against the duration one of the
        // boundary velocities already realizes, and that boundary duration is
        // itself accurate only to a handful of units in the last place. An
        // increment or decrement smaller than that error is indistinguishable
        // from zero, and the shape it selects would be located by noise.
        //
        // duration_and_scale_at chains fifteen arithmetic operations to produce
        // the boundary duration: a subtraction and a division for each of the two
        // ramp durations (four), an addition, a multiplication and a division for
        // each of the two swept ramp distances (six more, ten), two subtractions
        // and a division for the cruise term (three more, thirteen), and two
        // additions for the final sum (fifteen). Every operation is counted
        // whether or not it rounds -- the divisions by two are exact in a binary
        // radix -- so the count bounds the accumulated error from above rather
        // than describing it tightly, which is what a floor requires. Each is
        // worth one unit in the last place at the scale the sibling reports,
        // which is the largest operand that entered rather than the duration
        // alone.
        constexpr int boundary_duration_rounding_ops = 15;

        // Both solved shapes carry the DISTANCE FROM THEIR OWN BOUNDARY as the
        // unknown, never the cruise velocity itself. That distance is what the
        // ramp durations and the cruise term depend on, and it is routinely orders
        // of magnitude smaller than the boundary velocity it sits on, so forming
        // the cruise velocity first and subtracting the boundary velocity back out
        // of it would discard most of the distance's significant digits before the
        // phase durations ever see it. It is the parametrization the double-S
        // profile's general boundary-condition solve already uses for the peak's
        // rise above the larger boundary velocity, for the identical reason.
        //
        // It is also what makes these solves well conditioned rather than merely
        // better spelled, and the reason is one identity. Writing r for the
        // ramp-through residual computed above, the constant terms of the two
        // cruise-velocity forms satisfy, exactly,
        //
        //     [a h + (v0^2 + v1^2) / 2] - v_hi^2 =  a r        (plateau)
        //     [(v0^2 + v1^2) / 2 - a h] - v_lo^2 = -a r        (valley)
        //
        // so each form's discriminant at its own shape boundary is (a r / v)^2
        // while the two terms it subtracts are each about four times v^2. The
        // residual therefore governs both shapes exactly as it governs the
        // ramp-through shape between them, and the floor already derived for it
        // applies to all three unchanged. Neither shape is a milder case than the
        // other: as the residual closes, the plateau's whole velocity range
        // sqrt(v_hi^2 + a r) - v_hi and the valley's whole duration window close
        // with it, so the loss covers the entire branch rather than a
        // neighbourhood of one endpoint.
        //
        // Each rearrangement is derived in place by substituting the shifted
        // cruise velocity into the corresponding duration named in this
        // function's contract above and clearing the division by that velocity;
        // only the underlying three-phase parametrization is taken from the
        // reference cited there.
        //
        // The gates are the same on both, and none of them introduces a constant.
        // A non-positive linear coefficient means the quadratic's two roots are
        // both non-positive, so the shape admits no cruise velocity on the far
        // side of its boundary and no shape answers this request; a negative
        // discriminant means the request lies past the shape's reachable extreme.
        // Both say no admissible shape was found, which is what
        // unreachable_duration reports. The two resolution floors say something
        // different -- a shape may well exist, but the request is below the
        // arithmetic resolution of the expressions that would locate it, so the
        // caller should change the request rather than the limits.
        //
        // The floors are tested FIRST, and the order is load-bearing rather than
        // stylistic. The linear coefficient is built from the residual, so when
        // the residual is below its own floor the coefficient's SIGN is noise
        // too. Letting a noise-signed quantity answer a structural question would
        // report "no such shape exists" on the strength of a bit that carries no
        // information, which is a confident answer the arithmetic did not earn.
        // Below the floor the honest report is that the request cannot be told
        // apart from one the profile already realizes. Neither shape can reach
        // here with a decisively NEGATIVE residual, because a negative residual
        // puts the triangular velocity below the larger boundary velocity and the
        // reachability gate above has already rejected it -- so a residual that
        // fails its floor here is always noise, never a resolved negative.
        Scalar v{};
        if (selected == shape::plateau && v_hi > Scalar{0}) {
            // The unknown is the cruise velocity's RISE above the larger boundary
            // velocity, solved against the duration DECREMENT the caller asked
            // for: the total duration falls as the cruise velocity grows, and the
            // plateau lies above its boundary. Substituting v = v_hi + e gives
            //     e^2 - B e + C = 0,  B = a (r / v_hi - dT),  C = a v_hi dT
            // whose smaller positive root is the admissible one, taken in the form
            // that adds rather than subtracts. Both coefficients are built
            // directly from the residual and the decrement, so their rounding is
            // proportional to the answer rather than to v_hi^2.
            auto const dT = at_hi.T - T_new;
            if (!ramp_resolved) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            auto const boundary_floor = Scalar{boundary_duration_rounding_ops}
                                        * std::numeric_limits<Scalar>::epsilon() * at_hi.scale;
            if (!(dT > boundary_floor)) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            // The linear coefficient is a DIFFERENCE where the valley's mirror
            // is a sum, and it needs a floor of its own. Both quantities it
            // subtracts grow without bound as the larger boundary velocity
            // vanishes against the commanded displacement -- the residual
            // divided by that velocity, and a decrement measured from a boundary
            // duration which is itself that same quotient plus the two ramps.
            // Deep in the plateau the two agree to the full width of the
            // significand, and the coefficient that locates the root is left
            // with nothing. Solving on it returns a cruise velocity chosen by
            // noise, and the profile built from it realizes a duration the
            // caller never asked for while the call reports success.
            //
            // Six operations form the residual and a seventh divides it by the
            // boundary velocity; fifteen form the boundary duration and a
            // sixteenth subtracts the request from it. Six plus one plus fifteen
            // plus one is twenty-three, each worth one unit in the last place at
            // the larger of the two operands' own scales -- the residual's scale
            // carried through the same division, and the boundary duration's.
            //
            // The magnitude is tested, not the value, so that a coefficient
            // which is decisively negative still reaches the sign gate below and
            // is reported as no-such-shape. Only a coefficient that cannot be
            // told apart from zero is reported as below the arithmetic
            // resolution, which is the same distinction the two floors above
            // draw.
            constexpr int coefficient_rounding_ops = 23;
            auto const coefficient_scale = std::max(residual_scale / v_hi, at_hi.scale);
            auto const coefficient_floor = Scalar{coefficient_rounding_ops}
                                           * std::numeric_limits<Scalar>::epsilon()
                                           * coefficient_scale;
            auto const rise_coefficient = ramp_residual / v_hi - dT;
            if (!(std::abs(rise_coefficient) > coefficient_floor)) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            auto const B = a * rise_coefficient;
            auto const C = a * v_hi * dT;
            if (!(B > Scalar{0})) {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            auto const disc = B * B - Scalar{4} * C;
            // A discriminant that left the representable range says nothing
            // about reachability. The squared linear coefficient overflows
            // whenever that coefficient exceeds the square root of the type's
            // largest value, which a long move at a small cruise velocity
            // reaches long before anything about the request is unreasonable,
            // and an infinite discriminant then drives the root selection's
            // denominator to infinity and its root to zero -- a cruise velocity
            // sitting exactly on the shape boundary, inside its own validity
            // interval, with every phase duration nonnegative. The profile built
            // from it realizes the boundary's duration for a request that asked
            // for something else, and the call reports success. It is the same
            // fact the two resolution floors report: the expressions that would
            // locate the root cannot represent what they need to.
            if (!std::isfinite(disc)) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            if (!(disc >= Scalar{0})) {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            v = v_hi + Scalar{2} * C / (B + std::sqrt(disc));
        } else if (selected == shape::plateau) {
            // The larger boundary velocity is not positive, so it is not an
            // admissible cruise velocity and the shape has no boundary to sit
            // near: every admissible cruise velocity is separated from it by at
            // least its own magnitude, the rise is never small, and the two terms
            // of the discriminant do not approach each other. Solved directly.
            auto const b = (v0_ + v1_) + a * T_new;
            auto const c = a * h + v_sum_sq / Scalar{2};
            auto const disc = b * b - Scalar{4} * c;
            // A discriminant that left the representable range says nothing
            // about reachability. The squared linear coefficient overflows
            // whenever that coefficient exceeds the square root of the type's
            // largest value, which a long move at a small cruise velocity
            // reaches long before anything about the request is unreasonable,
            // and an infinite discriminant then drives the root selection's
            // denominator to infinity and its root to zero -- a cruise velocity
            // sitting exactly on the shape boundary, inside its own validity
            // interval, with every phase duration nonnegative. The profile built
            // from it realizes the boundary's duration for a request that asked
            // for something else, and the call reports success. It is the same
            // fact the two resolution floors report: the expressions that would
            // locate the root cannot represent what they need to.
            if (!std::isfinite(disc)) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
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
            // The mirror of the plateau solve, about the other boundary: the
            // unknown is the cruise velocity's DECREMENT below the smaller
            // boundary velocity, solved against the duration INCREMENT, because
            // the valley lies below its boundary where the plateau lies above its
            // own. Substituting v = v_lo - d gives
            //     d^2 - B d + C = 0,  B = a (dT + r / v_lo),  C = a v_lo dT
            // whose smaller positive root is the admissible one -- the decrement
            // nearest the boundary, which is the larger cruise velocity the
            // cruise-velocity form selected. The cruise velocity is formed once,
            // at the end.
            auto const dT = T_new - at_lo.T;
            if (!ramp_resolved) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            auto const boundary_floor = Scalar{boundary_duration_rounding_ops}
                                        * std::numeric_limits<Scalar>::epsilon() * at_lo.scale;
            if (!(dT > boundary_floor)) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            auto const B = a * (dT + ramp_residual / v_lo);
            auto const C = a * v_lo * dT;
            if (!(B > Scalar{0})) {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            auto const disc = B * B - Scalar{4} * C;
            // A discriminant that left the representable range says nothing
            // about reachability. The squared linear coefficient overflows
            // whenever that coefficient exceeds the square root of the type's
            // largest value, which a long move at a small cruise velocity
            // reaches long before anything about the request is unreasonable,
            // and an infinite discriminant then drives the root selection's
            // denominator to infinity and its root to zero -- a cruise velocity
            // sitting exactly on the shape boundary, inside its own validity
            // interval, with every phase duration nonnegative. The profile built
            // from it realizes the boundary's duration for a request that asked
            // for something else, and the call reports success. It is the same
            // fact the two resolution floors report: the expressions that would
            // locate the root cannot represent what they need to.
            if (!std::isfinite(disc)) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
            if (!(disc >= Scalar{0})) {
                return ctrlpp::unexpected(trajectory_error::unreachable_duration);
            }
            auto const d = Scalar{2} * C / (B + std::sqrt(disc));

            // Carry whichever of the two quantities is the smaller. Forming the
            // larger of them from the boundary velocity and the smaller costs at
            // most one bit; forming the smaller from the boundary velocity and the
            // larger is a cancellation, and it is the same cancellation this whole
            // rearrangement exists to avoid, merely moved to the last step.
            //
            // The two regimes are disjoint and each form is exact in the other's
            // blind spot. The decrement is the smaller one near the shape
            // boundary, which is where the cruise-velocity form's discriminant
            // collapses. The cruise velocity is the smaller one deep in the
            // valley, far below the boundary -- and there the cruise-velocity
            // form's constant term is negative, so its discriminant is a SUM of
            // two nonnegative terms and carries no cancellation at all. Measured
            // over the valley's full range rather than only the boundary
            // neighbourhood, the decrement alone is worse than the cruise-velocity
            // form in the deep valley by many orders of magnitude, so neither form
            // dominates and the choice between them is not cosmetic.
            //
            // The comparison introduces no constant: it asks only which of the
            // decrement and the cruise velocity it leaves behind is larger.
            if (d + d <= v_lo) {
                v = v_lo - d;
            } else {
                auto const b = a * T_new - (v0_ + v1_);
                auto const c = v_sum_sq / Scalar{2} - a * h;
                auto const disc_v = b * b - Scalar{4} * c;
                // A discriminant that left the representable range says nothing
                // about reachability. The squared linear coefficient overflows
                // whenever that coefficient exceeds the square root of the type's
                // largest value, which a long move at a small cruise velocity
                // reaches long before anything about the request is unreasonable,
                // and an infinite discriminant then drives the root selection's
                // denominator to infinity and its root to zero -- a cruise velocity
                // sitting exactly on the shape boundary, inside its own validity
                // interval, with every phase duration nonnegative. The profile built
                // from it realizes the boundary's duration for a request that asked
                // for something else, and the call reports success. It is the same
                // fact the two resolution floors report: the expressions that would
                // locate the root cannot represent what they need to.
                if (!std::isfinite(disc_v)) {
                    return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
                }
                if (!(disc_v >= Scalar{0})) {
                    return ctrlpp::unexpected(trajectory_error::unreachable_duration);
                }
                auto const root = std::sqrt(disc_v);
                // Larger root, in the form that avoids subtracting two nearly equal terms.
                v = (b > Scalar{0}) ? (Scalar{-2} * c / (b + root)) : ((root - b) / Scalar{2});
            }
        }

        // A root is accepted only inside its own shape's validity interval, with
        // every phase duration nonnegative and the cruise velocity within the
        // velocity limit. A root failing any of these is a rejection, never a
        // clamped value.
        //
        // The valley interval is OPEN at its upper end, and the strictness is
        // load-bearing rather than cosmetic. That shape is selected only when the
        // request exceeds the duration the smaller boundary velocity itself
        // realizes, and the total duration is strictly decreasing in the cruise
        // velocity, so the root that answers such a request lies strictly below
        // that boundary. A root landing exactly on it did not solve the equation:
        // it is what the closed form returns when its discriminant, a difference
        // of two nearly equal quantities, cancels to zero and leaves the root with
        // no significant digits at all. The profile built from it would realize
        // the boundary's own duration for every request past it, which is a
        // silently wrong retiming rather than a rejected one. The test costs no
        // constant, because the strictness comes from the branch condition that
        // selected the shape.
        auto const in_shape = (selected == shape::plateau)        ? (v >= v_hi)
                              : (selected == shape::ramp_through) ? (v >= v_lo && v <= v_hi)
                                                                  : (v < v_lo);
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
    trapezoidal_disposition<Scalar> disposition_{};
};

static_assert(trajectory_segment<trapezoidal_trajectory<double>, double, 1>);

}

#endif
