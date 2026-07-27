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
///
/// Construction goes through `create`, which returns
/// `ctrlpp::expected<double_s_trajectory, trajectory_error>` and is the only way
/// to obtain one. There is no non-validating public constructor: a profile
/// object either satisfies its contract -- a finite, nonnegative duration in
/// every segment and a command its seven segments actually traverse -- or it was
/// never built. Every rejection therefore reaches the caller as a value it has
/// to inspect, and `evaluate` can rely on the invariants that keep its own
/// arithmetic defined.
///
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

    /// @brief Construct a profile, or report why the command has none.
    ///
    /// The only construction path. Rejections, checked in order:
    ///  * a NaN or infinite position or boundary velocity ->
    ///    trajectory_error::non_finite_input
    ///  * a NaN, infinite, or non-positive velocity limit ->
    ///    trajectory_error::non_positive_velocity_limit
    ///  * a NaN, infinite, or non-positive acceleration limit ->
    ///    trajectory_error::non_positive_acceleration_limit
    ///  * a NaN, infinite, or non-positive jerk limit ->
    ///    trajectory_error::non_positive_jerk_limit
    ///  * a boundary velocity above the velocity limit in magnitude ->
    ///    trajectory_error::boundary_velocity_exceeds_limit
    ///  * a commanded displacement below the distance the fastest admissible
    ///    transition between the two boundary velocities already sweeps, of
    ///    which a zero displacement under a nonzero boundary velocity is the
    ///    extreme case -> trajectory_error::unreachable_boundary_velocity
    ///  * a negative segment duration or total duration ->
    ///    trajectory_error::unreachable_boundary_velocity
    ///  * a duration outside the representable range, or a total duration that
    ///    underflowed to zero on a nonzero commanded displacement ->
    ///    trajectory_error::unrepresentable_duration
    ///
    /// The velocity limit is a PRECONDITION on the boundary velocities, not a
    /// bound the construction raises to accommodate them. Raising it would
    /// return a profile that violates a limit the caller stated, which is the
    /// same silent-success defect a rejection exists to prevent.
    ///
    /// A seven-segment profile cannot sweep less ground than the fastest
    /// admissible transition from the larger of the two positive-frame boundary
    /// velocities to the smaller one: covering less would require overshooting
    /// the target and coming back, which is a different velocity profile shape.
    /// That minimum is the swept distance of the cruise-free profile whose peak
    /// rises by nothing above the larger boundary velocity, so it is computed
    /// from the same ramp expressions the profile is built out of rather than
    /// from a separate approximation of them. A zero commanded displacement is
    /// realizable by exactly one member of the family, the standstill, and only
    /// when both boundary velocities are zero as well.
    ///
    /// @cite biagiotti2009 -- Sec. 3.4.1, p.79-85
    [[nodiscard]] static auto create(config const& cfg)
        -> ctrlpp::expected<double_s_trajectory, trajectory_error>
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
        if (!std::isfinite(cfg.j_max) || cfg.j_max <= Scalar{0}) {
            return ctrlpp::unexpected(trajectory_error::non_positive_jerk_limit);
        }
        if (std::abs(cfg.v0) > cfg.v_max || std::abs(cfg.v1) > cfg.v_max) {
            return ctrlpp::unexpected(trajectory_error::boundary_velocity_exceeds_limit);
        }

        auto const h_signed = cfg.q1 - cfg.q0;
        if (h_signed == Scalar{0}) {
            if (cfg.v0 != Scalar{0} || cfg.v1 != Scalar{0}) {
                return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
            }
        } else {
            auto const sigma = (h_signed > Scalar{0}) ? Scalar{1} : Scalar{-1};
            auto const pv0 = sigma * cfg.v0;
            auto const pv1 = sigma * cfg.v1;
            auto const h_min = no_cruise_displacement(
                Scalar{0}, std::max(pv0, pv1), std::min(pv0, pv1), cfg.a_max, cfg.j_max);
            if (sigma * h_signed < h_min) {
                return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
            }
        }

        double_s_trajectory profile{unchecked_t{}, cfg};

        // The solved durations are checked as the solver stored them, not as
        // phase_durations() reports them: that accessor subtracts the two jerk
        // sub-segments out of each ramp and floors the remainder at zero to
        // absorb the rounding of a difference that is algebraically exact, and a
        // ramp duration that came out negative would disappear into that floor.
        for (auto const stored : std::array<Scalar, 6>{profile.T_j1_, profile.T_a_, profile.T_v_,
                                                       profile.T_d_, profile.T_j2_, profile.T_}) {
            if (stored < Scalar{0}) {
                return ctrlpp::unexpected(trajectory_error::unreachable_boundary_velocity);
            }
            if (!std::isfinite(stored)) {
                return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
            }
        }
        // A nonzero commanded displacement is not traversed in zero time at any
        // finite velocity, so a total duration that underflowed to zero under
        // one describes no motion the command asked for. The predicate is the
        // physical statement itself and carries no threshold.
        if (h_signed != Scalar{0} && !(profile.T_ > Scalar{0})) {
            return ctrlpp::unexpected(trajectory_error::unrepresentable_duration);
        }
        return profile;
    }

    /// @brief Evaluate trajectory at time t.
    ///
    /// The clamp below is well defined because no object of this type exists
    /// with a negative duration: `create` is the only construction path and it
    /// rejects one. `std::clamp` with a lower bound above its upper bound is
    /// undefined behavior, not a clamp that returns something unhelpful, so that
    /// construction-time guarantee is what keeps this line safe rather than a
    /// test performed here on every evaluation.
    ///
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

    /// @brief Rescale the profile to a longer duration for multi-axis synchronization.
    ///
    /// The profile is rebuilt under scaled kinematic limits, never patched. The
    /// velocity, acceleration, and jerk limits are multiplied by the first,
    /// second, and third power of one scale factor and the profile is constructed
    /// again from the same command; the two boundary velocities are left UNSCALED,
    /// because they are what the caller commanded the axis to enter and leave with
    /// and scaling them would land a synchronized axis at the wrong terminal
    /// velocity. Displacement, terminal velocity, and continuity then hold by
    /// construction rather than by repair. The stored duration stays whatever the
    /// rebuilt profile realizes and is never assigned the requested value; it lands
    /// within a few units in the last place of it.
    ///
    /// Rejections, checked in order:
    ///  * a duration below the current one -> trajectory_error::duration_shorter_than_current
    ///  * NaN, infinite, or non-positive T_new -> trajectory_error::non_positive_duration
    ///  * a duration no admissible scale realizes -> trajectory_error::unreachable_duration
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
        *this = solved.value();
        return {};
    }

    /// @brief Report whether rescale_to(T_new) would succeed, without mutating.
    ///
    /// Runs the identical solve rescale_to() runs and discards the result, so the
    /// two cannot disagree: identical inputs traverse identical code with no
    /// intervening state. That is what lets a multi-axis synchronization check
    /// every axis before it commits any of them. No closed-form reachability
    /// predicate exists for this family, so structural replay is the only way to
    /// make the check and the commit agree.
    [[nodiscard]] auto can_rescale_to(Scalar T_new) const -> ctrlpp::expected<void, trajectory_error>
    {
        auto const solved = solve_rescale(T_new);
        if (!solved.has_value()) {
            return ctrlpp::unexpected(solved.error());
        }
        return {};
    }

    /// @brief Phase durations for the 7 segments.
    ///
    /// Returns {T_j1, T_a - 2*T_j1, T_j1, T_v, T_j2, T_d - 2*T_j2, T_j2}
    /// representing jerk(+), const-accel, jerk(-), cruise, jerk(-), const-decel, jerk(+).
    ///
    /// The floor on the two constant-acceleration remainders is a rounding
    /// guard, not a stand-in for a rejected command. A ramp is triangular in
    /// acceleration exactly when its remainder is algebraically zero, and it
    /// carries a constant-acceleration segment exactly when the remainder is
    /// algebraically positive; the floor only absorbs the rounding of that
    /// difference at the boundary between the two. A ramp duration that came out
    /// negative is a different matter, and `create` rejects the command rather
    /// than letting this accessor hide it.
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

    /// @brief Tag selecting the non-validating constructor reserved for `create`.
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    /// @brief Solve the profile from a configuration `create` has checked.
    ///
    /// Non-validating by construction: it follows the B&M flowchart and stores
    /// what comes out, and it is private so the only way to reach it is through
    /// `create`, which decides both before and after whether what came out is a
    /// profile.
    ///
    /// The zero-displacement branch is keyed on exact equality rather than on a
    /// neighborhood of zero. A displacement small enough to underflow the rest
    /// of the algebra is a rejection `create` makes, not a command this
    /// constructor rounds down to a standstill; the two are different answers
    /// and only one of them is true.
    ///
    /// @cite biagiotti2009 -- Sec. 3.4, eq. (3.17)-(3.33), Fig. 3.18, p.79-96
    explicit double_s_trajectory(unchecked_t, config const& cfg)
        : q0_{cfg.q0}
        , q1_{cfg.q1}
        , v0_{cfg.v0}
        , v1_{cfg.v1}
        , v_max_{cfg.v_max}
        , a_max_{cfg.a_max}
        , j_max_{cfg.j_max}
    {
        auto const h_signed = cfg.q1 - cfg.q0;

        // Zero displacement: the standstill, which `create` admits only when
        // both boundary velocities are zero as well.
        if (h_signed == Scalar{0}) {
            sigma_ = Scalar{1};
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

    auto make_point(Scalar q, Scalar dq, Scalar ddq) const -> trajectory_point<Scalar, 1>
    {
        return {.position = Vector<Scalar, 1>{q},
                .velocity = Vector<Scalar, 1>{dq},
                .acceleration = Vector<Scalar, 1>{ddq}};
    }

    /// @brief Rebuild the profile with every kinematic limit scaled, or report that
    /// the rebuilt profile does not respect its own scaled limits.
    ///
    /// A time scaling that slows a profile by a factor divides its velocity by that
    /// factor, its acceleration by its square, and its jerk by its cube, which is
    /// why the three limits carry the first, second, and third power of the scale.
    /// The boundary velocities are deliberately NOT scaled: they are the command,
    /// not a limit.
    ///
    /// The rebuild goes back through `create`, so the scaled command is held to
    /// the same contract the original was and a scale that produces no profile
    /// is reported rather than applied. The rebuilt profile is admissible only
    /// when it also has a positive duration and peaks no higher than its own
    /// scaled velocity limit. The latter is what the scale's lower bound below
    /// expresses: at a scale where the commanded boundary velocities themselves
    /// reach the scaled velocity limit, no seven-segment shape can hold the peak
    /// underneath it.
    ///
    /// @cite biagiotti2009 -- Sec. 5.3, eq. (5.13)-(5.14) -- time scaling of a profile
    [[nodiscard]] static auto rebuild_scaled(config const& cfg, Scalar lambda)
        -> ctrlpp::expected<double_s_trajectory, trajectory_error>
    {
        auto scaled = cfg;
        scaled.v_max = cfg.v_max * lambda;
        scaled.a_max = cfg.a_max * lambda * lambda;
        scaled.j_max = cfg.j_max * lambda * lambda * lambda;

        auto rebuilt = create(scaled);
        if (!rebuilt.has_value()) {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        auto const& profile = rebuilt.value();
        if (!(profile.T_ > Scalar{0}) || !(profile.v_lim_ <= scaled.v_max)) {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        return rebuilt;
    }

    /// @brief Whether the scale is at or below the one that realizes T_new.
    ///
    /// A scale that yields no admissible profile lies below the admissible range
    /// and answers true, and inside that range the duration falls as the scale
    /// grows, so the answer is true on an interval reaching down from the crossing
    /// point and false above it. That is what makes bracket halving valid here.
    [[nodiscard]] static auto reaches_duration(config const& cfg, Scalar lambda, Scalar T_new) -> bool
    {
        auto const rebuilt = rebuild_scaled(cfg, lambda);
        return !rebuilt.has_value() || rebuilt.value().T_ >= T_new;
    }

    /// @brief Solve the rescaled profile, or report why the request is not realizable.
    ///
    /// This is the single routine behind both rescale_to() and can_rescale_to(),
    /// so the check and the commit replay the identical deterministic computation.
    ///
    /// With both boundary velocities at rest the duration is exactly proportional
    /// to the reciprocal of the scale, so the scale is the ratio of the current
    /// duration to the requested one and one rebuild settles it. With a nonzero
    /// boundary velocity that proportionality fails -- the duration is a cubic in
    /// the reciprocal of the scale within a fixed segment shape, and it carries
    /// real kinks where the shape flips -- so the scale is found by halving a
    /// bracket instead. A derivative step is not used: a kink can throw it out of
    /// the bracket, and bounding a safeguarded variant would need an iteration cap.
    ///
    /// The bracket runs from the scale at which the commanded boundary velocities
    /// themselves reach the scaled velocity limit up to the profile's own scale.
    /// Termination is bracket exhaustion, when the midpoint falls on an endpoint:
    /// no tolerance, no trip count. The worst case is derived rather than chosen.
    /// Halving an interval whose endpoints share a binary exponent reaches the
    /// spacing of the representable values after one step more than the significand
    /// width, which is 25 evaluations at single precision and 54 at double
    /// precision; a bracket spanning several exponents costs one further step per
    /// exponent spanned. Early exit is not a determinism problem: identical inputs
    /// exhaust the bracket at the identical step.
    ///
    /// @cite biagiotti2009 -- Sec. 5.3 -- time scaling for multi-axis synchronization
    [[nodiscard]] auto solve_rescale(Scalar T_new) const
        -> ctrlpp::expected<double_s_trajectory, trajectory_error>
    {
        // Exact equality first: a synchronized set passes the slowest axis a
        // bit-exact copy of its own duration, so that axis lands here rather than
        // in the shortening rejection below. There is no float-equality fragility
        // in that -- the value compared is the same object's own stored duration.
        if (T_new == T_) {
            return *this;
        }
        if (!std::isfinite(T_new) || T_new <= Scalar{0}) {
            return ctrlpp::unexpected(trajectory_error::non_positive_duration);
        }
        if (T_new < T_) {
            return ctrlpp::unexpected(trajectory_error::duration_shorter_than_current);
        }
        if (!(T_ > Scalar{0})) {
            // A standstill traverses nothing, so no scale gives it a positive
            // duration.
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        auto const cfg = config{.q0 = q0_,
                                .q1 = q1_,
                                .v_max = v_max_,
                                .a_max = a_max_,
                                .j_max = j_max_,
                                .v0 = v0_,
                                .v1 = v1_};

        if (v0_ == Scalar{0} && v1_ == Scalar{0}) {
            return rebuild_scaled(cfg, T_ / T_new);
        }

        auto const lambda_min = std::max(std::abs(v0_), std::abs(v1_)) / v_max_;
        if (!(lambda_min < Scalar{1})) {
            // A boundary velocity already reaches the velocity limit, so there is
            // no room left to slow the profile down.
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }

        auto lo = lambda_min;
        auto hi = Scalar{1};
        for (;;) {
            auto const mid = lo + (hi - lo) / Scalar{2};
            if (!(mid > lo) || !(mid < hi)) {
                break;
            }
            if (reaches_duration(cfg, mid, T_new)) {
                lo = mid;
            } else {
                hi = mid;
            }
        }

        // The bracket is exhausted, so lo and hi are neighboring values with the
        // crossing between them. Accepting lo requires it to be admissible and to
        // reach the requested duration; when it is not, the crossing lies inside
        // the range no admissible profile covers and the request is unreachable.
        auto const solved = rebuild_scaled(cfg, lo);
        if (!solved.has_value() || !(solved.value().T_ >= T_new)) {
            return ctrlpp::unexpected(trajectory_error::unreachable_duration);
        }
        return solved;
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
        // shape; `create` decides that before this constructor runs, so the
        // solve below can assume its bracket contains the root instead of
        // producing a profile that does not traverse its own command.
        //
        // The unknown carried through this step is the peak's RISE above the
        // larger boundary velocity, never the peak itself. The rise is what every
        // ramp duration depends on, and it is routinely orders of magnitude
        // smaller than the boundary velocity it sits on, so storing the peak and
        // subtracting the boundary velocity back out of it would throw away most
        // of the rise's significant digits before the ramps ever see it.

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
