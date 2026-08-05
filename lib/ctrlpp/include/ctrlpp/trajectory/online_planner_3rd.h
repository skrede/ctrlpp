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
/// and do NOT satisfy trajectory_segment. They are stateful filters.
///
/// @cite biagiotti2009 -- Biagiotti & Melchiorri, "Trajectory Planning for
/// Automatic Machines and Robots", 2009, Sec. 4.6.1
/// @cite lambrechts2005 -- Lambrechts, Boerlage & Steinbuch, "Trajectory Planning and Feedforward Design for Electromechanical Motion Systems", Control Engineering Practice 13(2), 2005 (jerk-limited online planning)

#include "ctrlpp/expected.h"

#include "ctrlpp/trajectory/trajectory_types.h"
#include "ctrlpp/trajectory/double_s_trajectory.h"
#include "ctrlpp/trajectory/online_planner_diagnostics.h"

#include "ctrlpp/util/concepts.h"

#include <cmath>
#include <limits>
#include <cstddef>
#include <algorithm>

namespace ctrlpp
{

/// @brief 3rd-order online planner producing double-S (jerk-limited) velocity profiles.
///
/// Generates bounded-velocity, bounded-acceleration, bounded-jerk trajectories
/// that can be replanned mid-motion when a new target arrives.
///
/// Construction goes through `create`, which validates the kinematic
/// limits and reports rejections through
/// `ctrlpp::expected<online_planner_3rd, trajectory_error>`.
///
/// @cite biagiotti2009 -- Sec. 4.6.1
template <ctrlpp_floating_scalar Scalar>
class online_planner_3rd
{
  public:
    /// @brief Kinematic limits, and the distances at which a sampled state
    /// counts as arrived.
    ///
    /// The three settle tolerances are policy, not rounding guards, and they are
    /// separate fields because they are compared against three quantities in
    /// three different units: a length, a speed and an acceleration. Their
    /// defaults sit roughly seven decades above `double` rounding on the
    /// quantities they test, so they state when the application considers the
    /// axis arrived rather than when the arithmetic can still tell two numbers
    /// apart. No derivation is offered for them because none exists: the
    /// distance at which a machine is "there" is a property of the machine.
    ///
    /// The defaults are the values the planner used before the fields existed,
    /// so omitting them reproduces that behavior exactly for every choice of
    /// limits. They are chosen for `double`; a `float` instantiation should set
    /// them, because the defaults sit below `float`'s own resolution near unity.
    ///
    /// No setting of them changes what the planner computes. They decide only
    /// when `is_settled()` reports true; the profile, its phase durations and
    /// every value `sample()` returns are identical under any value.
    struct config
    {
        Scalar v_max;
        Scalar a_max;
        Scalar j_max;

        /// @brief Position error below which the axis counts as arrived.
        Scalar position_settle_tol = static_cast<Scalar>(1e-9);

        /// @brief Speed below which the axis counts as stopped.
        Scalar velocity_settle_tol = static_cast<Scalar>(1e-9);

        /// @brief Acceleration below which the axis counts as at rest.
        Scalar acceleration_settle_tol = static_cast<Scalar>(1e-9);
    };

    /// @brief Validate the kinematic limits and construct a planner.
    ///
    /// All three limits are divisors in the planner math (cruise duration
    /// h / v_max, jerk-phase durations a_max / j_max and |a| / j_max, and the
    /// a_max-reached threshold a_max^2 / j_max), so the exact mathematical
    /// domain of each is finite and strictly positive. Rejections, checked in
    /// order:
    ///  * NaN/Inf or non-positive v_max -> trajectory_error::non_positive_velocity_limit
    ///  * NaN/Inf or non-positive a_max -> trajectory_error::non_positive_acceleration_limit
    ///  * NaN/Inf or non-positive j_max -> trajectory_error::non_positive_jerk_limit
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.1
    static auto create(config const& cfg)
        -> ctrlpp::expected<online_planner_3rd, trajectory_error>
    {
        if (!std::isfinite(cfg.v_max) || cfg.v_max <= Scalar{0}) {
            return ctrlpp::unexpected(trajectory_error::non_positive_velocity_limit);
        }
        if (!std::isfinite(cfg.a_max) || cfg.a_max <= Scalar{0}) {
            return ctrlpp::unexpected(trajectory_error::non_positive_acceleration_limit);
        }
        if (!std::isfinite(cfg.j_max) || cfg.j_max <= Scalar{0}) {
            return ctrlpp::unexpected(trajectory_error::non_positive_jerk_limit);
        }
        return online_planner_3rd{unchecked_t{}, cfg};
    }

    /// @brief Set new target position. Replans from current state.
    ///
    /// Computes time-optimal double-S profile from (q_, v_, a_) to (target, 0, 0)
    /// respecting v_max, a_max, and j_max. A same-direction move carries the
    /// current velocity through the profile; a velocity pointing away from the
    /// target (or too large to stop in the available distance) is braked to rest
    /// first, then replanned. The carry-velocity shape has a domain of its own,
    /// and a commanded state outside it is braked to rest and replanned as well.
    ///
    /// A non-finite target is rejected before any member changes. For a finite
    /// target, the commanded shape is not always the one realized: which
    /// profile was built is read back from `diagnostics()`, where `disposition`
    /// names the branch taken and
    /// `substitution_reason` names the condition that selected it. The motion
    /// respects every limit either way; what changes is the time it takes, which
    /// `planned_duration` and `brake_duration` are there to account for.
    ///
    /// @cite biagiotti2009 -- Sec. 4.6.1
    auto update(Scalar target) -> ctrlpp::expected<void, trajectory_error>
    {
        if(!std::isfinite(target))
            return ctrlpp::unexpected(trajectory_error::non_finite_input);

        target_ = target;
        t_ref_ = t_last_;

        // Snapshot current state
        q_ref_ = q_;
        v_ref_ = v_;
        a_ref_ = a_;

        compute_profile();
        return {};
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

        // Check settled state. Each dimension is tested against its own
        // configured distance, because the three quantities carry three
        // different units.
        settled_ = (std::abs(q - target_) < position_settle_tol_)
                   && (std::abs(v) < velocity_settle_tol_)
                   && (std::abs(a) < acceleration_settle_tol_);

        return {
            .position = Vector<Scalar, 1>{q},
            .velocity = Vector<Scalar, 1>{v},
            .acceleration = Vector<Scalar, 1>{a},
        };
    }

    /// @brief True when at target with zero velocity and zero acceleration.
    auto is_settled() const -> bool { return settled_; }

    /// @brief What the last update planned, against what it was commanded.
    ///
    /// Describes the profile built by the last `update` or `reset`, not the
    /// state reached since: `is_settled()` answers that.
    auto diagnostics() const -> online_planner_diagnostics<Scalar> const& { return diagnostics_; }

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
        diagnostics_ = online_planner_diagnostics<Scalar>{
            .disposition = online_planner_disposition::settled,
            .substitution_reason = online_planner_substitution_reason::none,
            .commanded_target = q0,
            .initial_velocity = Scalar{0},
            .planned_duration = Scalar{0},
            .brake_duration = Scalar{0},
            .replan_start_position = q0,
        };
    }

  private:
    /// @brief Tag selecting the non-validating constructor reserved for `create`.
    struct unchecked_t
    {
        explicit unchecked_t() = default;
    };

    /// @brief Construct from a configuration already validated by `create`.
    /// Initial state at rest at q=0.
    online_planner_3rd(unchecked_t, config const& cfg)
        : v_max_{cfg.v_max}
        , a_max_{cfg.a_max}
        , j_max_{cfg.j_max}
        , position_settle_tol_{cfg.position_settle_tol}
        , velocity_settle_tol_{cfg.velocity_settle_tol}
        , acceleration_settle_tol_{cfg.acceleration_settle_tol}
    {
    }

    Scalar v_max_{};
    Scalar a_max_{};
    Scalar j_max_{};

    // Settle policy, read only by sample(). Never a divisor and never a loop
    // bound, so no value of any of them can make the planner run unboundedly.
    Scalar position_settle_tol_{};
    Scalar velocity_settle_tol_{};
    Scalar acceleration_settle_tol_{};

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

    // Disposition of the profile built at the last update. Assigned as a whole
    // aggregate on every branch that completes a plan, which is what keeps a
    // field from one branch surviving into the report of another.
    online_planner_diagnostics<Scalar> diagnostics_{};

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
        n_phases_ = 0;

        // Is the commanded state a numerical no-op -- a command the arithmetic
        // cannot tell apart from the one already in force?
        //
        // This is NOT the settle policy the caller configures. That policy
        // answers "is the motion done", which is a statement about the machine
        // and has no derivation. This answers "is there anything left to
        // compute", which is a statement about the arithmetic and has one. The
        // two are deliberately separate, and the window between them is wide:
        // on a unit-scale axis the policy default sits some seven decades above
        // the floors below. Letting a policy tolerance decide what the planner
        // is allowed to compute would collapse them into one question, and they
        // are not one question.
        //
        // The position clause stays an exact comparison. A commanded target one
        // unit in the last place away from the reference position is a
        // different target, and no rounding argument says otherwise: both are
        // inputs, and neither is the result of a chain this planner ran.
        //
        // The velocity and acceleration clauses compare a caller-supplied
        // snapshot against the resolution of the limit that bounds it. This
        // function performs no arithmetic on either before the test, so the
        // honest local count is one operation each, the comparison itself.
        // Whatever error the caller's own sampling carried into the snapshot is
        // the caller's, and is out of scope here -- the same scoping the
        // trapezoidal retiming floor applies to its own.
        //
        // Below these floors the commanded velocity and acceleration are zero
        // to the precision the limits leave available, so there is nothing for
        // a profile to remove.
        constexpr int command_velocity_rounding_ops = 1;
        constexpr int command_acceleration_rounding_ops = 1;
        auto const command_velocity_resolution = Scalar{command_velocity_rounding_ops}
                                                 * std::numeric_limits<Scalar>::epsilon()
                                                 * v_max_;
        auto const command_acceleration_resolution = Scalar{command_acceleration_rounding_ops}
                                                     * std::numeric_limits<Scalar>::epsilon()
                                                     * a_max_;

        // Check if already settled
        if (target_ == q_ref_
            && std::abs(v_ref_) < command_velocity_resolution
            && std::abs(a_ref_) < command_acceleration_resolution) {
            T_ = Scalar{0};
            settled_ = true;
            diagnostics_ = online_planner_diagnostics<Scalar>{
                .disposition = online_planner_disposition::settled,
                .substitution_reason = online_planner_substitution_reason::none,
                .commanded_target = target_,
                .initial_velocity = v_ref_,
                .planned_duration = Scalar{0},
                .brake_duration = Scalar{0},
                .replan_start_position = q_ref_,
            };
            return;
        }

        settled_ = false;

        Scalar q_start = q_ref_;
        Scalar v_start = v_ref_;
        Scalar a_start = a_ref_;

        // Phase 0: bring the acceleration to zero, when doing so changes
        // anything.
        //
        // The question is not whether the starting acceleration is nonzero --
        // that comparison has no scale to be against, and a floor placed
        // directly on the acceleration is an absolute length borrowed from an
        // axis nobody named. The question is whether the phase that nulls it
        // moves the axis by anything the arithmetic can still see. The phase
        // lasts |a_start| / j_max, and over that time it
        //   * changes the velocity by exactly a_start^2 / (2 j_max), and
        //   * moves the axis by at most v_max * |a_start| / j_max, which is the
        //     leading term v_start * T_az with the speed replaced by the only
        //     bound the planner has for it.
        // Both follow from the limits alone, which is what makes them the
        // better question.
        //
        // Velocity side: three operations form the contribution -- the square,
        // the doubled jerk limit and the division -- compared against the
        // resolution of the velocity limit. Length side: two operations form
        // the bound -- the product and the division -- compared against the
        // resolution of the stopping distance from full speed, which three more
        // form: the squared velocity limit, the doubled acceleration limit and
        // the division. That distance is the right length scale here for the
        // same reason it is below: it is intrinsic to the two limits the
        // planner was given and needs no knowledge of the sample period, which
        // the planner is never told.
        //
        // The phase is emitted unless BOTH contributions are unresolvable.
        // Skipping it when only the velocity contribution vanishes would leave
        // the axis a resolvable distance from where the plan assumes it is: the
        // length contribution is linear in the starting acceleration where the
        // velocity contribution is quadratic, so on any ordinary axis the
        // length side binds first, by several decades. Below both, the phase
        // moves nothing that survives being written down.
        constexpr int nulling_velocity_rounding_ops = 3;
        constexpr int nulling_length_rounding_ops = 2;
        constexpr int nulling_scale_rounding_ops = 3;
        constexpr int nulling_displacement_rounding_ops = nulling_length_rounding_ops
                                                          + nulling_scale_rounding_ops;

        auto const abs_a_start = std::abs(a_start);
        auto const nulling_velocity_change = abs_a_start * abs_a_start / (Scalar{2} * j_max_);
        auto const nulling_displacement_bound = v_max_ * abs_a_start / j_max_;
        auto const nulling_length_scale = v_max_ * v_max_ / (Scalar{2} * a_max_);

        bool const nulling_moves_the_velocity =
            nulling_velocity_change > Scalar{nulling_velocity_rounding_ops}
                                          * std::numeric_limits<Scalar>::epsilon() * v_max_;
        bool const nulling_moves_the_axis =
            nulling_displacement_bound > Scalar{nulling_displacement_rounding_ops}
                                             * std::numeric_limits<Scalar>::epsilon()
                                             * nulling_length_scale;

        if (nulling_moves_the_velocity || nulling_moves_the_axis) {
            auto const T_az = abs_a_start / j_max_;
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
        auto const outcome = plan_from_zero_accel(q_start, v_start);

        // Compute total duration
        T_ = Scalar{0};
        for (int i = 0; i < n_phases_; ++i) {
            T_ += T_ph_[i];
        }

        if (outcome.reason == online_planner_substitution_reason::none) {
            diagnostics_ = online_planner_diagnostics<Scalar>{
                .disposition = online_planner_disposition::commanded_profile,
                .substitution_reason = online_planner_substitution_reason::none,
                .commanded_target = target_,
                .initial_velocity = v_ref_,
                .planned_duration = T_,
                .brake_duration = Scalar{0},
                .replan_start_position = q_ref_,
            };
        } else {
            diagnostics_ = online_planner_diagnostics<Scalar>{
                .disposition = online_planner_disposition::braked_and_replanned,
                .substitution_reason = outcome.reason,
                .commanded_target = target_,
                .initial_velocity = v_ref_,
                .planned_duration = T_,
                .brake_duration = outcome.brake_duration,
                .replan_start_position = outcome.replan_start_position,
            };
        }
    }

    /// @brief What `plan_from_zero_accel` did with the commanded shape.
    ///
    /// A reason of `none` leaves the two quantities unset: nothing was braked,
    /// and the plan starts where the command did.
    struct substitution_outcome
    {
        online_planner_substitution_reason reason{online_planner_substitution_reason::none};
        Scalar brake_duration{};
        Scalar replan_start_position{};
    };

    /// @brief Plan a profile from (q0, v0, a=0) to (target_, 0, 0).
    ///
    /// A same-direction move with room to stop carries v0 through the profile via
    /// the general nonzero-initial-velocity double-S. If v0 points away from the
    /// target or is too large to stop within the available distance, the planner
    /// first brakes to rest (jerk-limited) and then plans rest-to-rest.
    ///
    /// Reports which of the two it did, because the caller cannot tell them apart
    /// from the motion alone: both respect every limit and both reach the target.
    auto plan_from_zero_accel(Scalar q0, Scalar v0) -> substitution_outcome
    {
        auto const h_signed = target_ - q0;

        // The four decisions below ask four different questions about
        // quantities in two different units, and one number cannot be the right
        // answer to all of them. Each gets the scale of the quantity it tests
        // and a count of the operations that formed it.

        // Speed floor. The scale is the velocity limit: it is the only velocity
        // the planner is given, and every velocity the profile carries is
        // bounded by it. The count is the roundings this planner performed on
        // the speed being tested, taken over the branches that reach here and
        // maximized. On the branch that plans straight from the caller's
        // snapshot it performed none, and the snapshot's own error is the
        // caller's, out of scope exactly as the trapezoidal retiming floor
        // scopes its own. On the branch that nulled a starting acceleration
        // first, seven operations formed the speed: the phase duration's
        // division, the product of the starting acceleration with it, the three
        // products of the jerk term, and the two sums. The scaling by one half
        // inside that term is exact in a binary radix and is counted anyway,
        // which makes the total an upper bound rather than an estimate. Below
        // the floor the commanded state is at rest to the precision available.
        constexpr int speed_rounding_ops = 7;
        auto const speed_floor = Scalar{speed_rounding_ops}
                                 * std::numeric_limits<Scalar>::epsilon() * v_max_;

        // Length floor. The scale is the planner's own stopping distance from
        // full speed, v_max^2 / (2 a_max). That is the right length because it
        // is intrinsic to the two limits the planner was given and requires no
        // knowledge of the sample period, which the planner is never told. It
        // is also the largest operand that enters the comparison rather than
        // the cancelled result of it: the commanded displacement is a
        // difference of two positions and inherits their scale, not its own.
        //
        // The count breaks out as twelve for the position the
        // acceleration-nulling phase leaves behind -- the phase duration's
        // division, the two operations of the speed term, the four of the
        // squared term and the five of the cubed term -- plus one for the
        // subtraction that forms the commanded displacement from it, plus three
        // for the scale itself: the squared velocity limit, the doubled
        // acceleration limit and the division. On the branch that nulled no
        // acceleration the first term is zero, so the total is again an upper
        // bound. Below the floor the commanded displacement is zero to the
        // precision the limits leave available.
        constexpr int nulling_position_rounding_ops = 12;
        constexpr int displacement_rounding_ops = 1;
        constexpr int length_scale_rounding_ops = 3;
        constexpr int length_rounding_ops = nulling_position_rounding_ops
                                            + displacement_rounding_ops
                                            + length_scale_rounding_ops;
        auto const length_scale = v_max_ * v_max_ / (Scalar{2} * a_max_);
        auto const length_floor = Scalar{length_rounding_ops}
                                  * std::numeric_limits<Scalar>::epsilon() * length_scale;

        // If velocity is zero (or nearly), plan rest-to-rest directly
        if (std::abs(v0) < speed_floor) {
            plan_rest_to_rest(q0);
            return {};
        }

        // Compute stopping distance: distance to bring v0 to 0 using a_max, j_max
        auto const stop_info = compute_stop(v0);
        auto const stop_dist = stop_info.dist;

        // Direction. A length is compared against the length floor and a speed
        // against the speed floor, and the two never meet inside one
        // expression: a boolean that holds a distance and a speed to the same
        // number is a statement about neither.
        bool const wrong_way = (h_signed > length_floor && v0 < -speed_floor)
                               || (h_signed < -length_floor && v0 > speed_floor)
                               || (std::abs(h_signed) < length_floor);

        // Overshoot, compared RELATIVELY. Both sides are lengths this planner
        // has already computed, so the comparison needs no external scale at
        // all and its slack is the rounding the two chains carry rather than a
        // distance taken from somewhere else.
        //
        // The count is the sum along the two chains. The stopping distance's
        // deeper chain is the one that reaches the acceleration limit: the
        // jerk-phase duration is one, the constant-phase duration two, and the
        // three-phase distance twenty-nine -- four for the speed after the
        // first ramp, six for its distance, two for the speed after the
        // constant stretch, six for its distance and eleven for the final ramp
        // -- plus the sign multiplication, which is exact and is counted anyway:
        // thirty-three. The shorter chain, where the acceleration limit is not
        // reached, is four: a division, a square root and two products. The
        // deeper branch is the maximum and is the one written down. The
        // remaining distance carries the twelve of the nulled-acceleration
        // position plus one for the subtraction, and one alone on the branch
        // that nulled nothing.
        constexpr int stopping_distance_chain_rounding_ops = 33;
        constexpr int remaining_distance_chain_rounding_ops = nulling_position_rounding_ops
                                                              + displacement_rounding_ops;
        constexpr int overshoot_rounding_ops = stopping_distance_chain_rounding_ops
                                               + remaining_distance_chain_rounding_ops;
        auto const overshoot_slack = Scalar{overshoot_rounding_ops}
                                     * std::numeric_limits<Scalar>::epsilon();
        bool const overshoot = !wrong_way
                               && (std::abs(stop_dist)
                                   > std::abs(h_signed) * (Scalar{1} + overshoot_slack));

        if (wrong_way || overshoot) {
            // Velocity points away from the target or is too large to stop in the
            // available distance: the only feasible (and time-optimal) option is to
            // brake to rest first, then plan rest-to-rest from the stopping point.
            append_brake_phases(v0, stop_info);
            auto const q_after = q0 + stop_dist;
            plan_rest_to_rest(q_after);
            return {
                .reason = online_planner_substitution_reason::reversal_or_overshoot,
                .brake_duration = stop_info.T_a,
                .replan_start_position = q_after,
            };
        }

        if (append_incorporate_velocity(q0, v0)) {
            // Same-direction move with room to spare: carry the current velocity
            // through the profile rather than braking to rest first. The general
            // nonzero-initial-velocity double-S accelerates from v0 toward the
            // cruise velocity and decelerates to rest at the target, so the move is
            // time-optimal with no full-stop dip.
            return {};
        }

        // That shape has a domain of its own, and the tests above are a
        // planner-side approximation of it rather than a statement of it: the
        // current speed may sit above the velocity limit after a limit change,
        // and the remaining distance may fall below what the transition from
        // the current velocity to rest already sweeps. Where the shape does
        // not exist the planner brakes to rest and replans from the stopping
        // point, which is the same always-admissible fallback the two tests
        // above select and reaches the same target under the same limits. It
        // costs time, not correctness, and nothing is emitted from a profile
        // that was not built.
        //
        // The gap between the approximation and the domain is what the reported
        // reason names: this is the one substitution the two tests above do not
        // predict, so it is the one a caller has no way to anticipate.
        append_brake_phases(v0, stop_info);
        auto const q_after = q0 + stop_dist;
        plan_rest_to_rest(q_after);
        return {
            .reason = online_planner_substitution_reason::carry_velocity_shape_unavailable,
            .brake_duration = stop_info.T_a,
            .replan_start_position = q_after,
        };
    }

    /// @brief Append a nonzero-initial-velocity double-S from (q0, v0, 0) to
    /// (target_, 0, 0).
    ///
    /// Reuses the general Sec. 3.4.1 double-S formulation so the current velocity
    /// is carried through the profile. The resulting 7 constant-jerk phases share
    /// the sign structure of the rest-to-rest profile; only the durations differ.
    ///
    /// Returns whether the shape exists for this command. Nothing is appended
    /// when it does not, so the caller is free to plan the move a different way.
    ///
    /// @cite biagiotti2009 -- Sec. 3.4.1, eq. (3.19)-(3.27), p.79-85
    auto append_incorporate_velocity(Scalar q0, Scalar v0) -> bool
    {
        auto const created = double_s_trajectory<Scalar>::create({
            .q0 = q0,
            .q1 = target_,
            .v_max = v_max_,
            .a_max = a_max_,
            .j_max = j_max_,
            .v0 = v0,
            .v1 = Scalar{0},
        });
        if (!created.has_value()) {
            return false;
        }
        auto const& profile = *created;

        auto const durations = profile.phase_durations();
        auto const sigma = (target_ - q0 > Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const j_pos = sigma * j_max_;
        auto const j_neg = -sigma * j_max_;

        // 7-phase double-S: accel(3) + cruise(1) + decel(3).
        append_phase(durations[0], j_pos);      // jerk(+): build acceleration
        append_phase(durations[1], Scalar{0});  // constant acceleration
        append_phase(durations[2], j_neg);      // jerk(-): null acceleration at v_lim
        append_phase(durations[3], Scalar{0});  // cruise at v_lim
        append_phase(durations[4], j_neg);      // jerk(-): build deceleration
        append_phase(durations[5], Scalar{0});  // constant deceleration
        append_phase(durations[6], j_pos);      // jerk(+): null acceleration at rest
        return true;
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
            if (T_const > Scalar{0}) {
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

    /// @brief Append a single constant-jerk phase if its duration is non-negligible.
    void append_phase(Scalar duration, Scalar jerk)
    {
        if (duration > Scalar{0}) {
            T_ph_[n_phases_] = duration;
            j_ph_[n_phases_] = jerk;
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
        auto const h_signed = target_ - q0;

        if (h_signed == Scalar{0}) {
            return;
        }

        auto const sigma = (h_signed > Scalar{0}) ? Scalar{1} : Scalar{-1};
        auto const h = std::abs(h_signed);

        Scalar T_j1{}, T_a{}, T_v{}, T_d{}, T_j2{};
        compute_double_s_durations(h, T_j1, T_a, T_v, T_d, T_j2);

        auto const j_pos = sigma * j_max_;
        auto const j_neg = -sigma * j_max_;

        // 7-phase double-S: accel(3) + cruise(1) + decel(3)
        append_phase(T_j1, j_pos);
        append_phase(T_a - Scalar{2} * T_j1, Scalar{0});
        append_phase(T_j1, j_neg);
        append_phase(T_v, Scalar{0});
        append_phase(T_j2, j_neg);
        append_phase(T_d - Scalar{2} * T_j2, Scalar{0});
        append_phase(T_j2, j_pos);
    }

    /// @brief Assign doubly-degenerate durations (neither v_max nor a_max reached).
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.90
    static void assign_doubly_degenerate(Scalar h, Scalar j,
                                         Scalar& T_j1, Scalar& T_a,
                                         Scalar& T_d, Scalar& T_j2)
    {
        auto const T_j_dd = std::cbrt(h / (Scalar{2} * j));
        T_j1 = T_j_dd;
        T_j2 = T_j_dd;
        T_a = Scalar{2} * T_j_dd;
        T_d = Scalar{2} * T_j_dd;
    }

    /// @brief Solve no-cruise case where v_max is not reached.
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.89-91
    void solve_no_cruise(Scalar h, bool a_max_reached,
                         Scalar& T_j1, Scalar& T_a, Scalar& T_v,
                         Scalar& T_d, Scalar& T_j2) const
    {
        T_v = Scalar{0};
        auto const a = a_max_;
        auto const j = j_max_;

        if (a_max_reached) {
            auto const coeff_a = Scalar{1} / a;
            auto const coeff_b = a / j;
            auto const coeff_c = -h;
            auto const disc = coeff_b * coeff_b - Scalar{4} * coeff_a * coeff_c;
            auto const v_lim = (-coeff_b + std::sqrt(disc)) / (Scalar{2} * coeff_a);

            T_j1 = a / j;
            T_j2 = T_j1;
            T_a = T_j1 + v_lim / a;
            T_d = T_a;

            if (T_a < Scalar{2} * T_j1) {
                assign_doubly_degenerate(h, j, T_j1, T_a, T_d, T_j2);
            }
        } else {
            assign_doubly_degenerate(h, j, T_j1, T_a, T_d, T_j2);
        }
    }

    /// @brief Compute double-S phase durations for rest-to-rest displacement h.
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.88-91
    void compute_double_s_durations(Scalar h,
                                    Scalar& T_j1, Scalar& T_a, Scalar& T_v,
                                    Scalar& T_d, Scalar& T_j2) const
    {
        auto const v = v_max_;
        auto const a = a_max_;
        auto const j = j_max_;

        bool const a_max_reached = (v * j >= a * a);

        Scalar T_j_val = a_max_reached ? (a / j) : std::sqrt(v / j);
        Scalar T_a_val = a_max_reached ? (T_j_val + v / a) : (Scalar{2} * T_j_val);
        Scalar T_v_val = h / v - T_a_val;

        if (T_v_val > Scalar{0}) {
            T_j1 = T_j_val;
            T_j2 = T_j_val;
            T_a = T_a_val;
            T_d = T_a_val;
            T_v = T_v_val;
        } else {
            solve_no_cruise(h, a_max_reached, T_j1, T_a, T_v, T_d, T_j2);
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

}

#endif
