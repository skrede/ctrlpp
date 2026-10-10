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

#include "ctrlpp/trajectory/detail/jerk_limited_ramp.h"

#include "ctrlpp/util/concepts.h"

#include <array>
#include <cmath>
#include <limits>
#include <cstddef>
#include <cstdint>
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
    /// first, then replanned. "Too large to stop" is the carry-velocity shape's
    /// own displacement test, and the carried speed is exact by construction, so
    /// the shape refuses a command the planner admitted only through an
    /// arithmetic rejection. Such a state is braked to rest and replanned as well.
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
        , phase_end_velocity_{}
        , phase_end_acceleration_{}
        , phase_end_pinned_{}
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
    // A phase whose end state is known analytically records it, and evaluation
    // takes that state in place of the integrated one.
    static constexpr int max_phases_ = 11;
    int n_phases_{0};
    Scalar T_ph_[max_phases_]{};
    Scalar j_ph_[max_phases_]{};
    std::array<Scalar, max_phases_> phase_end_velocity_;
    std::array<Scalar, max_phases_> phase_end_acceleration_;
    std::array<bool, max_phases_> phase_end_pinned_;
    Scalar T_{};  ///< total duration

    // The roundings this planner performs on a speed it hands to the zero-acceleration planning step, taken over the
    // branches that reach it and maximized. On the branch that plans straight from the caller's snapshot it performed
    // none, and the snapshot's own error is the caller's, out of scope exactly as the trapezoidal retiming floor scopes
    // its own. On the branch that nulled a starting acceleration first, seven operations formed the speed: the phase
    // duration's division, the product of the starting acceleration with it, the three products of the jerk term, and
    // the two sums. The scaling by one half inside that term is exact in a binary radix and is counted anyway, which
    // makes the total an upper bound rather than an estimate.
    static constexpr std::int32_t speed_rounding_ops = 7;

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
        phase_end_pinned_.fill(false);

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

        // The evaluation starts from the snapshot's acceleration, so every nonzero one is nulled by its own phase,
        // however short. A phase too short to represent still pins its end state.
        if (a_start != Scalar{0}) {
            auto const T_az = std::min(std::abs(a_start) / j_max_, jerk_time_within_limit());
            auto const j_az = (a_start > Scalar{0}) ? -j_max_ : j_max_;

            // State after this phase:
            auto const v_after = snap_nulling_end_speed(v_start + a_start * T_az + Scalar{0.5} * j_az * T_az * T_az);
            auto const q_after = q_start + v_start * T_az
                                 + Scalar{0.5} * a_start * T_az * T_az
                                 + j_az * T_az * T_az * T_az / Scalar{6};

            auto const first = n_phases_;
            T_ph_[n_phases_] = T_az;
            j_ph_[n_phases_] = j_az;
            ++n_phases_;
            pin_group_end(first, v_after, Scalar{0});

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

        // The two floors below ask questions about quantities in two different
        // units, and one number cannot be the right answer to both. Each gets
        // the scale of the quantity it tests and a count of the operations that
        // formed it.

        // Speed floor, used only by the direction test. The scale is the velocity
        // limit: it is the only velocity the planner is given, and every
        // velocity the profile carries is bounded by it. The count is the
        // roundings this planner performed on the speed being tested. Below the
        // floor the speed's sign does not say which way the axis is moving.
        auto const speed_floor = static_cast<Scalar>(speed_rounding_ops)
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

        // The evaluation integrates from any nonzero speed, so only an exactly zero one is planned as rest.
        if (v0 == Scalar{0}) {
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

        // The carry shape's own displacement test on the operands create() sees, so the verdict and the
        // shape's refusal agree bit for bit.
        bool const overshoot = !wrong_way
                               && std::abs(h_signed) < detail::no_cruise_displacement(
                                      Scalar{0}, std::abs(v0), Scalar{0}, a_max_, j_max_);

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

        // The displacement test above is the shape's own and the carried speed is exact by construction, so
        // the shape refuses here only through an arithmetic rejection inside the construction. Where the
        // shape does not exist the planner brakes to rest and replans from the stopping point, the same
        // always-admissible fallback the two tests above select, reaching the same target under the same
        // limits. It costs time, not correctness, and nothing is emitted from a profile that was not built.
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
        append_double_s_phases(durations, sigma, sigma * profile.peak_velocity());
        return true;
    }

    /// @brief Append the seven phases of a double-S, accel(3) + cruise(1) + decel(3), pinning the ramp
    /// and cruise ends at (peak, 0) and the deceleration end at rest.
    void append_double_s_phases(std::array<Scalar, 7> const& durations, Scalar sigma, Scalar peak)
    {
        auto const j_pos = sigma * j_max_;
        auto const j_neg = -sigma * j_max_;

        auto const ramp = n_phases_;
        append_phase(durations[0], j_pos);
        append_phase(durations[1], Scalar{0});
        append_phase(durations[2], j_neg);
        pin_group_end(ramp, peak, Scalar{0});
        auto const cruise = n_phases_;
        append_phase(durations[3], Scalar{0});
        pin_group_end(cruise, peak, Scalar{0});
        auto const decel = n_phases_;
        append_phase(durations[4], j_neg);
        append_phase(durations[5], Scalar{0});
        append_phase(durations[6], j_pos);
        pin_group_end(decel, Scalar{0}, Scalar{0});
    }

    /// @brief Record (v, a) as the end state of the last phase appended since index `first`, if any.
    void pin_group_end(std::int32_t first, Scalar v, Scalar a)
    {
        if (n_phases_ > first) {
            auto const last = static_cast<std::size_t>(n_phases_ - 1);
            phase_end_velocity_[last] = v;
            phase_end_acceleration_[last] = a;
            phase_end_pinned_[last] = true;
        }
    }

    // Analytically the nulling-end speed of a state this planner produced is at most the ramp's peak, so an
    // excess within the counted roundings is that rounding and is placed on the limit; a larger excess is
    // left alone for the brake path.
    auto snap_nulling_end_speed(Scalar v_after) const -> Scalar
    {
        auto const excess = std::abs(v_after) - v_max_;
        auto const rounding = static_cast<Scalar>(speed_rounding_ops) * std::numeric_limits<Scalar>::epsilon() * v_max_;
        if (excess > Scalar{0} && excess <= rounding) {
            return std::copysign(v_max_, v_after);
        }
        return v_after;
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
            T_j = jerk_time_within_limit();
            auto const T_const = abs_v / a_max_ - T_j;
            T_total = T_j + T_const + T_j;
            dist = sign_v * compute_decel_distance(abs_v, T_j, T_const);
            a_reached = true;
        } else {
            T_j = std::min(std::sqrt(abs_v / j_max_), jerk_time_within_limit());
            T_total = Scalar{2} * T_j;
            dist = sign_v * abs_v * T_j;
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
        auto const first = n_phases_;
        append_phase(info.T_j, -sign_v * j_max_);
        if (info.a_max_reached) {
            append_phase(info.T_a - Scalar{2} * info.T_j, Scalar{0});
        }
        append_phase(info.T_j, sign_v * j_max_);
        pin_group_end(first, Scalar{0}, Scalar{0});
    }

    /// @brief Append a single constant-jerk phase if its duration is non-negligible.
    ///
    /// Every emitted jerk phase starts, or pinned ends, at an exactly zero acceleration, so bounding its duration bounds every
    /// acceleration it and the constant phase after it produce. Analytically none exceeds a_max / j_max: this trims only rounding.
    void append_phase(Scalar duration, Scalar jerk)
    {
        if (duration > Scalar{0}) {
            T_ph_[n_phases_] = (jerk == Scalar{0}) ? duration : std::min(duration, jerk_time_within_limit());
            j_ph_[n_phases_] = jerk;
            ++n_phases_;
        }
    }

    // The longest jerk time whose product with j_max does not round above a_max. One step suffices for a finite
    // quotient, because round-to-nearest leaves the quotient's predecessor below a_max / j_max.
    auto jerk_time_within_limit() const -> Scalar
    {
        auto const T_j = a_max_ / j_max_;
        if (j_max_ * T_j > a_max_) {
            return std::nextafter(T_j, Scalar{0});
        }
        return T_j;
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

        Scalar T_j1{}, T_a{}, T_v{}, T_d{}, T_j2{}, v_lim{};
        compute_double_s_durations(h, T_j1, T_a, T_v, T_d, T_j2, v_lim);

        append_double_s_phases({T_j1, T_a - Scalar{2} * T_j1, T_j1, T_v, T_j2, T_d - Scalar{2} * T_j2, T_j2},
                               sigma, sigma * v_lim);
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
    ///
    /// The peak `v_lim` is capped at v_max: the branch is entered only when the cruise duration is not positive,
    /// which is equivalent to a peak at most v_max, so a computed excess is rounding.
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.89-91
    void solve_no_cruise(Scalar h, bool a_max_reached, Scalar& T_j1, Scalar& T_a, Scalar& T_v,
                         Scalar& T_d, Scalar& T_j2, Scalar& v_lim) const
    {
        T_v = Scalar{0};
        auto const a = a_max_;
        auto const j = j_max_;
        if (a_max_reached) {
            auto const coeff_a = Scalar{1} / a;
            auto const coeff_b = a / j;
            auto const disc = coeff_b * coeff_b + Scalar{4} * coeff_a * h;
            auto const root = (-coeff_b + std::sqrt(disc)) / (Scalar{2} * coeff_a);

            T_j1 = a / j;
            T_j2 = T_j1;
            T_a = T_j1 + root / a;
            T_d = T_a;
            v_lim = std::min(root, v_max_);

            if (!(T_a < Scalar{2} * T_j1)) {
                return;
            }
        }
        assign_doubly_degenerate(h, j, T_j1, T_a, T_d, T_j2);
        v_lim = std::min(j * T_j1 * T_j1, v_max_);
    }

    /// @brief Compute double-S phase durations for rest-to-rest displacement h.
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.88-91
    void compute_double_s_durations(Scalar h,
                                    Scalar& T_j1, Scalar& T_a, Scalar& T_v,
                                    Scalar& T_d, Scalar& T_j2, Scalar& v_lim) const
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
            v_lim = v_max_;
        } else {
            solve_no_cruise(h, a_max_reached, T_j1, T_a, T_v, T_d, T_j2, v_lim);
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
            auto const t_remaining = dt - t_elapsed;

            if (t_remaining <= Scalar{0}) {
                break;
            }

            if (t_remaining < T_i) {
                advance_phase(i, t_remaining, q, v, a);
                break;
            }

            advance_phase(i, T_i, q, v, a);
            t_elapsed += T_i;
        }
    }

    // Inside a pinned phase the sample is formed backward from the pinned end state, the backward-time form
    // double_s_trajectory::eval_accel_phase uses for its last segment, because a forward sum can round past the
    // end value just before the boundary.
    void advance_phase(std::int32_t i, Scalar tau, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const j_i = j_ph_[i];
        q += v * tau + Scalar{0.5} * a * tau * tau + j_i * tau * tau * tau / Scalar{6};
        auto const index = static_cast<std::size_t>(i);
        if (!phase_end_pinned_[index]) {
            v += a * tau + Scalar{0.5} * j_i * tau * tau;
            a += j_i * tau;
            return;
        }
        auto const s = T_ph_[i] - tau;
        v = phase_end_velocity_[index] - phase_end_acceleration_[index] * s + Scalar{0.5} * j_i * s * s;
        a = phase_end_acceleration_[index] - j_i * s;
    }
};

}

#endif
