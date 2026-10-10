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
        , n_phases_{0}
        , target_anchor_first_{0}
        , T_ph_{}
        , j_ph_{}
        , phase_end_position_{}
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

    // At most nulling 1 + brake 3 + ramp 3 + cruise 1 + deceleration 3; an emptied group stores one pinned phase.
    // A phase whose end state is known analytically records it, and evaluation takes that state in place of the integrated one.
    // From the cruise end on, or over the last phase when there is no double-S, end positions are formed backward from the target.
    static constexpr std::int32_t max_phases_ = 11;
    std::int32_t n_phases_;
    std::int32_t target_anchor_first_;
    std::array<Scalar, max_phases_> T_ph_;
    std::array<Scalar, max_phases_> j_ph_;
    std::array<Scalar, max_phases_> phase_end_position_;
    std::array<Scalar, max_phases_> phase_end_velocity_;
    std::array<Scalar, max_phases_> phase_end_acceleration_;
    std::array<bool, max_phases_> phase_end_pinned_;
    Scalar T_{};  ///< total duration

    // Upper bound on the roundings in the nulling-end speed: the duration's division, the product a T, the three products of
    // the jerk term (its exact halving counted anyway) and the two sums. A plan straight from the snapshot performs none.
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
        target_anchor_first_ = max_phases_;
        phase_end_pinned_.fill(false);

        // A command at the target is a no-op only from rest: removing any nonzero velocity or acceleration takes time
        // at the limits, and no magnitude exists below which dropping it would respect them.
        if (target_ == q_ref_ && v_ref_ == Scalar{0} && a_ref_ == Scalar{0}) {
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
        if (a_ref_ != Scalar{0}) {
            append_nulling_phase(q_start, v_start);
        }

        // Now state is (q_start, v_start, 0). Plan double-S from here.
        auto const outcome = plan_from_zero_accel(q_start, v_start);
        anchor_target_positions();

        T_ = Scalar{0};
        for (std::int32_t i = 0; i < n_phases_; ++i) {
            T_ += T_ph_[static_cast<std::size_t>(i)];
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

    // The evaluation starts from the snapshot's acceleration, so every nonzero one is nulled by its own phase,
    // however short. A phase too short to represent still pins its end state.
    void append_nulling_phase(Scalar& q_start, Scalar& v_start)
    {
        auto const a_start = a_ref_;
        auto const T_az = std::min(std::abs(a_start) / j_max_, detail::jerk_time_within(a_max_, j_max_));
        auto const j_az = (a_start > Scalar{0}) ? -j_max_ : j_max_;
        auto const v_after = snap_nulling_end_speed(v_start + a_start * T_az + Scalar{0.5} * j_az * T_az * T_az);
        q_start = q_start + v_start * T_az + Scalar{0.5} * a_start * T_az * T_az + j_az * T_az * T_az * T_az / Scalar{6};
        auto const first = n_phases_;
        append_phase(T_az, j_az);
        pin_group_end(first, v_after, Scalar{0});
        v_start = v_after;
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
        constexpr std::int32_t nulling_position_rounding_ops = 12;
        constexpr std::int32_t displacement_rounding_ops = 1;
        constexpr std::int32_t length_scale_rounding_ops = 3;
        constexpr std::int32_t length_rounding_ops = nulling_position_rounding_ops
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
        append_hold_phase(durations[1], durations[2], peak, j_neg);
        append_phase(durations[2], j_neg);
        pin_group_end(ramp, peak, Scalar{0});
        auto const cruise = n_phases_;
        append_phase(durations[3], Scalar{0});
        pin_group_end(cruise, peak, Scalar{0});
        target_anchor_first_ = n_phases_ - 1;
        auto const decel = n_phases_;
        append_phase(durations[4], j_neg);
        append_hold_phase(durations[5], durations[6], Scalar{0}, j_pos);
        append_phase(durations[6], j_pos);
        pin_group_end(decel, Scalar{0}, Scalar{0});
    }

    // Pinned where the closing jerk phase's backward form from (v_end, 0) starts, so a sample in its second half
    // never passes v_end however short that phase is.
    void append_hold_phase(Scalar duration, Scalar closing_duration, Scalar v_end, Scalar j_close)
    {
        auto const hold = n_phases_;
        append_phase(duration, Scalar{0});
        if (duration > Scalar{0}) {
            auto const T_c = std::min(closing_duration, detail::jerk_time_within(a_max_, j_max_));
            pin_group_end(hold, v_end + Scalar{0.5} * j_close * T_c * T_c, -j_close * T_c);
        }
    }

    // Walks back from (target_, 0, 0) through the pinned or back-integrated end states, so on the final approach no
    // forward sum of earlier roundings can carry a position past the target.
    void anchor_target_positions()
    {
        target_anchor_first_ = std::max(std::min(target_anchor_first_, n_phases_ - 1), std::int32_t{0});
        auto q = target_;
        Scalar v{};
        Scalar a{};
        for (auto i = n_phases_ - 1; i >= target_anchor_first_; --i) {
            auto const index = static_cast<std::size_t>(i);
            if (phase_end_pinned_[index]) {
                v = phase_end_velocity_[index];
                a = phase_end_acceleration_[index];
            }
            phase_end_position_[index] = q;
            auto const T = T_ph_[index];
            auto const j = j_ph_[index];
            q -= travel(v, -a, j, T);
            v = v - a * T + Scalar{0.5} * j * T * T;
            a = a - j * T;
        }
    }

    static auto travel(Scalar v, Scalar a, Scalar j, Scalar t) -> Scalar
    {
        return v * t + Scalar{0.5} * a * t * t + j * t * t * t / Scalar{6};
    }

    /// @brief Record (v, a) as the end state of the last phase appended since index `first`, storing a zero-length
    /// phase to carry it when the group stored none.
    void pin_group_end(std::int32_t first, Scalar v, Scalar a)
    {
        if (n_phases_ == first) {
            T_ph_[static_cast<std::size_t>(n_phases_)] = Scalar{0};
            j_ph_[static_cast<std::size_t>(n_phases_)] = Scalar{0};
            ++n_phases_;
        }
        auto const last = static_cast<std::size_t>(n_phases_ - 1);
        phase_end_velocity_[last] = v;
        phase_end_acceleration_[last] = a;
        phase_end_pinned_[last] = true;
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
            T_j = detail::jerk_time_within(a_max_, j_max_);
            auto const T_const = abs_v / a_max_ - T_j;
            T_total = T_j + T_const + T_j;
            dist = sign_v * compute_decel_distance(abs_v, T_j, T_const);
            a_reached = true;
        } else {
            T_j = std::min(std::sqrt(abs_v / j_max_), detail::jerk_time_within(a_max_, j_max_));
            T_total = Scalar{2} * T_j;
            dist = sign_v * abs_v * T_j;
            a_reached = false;
        }

        return {dist, T_j, T_total, a_reached};
    }

    /// @brief Compute distance during jerk-limited deceleration (a_max reached case).
    auto compute_decel_distance(Scalar abs_v, Scalar T_j, Scalar T_const) const -> Scalar
    {
        auto const v1 = abs_v - j_max_ * T_j * T_j / Scalar{2};
        auto const q1 = abs_v * T_j - j_max_ * T_j * T_j * T_j / Scalar{6};
        auto const v2 = v1 - a_max_ * T_const;
        auto const q2 = q1 + v1 * T_const - a_max_ * T_const * T_const / Scalar{2};
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

    /// @brief Append a constant-jerk phase if its duration is positive. A jerk phase starts or pinned ends at zero
    /// acceleration, so bounding its duration, which analytically never exceeds a_max / j_max, bounds every acceleration.
    void append_phase(Scalar duration, Scalar jerk)
    {
        if (duration > Scalar{0}) {
            auto const index = static_cast<std::size_t>(n_phases_);
            T_ph_[index] = (jerk == Scalar{0}) ? duration : std::min(duration, detail::jerk_time_within(a_max_, j_max_));
            j_ph_[index] = jerk;
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

        Scalar T_j1{}, T_a{}, T_v{}, T_d{}, T_j2{}, v_lim{};
        compute_double_s_durations(h, T_j1, T_a, T_v, T_d, T_j2, v_lim);

        append_double_s_phases({T_j1, T_a - Scalar{2} * T_j1, T_j1, T_v, T_j2, T_d - Scalar{2} * T_j2, T_j2},
                               sigma, sigma * v_lim);
    }

    /// @brief Assign doubly-degenerate durations (neither v_max nor a_max reached).
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.90
    void assign_doubly_degenerate(Scalar h, Scalar& T_j1, Scalar& T_a, Scalar& T_d, Scalar& T_j2, Scalar& v_lim) const
    {
        auto const T_j_dd = detail::doubly_degenerate_jerk_time(h, v_max_, a_max_, j_max_);
        T_j1 = T_j_dd;
        T_j2 = T_j_dd;
        T_a = Scalar{2} * T_j_dd;
        T_d = Scalar{2} * T_j_dd;
        v_lim = j_max_ * T_j_dd * T_j_dd;
    }

    /// @brief Solve the no-cruise case, where c = v_max T_a - h >= 0 is the distance the cruise-free profile at
    /// v_max lacks. Each peak is formed at or below v_max, by the same arithmetic as double_s_trajectory.
    /// @cite biagiotti2009 -- Sec. 3.4.3, p.89-91
    void solve_no_cruise(Scalar h, Scalar c, bool a_max_reached, Scalar& T_j1, Scalar& T_a, Scalar& T_d, Scalar& T_j2,
                         Scalar& v_lim) const
    {
        if (a_max_reached) {
            v_lim = detail::no_cruise_peak(c, h, v_max_, a_max_, j_max_);
            T_j1 = a_max_ / j_max_;
            T_j2 = T_j1;
            T_a = T_j1 + v_lim / a_max_;
            T_d = T_a;
            if (!(T_a < Scalar{2} * T_j1)) {
                return;
            }
        }
        assign_doubly_degenerate(h, T_j1, T_a, T_d, T_j2, v_lim);
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

        auto const ramp = detail::make_velocity_ramp(v, a, j);
        auto const c = v * ramp.T - h;

        if (c < Scalar{0}) {
            T_j1 = ramp.T_j;
            T_j2 = ramp.T_j;
            T_a = ramp.T;
            T_d = ramp.T;
            T_v = -c / v;
            v_lim = v_max_;
        } else {
            T_v = Scalar{0};
            solve_no_cruise(h, c, a_max_reached, T_j1, T_a, T_d, T_j2, v_lim);
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

        for (std::int32_t i = 0; i < n_phases_; ++i) {
            auto const T_i = T_ph_[static_cast<std::size_t>(i)];
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

    // Inside a pinned phase the sample is formed backward from the pinned end state, because a forward sum can round past the end value just before the boundary.
    // A pinned zero-jerk phase is sampled from its nearer end, since near its start the pin's magnitude would quantize a much smaller velocity.
    // Positions take the backward form only where they are anchored to the target; before that they are summed forward.
    void advance_phase(std::int32_t i, Scalar tau, Scalar& q, Scalar& v, Scalar& a) const
    {
        auto const index = static_cast<std::size_t>(i);
        auto const j_i = j_ph_[index];
        auto const s = T_ph_[index] - tau;
        bool const forward = !phase_end_pinned_[index] || (j_i == Scalar{0} && tau < s);
        if (forward || i < target_anchor_first_) {
            q += travel(v, a, j_i, tau);
        } else {
            q = phase_end_position_[index] - travel(phase_end_velocity_[index], -phase_end_acceleration_[index], j_i, s);
        }
        if (forward) {
            v += a * tau + Scalar{0.5} * j_i * tau * tau;
            a += j_i * tau;
            return;
        }
        v = phase_end_velocity_[index] - phase_end_acceleration_[index] * s + Scalar{0.5} * j_i * s * s;
        a = phase_end_acceleration_[index] - j_i * s;
    }
};

}

#endif
