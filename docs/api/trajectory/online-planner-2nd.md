# online_planner_2nd

2nd-order online trajectory planner that generates trapezoidal-like velocity profiles in real time. On each `update(target)`, the planner computes a time-optimal profile from the current state to the target position, respecting `v_max` and `a_max` constraints. Suitable for real-time control loops where targets change dynamically.

Unlike pre-computed trajectory segments, online planners are stateful filters with no fixed duration and do not satisfy `trajectory_segment`.

## Header

```cpp
#include "ctrlpp/trajectory/online_planner_2nd.h"
```

## Template Parameters

| Parameter | Constraint | Description |
|-----------|------------|-------------|
| `Scalar` | arithmetic type | Numeric type (`float`, `double`, `long double`) |

## Config

```cpp
struct config {
    Scalar v_max;   // Maximum velocity magnitude
    Scalar a_max;   // Maximum acceleration magnitude

    // Settle policy. Omitting these reproduces the behavior the planner had
    // before they existed, for every choice of limits.
    Scalar position_settle_tol = 1e-9;  // position error counted as arrived
    Scalar velocity_settle_tol = 1e-9;  // speed counted as stopped
};
```

Both limits divide in the planner math (stopping distance `v^2 / (2 * a_max)`, phase durations `v_v / a_max` and `h / v_v`), so the domain of each is finite and strictly positive.

### Tuning the settle policy

The two tolerances decide when [`is_settled`](#is_settled) starts reporting true. Each is compared against exactly one sampled quantity: `position_settle_tol` against the distance still to run, `velocity_settle_tol` against the speed. They are two fields rather than one because they are compared against two quantities in two different units, and a single number compared against both a length and a speed is making two unrelated claims at once.

```cpp
// An axis that counts as arrived within a tenth of a millimeter, whose encoder
// resolves nothing finer than a millimeter per second.
auto result = ctrlpp::online_planner_2nd<double>::create({
    .v_max = 1.0,
    .a_max = 5.0,
    .position_settle_tol = 1e-4,
    .velocity_settle_tol = 1e-3,
});
```

There is no acceleration tolerance here. This planner bounds no jerk and carries no acceleration state to settle, so a third field would be surface with nothing behind it. The [3rd-order planner](online-planner-3rd.md), which does carry one, has three.

They are knobs rather than derived constants because no derivation exists for them. The distance at which an axis counts as arrived is a property of the machine, not of the arithmetic. The default's provenance is stated rather than implied: `1e-9` is the value the planner compared both residuals against before the fields existed, and it sits roughly seven decades above `double` rounding on the quantities it tests. A number that far above the resolution of its own operands is an application statement, not a rounding guard, so it belongs to the caller. Omitting the fields reproduces that earlier behavior exactly, for every choice of limits.

**The defaults are chosen for `double`. A `float` user should set the fields.** On `float`, `1e-9` sits about two decades below the type's own resolution near unity, so the default is a condition a sampled value near unity effectively never satisfies unless it is exactly zero. The defaults are deliberately not type-dependent: a type-dependent default would stop omission reproducing the earlier behavior.

**The knob cannot break the result contract.** No setting of it changes what the planner computes. The profile, its phase durations and every value [`sample`](#sample) returns are identical under any value; only the moment `is_settled` flips is moved.

#### The window this leaves open

Two different distance questions live in the planner, and they have different answers on purpose.

1. **"Is the motion done?"** is asked by `sample`, at the policy distance above. It is the caller's to choose.
2. **"Is this command a numerical no-op?"** is asked when a new target arrives. It is not policy and is not tunable. It holds only when the commanded displacement is exactly zero and the current velocity is exactly zero. A target differing from the current position at all is planned in full however wide the settle policy is.

Take the axis `v_max = 1`, `a_max = 5` and `double` arithmetic, where `eps = 2.22e-16`.

| Question | Compared against | Value on that axis | Where it comes from |
|---|---|---|---|
| Is the motion done? | `velocity_settle_tol` | `1e-9` | Policy. The default is the value the planner used before the field existed. No derivation, because none exists: the speed at which a machine counts as stopped is a property of the machine. |
| Is this command a numerical no-op? | exactly zero | `0` | Derived. Removing any nonzero speed takes `abs(v) / a_max` at the acceleration limit, however small the speed is. The planner is never told the sample period, so no magnitude exists below which dropping it would respect `a_max`. |

**The gap is deliberate.** A policy tolerance must not decide what the planner is allowed to compute. The settle policy states when the application considers the axis arrived, which is a statement about the machine; the no-op test states when there is nothing left to remove, which is a statement about the limits. Collapsing them would let a caller who widens the arrival distance silently stop the planner from computing moves it can represent perfectly well.

The window between the two is the range in which the planner reports the axis arrived while remaining willing to plan a move to a nearer target. Widening the settle policy widens that window. It never narrows what the planner will plan.

A command at the current position that still carries a speed is not `settled`. The planner brings the speed to rest at `a_max` and ends at the target, which takes a short time. A speed above the speed floor below reports `braked_and_replanned` with `reversal_or_overshoot`. A speed at or below it, at most `eps * v_max`, is carried rather than braked first and reports `commanded_profile`, with the velocity continuous at the command and within `v_max` and `a_max`. A negative speed is reversed and returned to the target. A positive speed is brought to rest with its stopping point placed at the target, up to `(eps * v_max)^2 / (2 a_max)` short of where braking from that speed ends. That distance is below the target's own resolution unless the target lies within `eps * v_max^2 / a_max` of zero.

#### Constants the planner derives rather than fixes

Nothing in this group is tunable, and none of it is an absolute number. Each floor is a counted number of rounding operations at the scale of the quantity the comparison is about, with the count enumerated in the header beside the code it guards. The overshoot verdict carries no count at all: it is the braking distance the planner also places its stop with.

| Decision | Bound | Scale, and why it is that scale |
|---|---|---|
| Is the commanded speed zero? | `1 * eps * v_max` | The velocity limit. The count is one because this planner performs no arithmetic on the speed before testing it: it bounds no jerk, so it has no acceleration-nulling phase to run first. |
| Is the commanded displacement zero? | `4 * eps * v_max^2 / (2 a_max)` | The planner's own stopping distance from full speed. It is intrinsic to the limits and needs no knowledge of the sample period, which the planner is never told. One operation forms the displacement and three form the scale. |
| Would the move overshoot? | `\|v0\| * (\|v0\| / a_max) / 2 > \|h\|` | Neither. The left side is the braking distance, the mean speed over the braking duration, and it is the same value the planner adds to the start position to place its stopping point. Both sides are lengths, so the comparison needs no external scale and carries no slack. |

**These counts are not the jerk-limited planner's counts, and copying them across would be wrong.** A count that is not counted against the code it guards is an unexplained constant wearing a different name.

The overshoot verdict here is scale-invariant for the same reason as the jerk-limited planner's: the same relative shortfall is reported as an overshoot on an axis whose stopping distance is meters and on one whose stopping distance is picometers.

#### One braking distance for the verdict and the stopping point

The overshoot test compares the braking distance against the distance that remains, and a braked command places its stopping point at the start position plus that same braking distance, formed once as the mean speed times the braking duration. The verdict and the placement therefore cannot disagree. The mean-speed form is the one used because it forms no square of the speed.

**Commands just short of the braking distance are reported as overshoots.** A commanded target whose remaining distance falls below the braking distance by as little as one unit in the last place reports `reversal_or_overshoot`: the axis brakes to rest one to four units in the last place past the target and returns to it. The motion respects every limit either way; the report is what tells the caller the move took longer than commanded.

The agreement is measured as well, in `double` and in `float`: 13,499 commanded targets in `double` and 9,634 in `float`, over 60 limit-set, start-position and initial-speed configurations, find zero disagreements between the planner's verdict and the braking distance formula above. Each boundary is located by bisection in the space of representable targets, and every representable target within 64 units in the last place of it is walked one at a time. On all 60, in both types, the two boundaries land on the same representable target. In `float` a unit in the last place of a position near `1e4` is about `1e-3`, so on some configurations the walk reaches behind the start position. Those 415 targets ask the direction test's question rather than the overshoot test's, so they are left out of the 9,634, and each one was checked against the direction test's answer: every one reports `reversal_or_overshoot`.

**What the sweep covers, and what it does not.** The reference here is the planner's own documented braking distance rather than a shape constructor's rejection, because [`double_s_trajectory`](double-s-trajectory.md) is not reachable from this planner at all, so the sweep checks the implementation against the formula. States on the acceleration ramp at four fractions of `v_max`; this planner carries no acceleration state, so a ramp state is any state a caller can present. Forward commands only, since a reversing command is decided by the direction test rather than the overshoot test. `double` and `float`; `long double` is not measured. The report is byte-identical across eight builds, all on x86-64 Linux against Eigen 3.4.1: g++ 16.2.1 at `-O0`, `-O2` and `-O3`, and at `-O2` with `-ffp-contract=off` and with `-ffp-contract=fast`; clang 23.1.1 at `-O2` with and without `-ffp-contract=fast`; and clang 18.1.8 at `-O2` against libstdc++ 15.3.0.

## Construction

```cpp
static auto create(config const& cfg)
    -> ctrlpp::expected<online_planner_2nd, trajectory_error>;
```

`create` is the only construction path. There is no non-fallible constructor: a rejected configuration is a value the caller has to inspect, never an object that quietly stands in for one. It validates the kinematic limits and constructs a planner with initial state at rest at q = 0. Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| NaN/Inf or non-positive `v_max` | `trajectory_error::non_positive_velocity_limit` |
| NaN/Inf or non-positive `a_max` | `trajectory_error::non_positive_acceleration_limit` |

## Methods

### update

```cpp
auto update(Scalar target) -> ctrlpp::expected<void, trajectory_error>;
```

Set a new target position and replan from the current state. Computes a time-optimal trapezoidal profile from (q, v) to (target, 0) respecting `v_max` and `a_max`. A velocity pointing away from the target, or too large to stop in the available distance, is braked to rest first and the move is replanned from the stopping point.

The command is rejected with `trajectory_error::non_finite_input` before any
planner state changes when the target is NaN or infinite. For an accepted
target, which profile was built is read back from
[`diagnostics()`](#diagnostics). The motion respects every limit either way;
what changes is the time it takes.

### diagnostics

```cpp
auto diagnostics() const -> online_planner_diagnostics<Scalar> const&;
```

Report what the last `update` (or `reset`) planned, against what it was commanded. This describes the plan, not the state reached since; `is_settled()` answers that.

### sample

```cpp
auto sample(Scalar t) const -> trajectory_point<Scalar, 1>;
```

Evaluate the planned trajectory at absolute time `t`. Returns position, velocity, and acceleration. Updates internal state for future `update()` calls. Call this once per control loop timestep.

### is_settled

```cpp
auto is_settled() const -> bool;
```

Returns true when the planner has reached the target with zero velocity.

### reset

```cpp
void reset(Scalar q0);
```

Reset state to position `q0` with zero velocity.

## Profile Behavior

The planner produces trapezoidal velocity profiles with three phases:

1. **Acceleration**<br/>ramp velocity toward `v_max` at rate `a_max`
2. **Cruise**<br/>hold at `v_max` (may be zero duration for short moves)
3. **Deceleration**<br/>ramp velocity to zero at rate `a_max`

For short displacements where `v_max` cannot be reached, the profile degenerates to a triangular velocity shape. Mid-motion target changes trigger replanning from the current state, with automatic brake-and-replan for reversal and overshoot scenarios.

`sample` uses the deceleration ramp's formula only at times strictly inside that ramp, measured backward from the end of the plan. So at the switch from cruise it reports the cruise velocity itself, never a rounded value above `v_max`. A carried acceleration ramp, which a retarget starts from a nonzero velocity, is sampled from whichever of its two exact end states is nearer in time, so it cannot round past `v_max` either.

This is measured in `double` and in `float` by sampling every representable instant within 64 steps of each phase switch. The planner is retargeted further along its direction from a seeded instant in the first half of its acceleration ramp, so the new plan starts from a carried speed and a carried ramp, in both directions and over limit pairs that include one with full-width significands: 9,216 retargets per type. Over 1,188,864 sampled instants per type, zero samples exceed `v_max` or `a_max`, and every window shows the switch it was placed on. The output is byte-identical from the eight builds listed above.

## Substitution Reporting

The brake-and-replan is a substitution: the caller commanded a move from the current velocity, and the planner realized a different, longer profile that respects the same limits and reaches the same target. The motion alone does not distinguish the two. The disposition does.

```cpp
enum class online_planner_disposition
{
    commanded_profile,     // the commanded shape was planned as asked
    braked_and_replanned,  // braked to rest, then replanned from the stopping point
    settled,               // exactly at the target with zero velocity; a zero-duration profile
};

enum class online_planner_substitution_reason
{
    none,
    reversal_or_overshoot,             // the commanded motion reverses, or would overshoot
    carry_velocity_shape_unavailable,  // reported only by the 3rd-order planner
};

template <typename Scalar>
struct [[nodiscard]] online_planner_diagnostics
{
    online_planner_disposition          disposition;
    online_planner_substitution_reason  substitution_reason;
    Scalar commanded_target;
    Scalar initial_velocity;
    Scalar planned_duration;
    Scalar brake_duration;
    Scalar replan_start_position;
};
```

| Field | Meaning |
|-------|---------|
| `disposition` | Which profile was built |
| `substitution_reason` | Which condition selected the substitution; `none` when nothing was substituted |
| `commanded_target` | Target position the update was given |
| `initial_velocity` | Velocity the planner was carrying when it planned |
| `planned_duration` | Total duration the plan realizes |
| `brake_duration` | Duration spent braking before the replan; zero when nothing was substituted |
| `replan_start_position` | Position the replan starts from; the commanded start position when nothing was substituted |

The type is shared with [online-planner-3rd](online-planner-3rd.md). This planner bounds no jerk, so it has no carry-velocity shape whose domain a commanded state can fall outside of: it reports `none` or `reversal_or_overshoot` and nothing else.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include <ctrlpp/trajectory/online_planner_2nd.h>

#include <iostream>

int main()
{
    auto result = ctrlpp::online_planner_2nd<double>::create({.v_max = 1.0, .a_max = 5.0});
    if (!result.has_value())
        return 1;
    auto& planner = *result;
    if (!planner.update(10.0).has_value())
        return 1;

    constexpr double dt = 0.001;
    for (double t = 0.0; !planner.is_settled(); t += dt) {
        auto pt = planner.sample(t);
        std::cout << t << "," << pt.position(0) << "," << pt.velocity(0) << "\n";
    }
}
```

## See Also

- [online-planner-3rd](online-planner-3rd.md)<br/> 3rd-order variant with jerk limiting
- [trapezoidal-trajectory](trapezoidal-trajectory.md)<br/> Pre-computed trapezoidal profile
- [Trajectory Generation Theory](../../background/trajectory-generation.md)<br/> Online trajectory generation algorithms
