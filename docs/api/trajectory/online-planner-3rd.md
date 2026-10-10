# online_planner_3rd

3rd-order online trajectory planner that generates double-S (jerk-limited) velocity profiles in real time. On each `update(target)`, the planner computes a time-optimal profile from the current state to the target position, respecting `v_max`, `a_max`, and `j_max` constraints. Produces smoother motion than the 2nd-order planner at the cost of slightly longer move times.

Unlike pre-computed trajectory segments, online planners are stateful filters with no fixed duration and do not satisfy `trajectory_segment`.

## Header

```cpp
#include "ctrlpp/trajectory/online_planner_3rd.h"
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
    Scalar j_max;   // Maximum jerk magnitude

    // Settle policy. Omitting these reproduces the behavior the planner had
    // before they existed, for every choice of limits.
    Scalar position_settle_tol     = 1e-9;  // position error counted as arrived
    Scalar velocity_settle_tol     = 1e-9;  // speed counted as stopped
    Scalar acceleration_settle_tol = 1e-9;  // acceleration counted as at rest
};
```

All three limits divide in the planner math (cruise duration `h / v_max`, jerk-phase durations `a_max / j_max` and `|a| / j_max`, the `a_max`-reached threshold `a_max^2 / j_max`), so the domain of each is finite and strictly positive.

### Tuning the settle policy

The three tolerances decide when [`is_settled`](#is_settled) starts reporting true. Each is compared against exactly one sampled quantity: `position_settle_tol` against the distance still to run, `velocity_settle_tol` against the speed, `acceleration_settle_tol` against the acceleration. They are three fields rather than one because they are compared against three quantities in three different units, and a single number compared against a length, a speed and an acceleration is making three unrelated claims at once.

```cpp
// An axis that counts as arrived within a tenth of a millimeter, whose encoder
// resolves nothing finer than a millimeter per second.
auto result = ctrlpp::online_planner_3rd<double>::create({
    .v_max = 1.0,
    .a_max = 5.0,
    .j_max = 50.0,
    .position_settle_tol = 1e-4,
    .velocity_settle_tol = 1e-3,
    .acceleration_settle_tol = 1e-2,
});
```

They are knobs rather than derived constants because no derivation exists for them. The distance at which an axis counts as arrived is a property of the machine, not of the arithmetic. The default's provenance is stated rather than implied: `1e-9` is the value the planner compared all three residuals against before the fields existed, and it sits roughly seven decades above `double` rounding on the quantities it tests. A number that far above the resolution of its own operands is an application statement, not a rounding guard, so it belongs to the caller. Omitting the fields reproduces that earlier behavior exactly, for every choice of limits.

**The defaults are chosen for `double`. A `float` user should set the fields.** On `float`, `1e-9` sits about two decades below the type's own resolution near unity, so the default is a condition a sampled value near unity effectively never satisfies unless it is exactly zero. The defaults are deliberately not type-dependent: a type-dependent default would stop omission reproducing the earlier behavior.

**The knob cannot break the result contract.** No setting of it changes what the planner computes. The profile, its phase durations and every value [`sample`](#sample) returns are identical under any value; only the moment `is_settled` flips is moved.

#### The window this leaves open

Two different distance questions live in the planner, and they have different answers on purpose.

1. **"Is the motion done?"** is asked by `sample`, at the policy distance above. It is the caller's to choose.
2. **"Is this command a numerical no-op?"** is asked when a new target arrives. It is not policy and is not tunable. Its position side is an exact comparison, so a target differing from the current position at all is planned in full however wide the settle policy is; its velocity and acceleration sides are answered against the resolution of the limit that bounds each, derived below.

Take the axis `v_max = 1`, `a_max = 5`, `j_max = 50` and `double` arithmetic, where `eps = 2.22e-16`.

| Question | Compared against | Value on that axis | Where it comes from |
|---|---|---|---|
| Is the motion done? | `velocity_settle_tol` | `1e-9` | Policy. The default is the value the planner used before the field existed. No derivation, because none exists: the speed at which a machine counts as stopped is a property of the machine. |
| Is this command a numerical no-op? | `1 * eps * v_max` | `2.22e-16` | Derived. One operation forms the compared quantity, the comparison itself, because the planner does no arithmetic on the caller's snapshot before testing it; the scale is the velocity limit, which is the only speed the planner is given and the bound on every speed the profile carries. |

**The gap is about 6.7 decades, and it is deliberate.** A policy tolerance must not decide what the planner is allowed to compute. The settle policy states when the application considers the axis arrived, which is a statement about the machine; the no-op floor states when the arithmetic can no longer tell one command from another, which is a statement about the type and the limits. Collapsing them would let a caller who widens the arrival distance silently stop the planner from computing moves it can represent perfectly well.

The acceleration side reads the same way: the policy is `acceleration_settle_tol`, and the no-op floor is `1 * eps * a_max`, which is `1.11e-15` on that axis.

The window between the two is the range in which the planner reports the axis arrived while remaining willing to plan a move to a nearer target. Widening the settle policy widens that window. It never narrows what the planner will plan.

#### Constants the planner derives rather than fixes

Nothing in this group is tunable, and none of it is an absolute number. Each floor is a counted number of rounding operations at the scale of the quantity the comparison is about, with the count enumerated in the header beside the code it guards. The overshoot verdict carries no count at all: it is the carry-velocity shape's own test.

| Decision | Bound | Scale, and why it is that scale |
|---|---|---|
| Does the commanded speed point away from the target? | `7 * eps * v_max` | The velocity limit. Seven operations form the speed on the branch that first nulls a starting acceleration; on the branch that does not, the speed is the caller's snapshot and its error is the caller's. Below the bound the speed's sign is not read as a direction. Only an exactly zero speed is planned as rest: any nonzero speed is carried into the plan, through the carry-velocity shape or a brake. |
| Is the commanded displacement zero? | `16 * eps * v_max^2 / (2 a_max)` | The planner's own stopping distance from full speed, which is `0.1` m on the axis above. It is intrinsic to the limits and needs no knowledge of the sample period, which the planner is never told. The count is twelve for the position the acceleration-nulling phase leaves behind, one for the subtraction that forms the displacement, and three for the scale itself. |
| Would the move overshoot? | `\|h\| < h_min(\|v0\|)` | Neither. `h_min(\|v0\|)` is the distance the fastest admissible transition from the current speed to rest sweeps, and the test is the minimum-displacement test of [`double_s_trajectory`](double-s-trajectory.md) itself, evaluated by the same function on the same operands its constructor sees. Both sides are lengths, so the comparison needs no external scale and carries no slack. |
| Is the speed after nulling the limit itself? | `7 * eps * v_max` above `v_max` | The velocity limit, with the same seven operations as the speed floor. Once its acceleration is nulled, a state this planner produced is at most its ramp's peak, so an excess within the count is that rounding and the speed is placed on `v_max`. A larger excess is left alone for the brake path. |

A nonzero starting acceleration is always brought to zero by its own jerk phase first, so the profile is planned from exactly the state it is evaluated from.

Two consequences worth knowing about. First, **the same absolute displacement is treated differently on axes with different limits**, which is the point of scaling the floor: an axis whose stopping distance is 1.25 m treats a commanded femtometer as zero, and one whose stopping distance is 5 mm plans it as a move. Second, **the overshoot verdict is scale-invariant**: the same relative shortfall is reported as an overshoot on an axis whose stopping distance is meters and on one whose stopping distance is picometers.

#### Where the overshoot verdict and the carry-velocity shape disagree

The overshoot test above is the carry-velocity shape's own displacement test, so where the planner admits a command the shape's displacement test admits it too, and where the shape refuses a displacement the planner has already reported `reversal_or_overshoot`. The two agree at that boundary by construction. The agreement is measured as well: 4,048 commanded targets over 18 limit-set and start-position configurations, `double` throughout, each boundary located by bisection in the space of representable targets and every representable target within 64 units in the last place of it walked one at a time, find zero disagreements. On all 18 the planner's boundary and the shape's land on the same representable target, and on all 18 the cruise speed the planner samples is exactly the profile's peak, `v_max` itself, so on none of them is the carry-velocity shape handed a boundary velocity above its limit.

**What the sweep covers, and what it does not.** Cruise states only, where the sampled acceleration is exactly zero and the command reaches the overshoot test carrying the state the planner itself reported; a command issued from a state with nonzero acceleration routes through the acceleration-nulling phase first and is not covered. Forward commands only, since a reversing command is decided by the direction test rather than the overshoot test. `double` only. The report was produced by one build, g++ 16.2.1 at `-O2` against Eigen 3.4.1.

## Construction

```cpp
static auto create(config const& cfg)
    -> ctrlpp::expected<online_planner_3rd, trajectory_error>;
```

`create` is the only construction path. There is no non-fallible constructor: a rejected configuration is a value the caller has to inspect, never an object that quietly stands in for one. It validates the kinematic limits and constructs a planner with initial state at rest at q = 0 with zero acceleration. Rejections, checked in order:

| Condition | Error |
|-----------|-------|
| NaN/Inf or non-positive `v_max` | `trajectory_error::non_positive_velocity_limit` |
| NaN/Inf or non-positive `a_max` | `trajectory_error::non_positive_acceleration_limit` |
| NaN/Inf or non-positive `j_max` | `trajectory_error::non_positive_jerk_limit` |

## Methods

### update

```cpp
auto update(Scalar target) -> ctrlpp::expected<void, trajectory_error>;
```

Set a new target position and replan from the current state. Computes a time-optimal double-S profile from (q, v, a) to (target, 0, 0) respecting `v_max`, `a_max`, and `j_max`. A same-direction move carries the current velocity through the profile (no full-stop dip); a velocity pointing away from the target, or too large to stop in the available distance, is braked to rest first and then replanned. "Too large to stop" is the carry-velocity shape's own displacement test, and the carried speed is exact by construction, so the shape refuses a command the planner admitted only through an arithmetic rejection. Such a state is braked to rest and replanned with `carry_velocity_shape_unavailable` reported.

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

Evaluate the planned trajectory at absolute time `t`. Returns position, velocity, and acceleration. Updates internal state for future `update()` calls.

### is_settled

```cpp
auto is_settled() const -> bool;
```

Returns true when the planner has reached the target with zero velocity and zero acceleration.

### reset

```cpp
void reset(Scalar q0);
```

Reset state to position `q0` with zero velocity and zero acceleration.

## Profile Behavior

The planner produces double-S velocity profiles composed of constant-jerk phases (up to 11 phases total):

1. **Acceleration ramp-up**<br/>apply +j_max until a_max reached
2. **Constant acceleration**<br/>hold at a_max (may be zero duration)
3. **Acceleration ramp-down**<br/>apply -j_max to bring acceleration to zero
4. **Cruise**<br/>hold at v_max with zero acceleration (may be zero duration)
5. **Deceleration ramp-up**<br/>apply -j_max to build deceleration
6. **Constant deceleration**<br/>hold at -a_max (may be zero duration)
7. **Deceleration ramp-down**<br/>apply +j_max to bring acceleration and velocity to zero

Degenerate cases (v_max or a_max not reached) automatically reduce the number of phases.

The sampled velocity and acceleration never exceed `v_max` and `a_max`, including on a replan from a ramp or from cruise. The end of each acceleration and deceleration ramp, of the cruise, of the acceleration-nulling phase and of the brake is evaluated from its analytic end state rather than from the rounded sum of the phases before it, and a sample inside such a phase is formed backward from that end state, so no sum that could round past the limit just before the boundary is ever formed. Every jerk phase the planner emits lasts at most the longest time whose product with `j_max` does not round above `a_max`, one representable step below `a_max / j_max` where that quotient rounds up. So the acceleration a jerk phase reaches, and a constant-acceleration phase holds, never exceeds `a_max`, including at the speed `a_max^2 / j_max` where a brake's two shapes meet. No sampled value is clamped.

## Substitution Reporting

A mid-motion replan keeps the current velocity where the carry-velocity shape exists, and brakes to rest and replans from the stopping point where it does not. Both outcomes respect every limit and both reach the target, so the motion alone does not tell them apart. The disposition does.

```cpp
enum class online_planner_disposition
{
    commanded_profile,     // the commanded shape was planned as asked
    braked_and_replanned,  // braked to rest, then replanned from the stopping point
    settled,               // already within the settle tolerance; a zero-duration profile
};

enum class online_planner_substitution_reason
{
    none,
    reversal_or_overshoot,             // the commanded motion reverses, or would overshoot
    carry_velocity_shape_unavailable,  // the carry-velocity shape does not exist for this state
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

The quantitative fields are what let a supervisory layer choose between axes, or log why a move took longer than it commanded. A boolean could not.

```cpp
if (!planner.update(target).has_value()) {
    // Reject the non-finite command without changing the active profile.
}

if (planner.diagnostics().disposition
    == ctrlpp::online_planner_disposition::braked_and_replanned)
{
    // The axis is on a longer profile than the one commanded. Its extra time is
    // brake_duration, and it restarts from replan_start_position.
}
```

The type is shared with [online-planner-2nd](online-planner-2nd.md), which bounds no jerk and therefore never reports `carry_velocity_shape_unavailable`.

## Usage Example

```cpp
// Usage: ./program | gnuplot -p -e "set datafile separator ','; plot '-' using 1:2 with lines title 'pos', '' using 1:3 with lines title 'vel'"

#include <ctrlpp/trajectory/online_planner_3rd.h>

#include <iostream>

int main()
{
    auto result = ctrlpp::online_planner_3rd<double>::create({
        .v_max = 1.0,
        .a_max = 5.0,
        .j_max = 50.0,
    });
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

- [online-planner-2nd](online-planner-2nd.md)<br/> 2nd-order variant (faster, acceleration-limited only)
- [double-s-trajectory](double-s-trajectory.md)<br/> Pre-computed double-S profile
- [Trajectory Generation Theory](../../background/trajectory-generation.md)<br/> Double-S profile mathematics and jerk limitation
