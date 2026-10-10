# Planner overshoot sweep

Probe: do the online planners' overshoot verdicts and the chain each
one reports against still disagree anywhere, and over what population? Does
every retargeted plan stay inside its limits, and does each shape
constructor honor its boundary-velocity precondition exactly?

The program runs two families of legs. The agreement legs, described first,
compare two booleans the library returns. The limit legs, described after
them, compare sampled values against the limits the caller stated.

Each planner decides whether a command overshoots, and then either carries
its current velocity into a profile or brakes to rest first. Each verdict is
meant to share one arithmetic chain with what it reports against:

  * the jerk-limited planner's overshoot verdict calls the double-S shape's
    own minimum-displacement predicate on the operands the shape's
    constructor sees, so the planner and the constructor answer the same
    question with the same arithmetic. A target the planner admits that the
    shape still refuses is a disagreement, and the planner reports it as an
    unavailable carry-velocity shape;
  * the trapezoidal planner's overshoot verdict and the stopping point it
    places are one braking distance, the mean speed times the braking
    duration. This leg's reference is that documented formula, so it checks
    the implementation against the formula rather than against a second
    chain.

Where a verdict and its reference differ lies a band of commanded targets
the two classify differently. The planner reports every substitution, so a
band would be public behavior rather than an internal detail.

This program measures that band. It compares no quantity against a tolerance
of its own: both verdicts are booleans the library itself returns, and the
comparison between them is exact. A threshold inside the instrument would be
exactly the kind of underived constant this measurement exists to locate.

**What is measured, per configuration**

  1. Each boundary is LOCATED rather than assumed, by bisecting the
     commanded target in the space of values the scalar type represents. It
     comes out as the smallest representable target at which its verdict
     flips, so the distance between the two is a count of representable
     targets -- units in the last place of the commanded target -- and not a
     length that has to be converted into one.
  2. Every representable target from a fixed number of units in the last
     place below the lower boundary to the same number above the upper one
     is swept exhaustively. Nothing between the two is stepped over: the
     band is walked one representable value at a time. A walked target at or
     behind the start position asks the direction test's question, not the
     overshoot test's, so it is left out of the agreement count, counted on
     its own line, and required to report reversal_or_overshoot.
  3. A geometric ladder of offsets, out to 2^(digits - 9) units in the last
     place (about thirteen decades in double, four and a half in float) and
     eleven decades of the stopping distance in both directions, establishes
     whether the disagreement is confined to the band or is a property of
     the whole domain.

Both agreement legs run in double and in float. In float a unit in the last
place of an axis position near 1e4 is about 1e-3, coarse enough that the
walk's margin reaches behind the start position; in double it does not.

A verdict that does not flip anywhere in the widened bracket is REPORTED as
not flipping rather than dropped. A shape that refuses every commanded
displacement is a disagreement over the entire admissible domain, which is a
larger finding than a narrow band and would be invisible to an instrument
that required both boundaries to exist before it recorded anything.

The swept and disagreeing counts are exact numbers, and every disagreeing
point is printed with the configuration that produced it.

**The planner's verdict is read from the reason it reports, not from whether
it substituted.** A substitution reported as an unavailable carry-velocity
shape IS the shape's refusal showing through the planner, so reading the
disposition would compare the shape's verdict against itself and find
perfect agreement by construction.

**The two legs are not the same comparison.** The jerk-limited planner
calls the double-S constructor directly, so its leg compares the planner's
verdict against that constructor's own rejection on exactly the command the
planner hands it. The trapezoidal planner never calls that constructor -- it
bounds no jerk and has no shape with a domain to fall outside of -- so its
leg compares the planner's verdict against the braking distance formula the
planner documents for both its decision and its stopping point.

**Every state the sweep probes is one a caller can present.** The jerk
limited planner is driven to its own cruise phase and probed there, where
the sampled acceleration is exactly zero, so the acceleration-nulling phase
does not run and the command reaches the overshoot test carrying the
position and speed the planner itself reported. The shape constructor is
then handed those same two numbers. The trapezoidal planner carries no
acceleration state, so it is probed part way up its acceleration ramp and
its initial speed spans the reachable range. Neither leg reaches into a
private member or reconstructs one.

**Retargets from nonzero acceleration.** A third agreement leg, in double
and float, retargets the jerk-limited planner forward from five instants in
each of its acceleration-ramp and acceleration-nulling windows, on the eight
limit sets of the limit legs. There the nulling phase runs first and the
state it leaves behind is private, so the leg compares no second chain: it
locates the planner's overshoot boundary on the reason it reports, walks and
ladders around it, and counts the shape's refusal as the planner reports it
(carry_velocity_shape_unavailable) and any substitution at or above the
boundary. Every probe drives a planner afresh from rest to the same instant.

**The limit legs**

  * Retarget envelope, in double and float. The jerk-limited planner is
    driven from rest toward a far target and retargeted further along the
    same direction from 40 instants in each of its cruise, acceleration-ramp
    and acceleration-nulling windows, on eight limit sets. Every sample of
    the new plan is checked against v_max and a_max: densely over its first
    and last twenty seconds, where every phase but the cruise lives, and once
    a second between. Every dense pair of samples is checked against j_max
    times their spacing plus an allowance built from two counted rounding
    chains, one on elapsed time and one on acceleration. That allowance is
    printed beside j_max dt for each row, so the resolution of the jerk check
    is visible: in float its time term dominates once elapsed time reaches
    the order of a thousand seconds. A same-direction retarget that reports
    any substitution is counted as a brake to rest that was not needed.
  * Reversing, in double and float. The same planner is retargeted behind
    its current position from the cruise and acceleration-ramp windows. It
    must report a reversal and keep the same envelope while it brakes and
    replans.
  * Boundary-velocity precondition, for both the double-S and the
    trapezoidal constructor in double and float. A boundary speed is
    stepped from three representable values below v_max to six above it, as
    the initial or the final velocity, in both directions of travel, under a
    boundary-velocity rounding count of 0, 1, 2 and 4. A speed at or below
    the limit must build. A speed above it must build only within the
    admitted excess, and the clamped end must then report exactly +-v_max.
    The realized peak, the double-S shape's reported peak or the trapezoidal
    shape's peak sampled on a millisecond grid, must not exceed the limit.

Like the agreement legs, the limit legs reach no private member: they use
only create, update, sample, diagnostics, evaluate and peak_velocity.

**The phase-switch legs**

These legs sample every representable instant within 64 steps of each phase
switch, for the trapezoidal and double-S shapes and for both online
planners, in double and float. Each sampled speed is compared against
v_max. Each sampled acceleration is compared against a_max, or for the
trapezoid against the realized acceleration its disposition reports.
Windows that overlap are merged, so every instant is sampled once.

Rounding lands a sample above a limit only in particular regimes:

  * a boundary velocity that is negative, or below half the limit;
  * an end velocity at the limit, pointing away from the target;
  * a jerk phase near or below the time resolution of a long profile;
  * a ramp at the threshold where it stops reaching a_max.

The populations aim at those regimes rather than at fixed fractions.

  * Trapezoidal: the eight limit sets, both directions, seven fixed and 512
    seeded pairs of boundary fractions in [-1, 1) over a seeded displacement.
    Every profile is also rescaled to 1.5 and 3 times its duration where
    can_rescale_to allows it.
  * Double-S: the same pairs over a displacement that admits both cruise and
    no-cruise shapes. It adds long profiles of about 2^m seconds, for three m
    taken from the type's significand width, with a jerk phase of 1/4 to 2
    units of the time resolution there and an end velocity of +-v_max. It
    also adds a boundary speed walked one representable value at a time across
    a_max^2 / j_max.
  * Second-order planner: same-direction retargets from a seeded instant in
    the first half of the first plan's ramp, where the carried speed is below
    v_max / 2. The limits are the eight pairs plus one pair with full-width
    significands; on the short ones v_max - v_r is often exact. The switch is
    placed at t_r + (v_max - |v_r|) / a_max. A window in which the sampled
    acceleration does not go from +-a_max to 0 counts as missing the switch.
  * Third-order planner: the eight sets with j_max set to a_max over 1/4 to
    4 units of the time resolution at the end of the ramp. Each set runs
    three plans, all driven from rest:
      - the move from rest itself;
      - a retarget from a seeded cruise instant;
      - a retarget from a seeded instant on the acceleration ramp.
    Switches come from a double-S shape built on the same command from the
    sampled state. After a ramp retarget, the acceleration-nulling phase is
    first rebuilt in closed form. The planner trims that phase in private
    arithmetic the closed form does not mirror, so a window placed this way
    that shows no switch is counted as not located. That count is reported,
    and is not a zero claim. A window placed from public values alone that
    shows no switch counts as missing it, which is a zero claim.

Every boundary fraction and magnitude is formed from the raw output of a
seeded std::mt19937, whose sequence the standard fixes, as a 24-bit dyadic
value that float and double both hold exactly. The standard distributions
are not used, because their output differs between standard libraries. Each
block prints its seed. Full mode prints every violating sample with its
configuration and instant in hexadecimal floating point, its value, the
limit, and the excess in units in the last place. It also prints every
window that missed or was not located.

Building: standalone, no test framework and no build system. Eigen enters as
a system include, the way the library's own build treats it, so the only
diagnostics a build reports are the tool's own.

```sh
g++ -std=c++20 -O2 -Wall -Wextra -fno-exceptions -fno-rtti -I lib/ctrlpp/include -isystem /usr/include/eigen3 tools/planner_overshoot_sweep/main.cpp -o build/tools/planner_overshoot_sweep
```

Usage: no argument prints the full report. `--summary` prints the per
configuration measurements and the totals without the per-point lines.

**Where the boundaries this program locates were located.** A crossing that
moves with instruction selection is a crossing that has to name the
toolchain it was found on. The full report, the agreement, limit and
phase-switch legs alike, was taken from eight builds on x86-64 Linux against
Eigen 3.4.1, each with -Werror and without a warning, and came out
byte-identical from all of them, in summary and in full mode: g++ 16.2.1 at
-O0, -O2 and -O3, and at -O2 with -ffp-contract=off and with
-ffp-contract=fast; clang 23.1.1 at -O2 with and without -ffp-contract=fast;
and clang 18.1.8 at -O2 against libstdc++ 15.3.0. On every build every
zero claim reads zero: the disagreeing totals, the shape refusals and
substitutions at or above a located boundary, the precondition violations,
the windows missing their switch, every sample above a limit and every jerk
pair above its allowance. The third-order count of windows not located is
reported, not claimed zero, and reads 9 in double and 5 in float. Eigen
takes no part in any verdict, since neither planner nor any shape
constructor forms a matrix.
