# Curated fuzz seed corpus

A small, hand-curated regression set for the fuzz targets that have one. One directory per target,
named to match the target binary exactly, which is libFuzzer's own corpus-directory convention:
the directory is handed to the binary with no translation layer. A target without a directory here
is explored and not replayed; adding one is what pins a recorded counterexample.

```
tests/fuzz/corpus/
  fuzz_trapezoidal/
  fuzz_double_s/
  fuzz_ekf/
  fuzz_care/
  fuzz_dare/
```

This is a curated set, not a captured campaign. Each seed is a configuration that pins one named
behavior: a recorded counterexample, a branch boundary, or a typed rejection. Growing coverage
corpora are unbounded binary churn and do not belong in the tree.

## How it is replayed

`scripts/fuzz_smoke.sh` and the fuzz workflow both replay these directories with a zero-run pass
before any exploratory run. The pass is zero-run on purpose: handing a corpus directory to a fuzzing
run makes libFuzzer write newly discovered inputs into it, which would turn a versioned directory
into a corpus dump and silently change the regression set. The exploratory run is never given the
directory.

```sh
scripts/fuzz_smoke.sh build-fuzz
```

## Why the source revision is recorded

A stored byte sequence means nothing on its own. It is meaningful only against the decoder that maps
it to a problem, so a change to a target's input layout silently reinterprets every stored entry
that predates it: the bytes still load, still run, and no longer test what they were written to test.
The provenance table below therefore records, for every seed, the revision it was authored against
and the decoder contract it assumes. Because the corpus lives in the tree it versions with
`LLVMFuzzerTestOneInput` automatically, and the table makes that relationship explicit rather than
leaving it implicit.

The source revision column is the commit that last changed the target's decoder, which is the commit
these bytes were authored against. If a target's `LLVMFuzzerTestOneInput` gains, drops, or reorders
a field, every one of its seeds must be regenerated and this table updated in the same change.

## A decoder revision is not a filter revision, and both can retire a seed

The column above tracks the byte layout, because that is what silently reinterprets a stored entry.
A target's input FILTERS are a second, independent axis: a filter change leaves every byte meaning
exactly what it meant before and can still stop a seed reaching the branch it was written for, which
is the corpus form of a check that cannot fire. A seed that is filtered out replays clean while
covering nothing, and no verdict in the replay log distinguishes the two.

So the rule for a filter change is the same as for a layout change with a different remedy: every
seed of that target is re-verified to still reach its branch, and the ones that no longer do are
re-authored. The two Riccati targets are where these axes currently differ. Both decoders were last
changed in `ef911e5`. The continuous target's entitlement filter was rewritten later, in `1bb2b53`,
which replaced an invertibility test on the weight factor with a detectability test on the pair and
so admits a family the target never explored before; the discrete target's filters have not moved
since `ef911e5`.

## The Riccati field-order contract, stated once

Both algebraic Riccati targets read the same eleven fields in the same order, so one encoder serves
both and the order is a single contract rather than two:

```
11 little-endian binary64, in decode order:
    A(0,0), A(0,1), A(1,0), A(1,1),   the state matrix by rows
    B(0), B(1),                       the input matrix
    Q(0,0), Q(0,1), Q(1,0), Q(1,1),   the raw weight factor by rows
    R                                 the raw input weight
minimum 88 bytes
```

`tools/fuzz_corpus_encode.py` writes exactly this order and carries the same wording in its own
header, so the tool and this document cannot drift apart without one of them being visibly wrong.

A stored field is not the value the target solves with, in two separate ways, and a reader of the
tables below needs both. Each of the first ten is clamped to the bound the decoder states and then
flushed to exact zero below the floor the decoder derives from that bound, so a value under the
floor reaches the solve as nothing at all. The eleventh is clamped and then SQUARED and added to a
derived floor, so the input weight the solve receives is never the number stored in the file.

## Provenance

| target | seed file | what it reproduces | source revision | decoder contract assumed |
|---|---|---|---|---|
| `fuzz_trapezoidal` | `plateau_cruise_above_both.bin` | retiming that lands on a cruise velocity at or above both boundary velocities | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `ramp_through_cruise_between.bin` | retiming that lands on a cruise velocity strictly between the two boundary velocities, where the branch is a hyperbola and not a quadratic | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `valley_cruise_below_both.bin` | retiming that lands on a cruise velocity at or below both boundary velocities, so the axis dips and climbs back | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `valley_equal_boundary_velocities.bin` | the recorded counterexample: equal boundary velocities where the plateau branch's smaller root falls outside its own validity interval and produces a negative acceleration phase, and the valley root is the admissible one | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `rejection_displacement_below_kinetic.bin` | displacement below the distance the two boundary velocities already sweep, so the reachable durations have a finite supremum and a request past it must be a typed rejection that leaves the profile untouched | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `zero_displacement_nonzero_boundary.bin` | nothing to traverse and opposing boundary velocities, asked for a longer duration: a typed rejection, not a duration the profile does not realize | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `valley_reparametrized_boundary_window.bin` | the recorded valley counterexample: an increment several times wider than the whole duration window the valley's own ramp residual leaves it, which the retired parametrization accepted and then realized 7.7e5 units in the last place short of. Now a typed rejection that leaves the profile untouched | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `valley_equal_boundary_long_realization.bin` | the counterexample's equal-boundary neighbor, which IS reachable: the retired parametrization accepted it and realized a duration 3.3e6 units in the last place LONG, where the shape-boundary parametrization lands within a fraction of one. The two seeds together are why an acceptance verdict alone was never evidence -- the retired form got the verdict wrong in one direction here and in the other on its neighbor, reporting success both times | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `valley_ramps_sweep_displacement_exactly.bin` | equal boundary velocities whose two ramps sweep the commanded displacement exactly, so the cruise-velocity form's constant term cancels to a representable zero while its linear coefficient stays negative. The library takes the root selection that ADDS two magnitudes there and resolves the retiming to the last bit; the retired conditioning model charged that selection the constant term's relative error regardless, and an exactly cancelled term makes that error unbounded, so the model called an exact answer a defect | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `plateau_rise_vanishing_boundary_velocity.bin` | a displacement of 1e143 against a larger boundary velocity of 1e-11, which puts the shifted solve's linear coefficient and that coefficient's own error budget 149 and 160 decades out. Their product leaves the representable range while every quantity the model is built from is still ordinary, so this seed pins the conditioning model's divide-before-multiply structure: a model that FORMS the discriminant's error instead of dividing into it reports an unbounded amplification on a retiming the library lands within four parts per million of | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_double_s` | `no_cruise_unequal_boundary.bin` | the recorded counterexample, construction only: a shape with no cruise segment and unequal boundary velocities, the family whose earlier back-off loop could not terminate | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `no_cruise_unequal_boundary_stretched.bin` | the same configuration retimed halfway to its reachable supremum, so the seed also carries the time-scaling leg | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `no_cruise_negative_boundary.bin` | the recorded counterexample whose position was wrong by three orders of magnitude, construction only: no cruise segment, both boundary velocities negative against a negative displacement | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `no_cruise_negative_boundary_stretched.bin` | the same configuration retimed halfway to its reachable supremum | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `rest_to_rest_closed_form.bin` | rest-to-rest retiming, the closed-form path where the scale is a single quotient of two durations | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `equal_boundary_velocities.bin` | equal nonzero boundary velocities, the family the earlier construction happened to survive and which must therefore keep passing | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `rejection_past_reachable_supremum.bin` | a request at twice the reachable supremum: a typed rejection that leaves the profile untouched | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_care` | `marginal_mode_rank_deficient_weight.bin` | the ordinary output-map regulator: a mode exactly on the imaginary axis observed through a rank-one weight factor. The retired invertibility guard rejected every such factor outright, so this whole family was unexplored; the detectability test admits it, and this seed pins that the entitlement test runs at a mode of zero real part and admits there | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `unstable_mode_rank_deficient_weight.bin` | the same rank-one weight factor against a mode of strictly positive real part, so the entitlement test is exercised in the open right half-plane and not only on its boundary | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `stable_pair_full_rank_weight.bin` | every mode strictly stable and the factor full rank, so the entitlement loop tests no mode at all. The complementary branch to the two above: a pose whose filter is vacuous by construction and which must still be judged by the forward error | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `zero_weight_zero_solution_no_verdict.bin` | a weight factor flushed to exact zero against a stable pair, so the returned solution is the exact zero matrix and a RELATIVE forward error has no scale to be relative to. The reference converges and the oracle then declines, which is the largest single abstention class this target's population carries and is neither ill-conditioning nor an exhausted budget | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `refusal_unobservable_unstable_mode.bin` | a mode of positive real part invisible through the weight factor, so the pair is undetectable and no stabilizing solution is unique. The target refuses before it is entitled to demand an answer, which is a typed refusal and not an oracle verdict | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `unstable_mode_weakly_controllable.bin` | the one-parameter family the error enumerator's own reliance statement is measured on, at the weakest input coupling this target's controllability-conditioning filter admits. It pins a pose reaching toward the resolution boundary rather than only well-conditioned ones, and it is 2.4 decades above that boundary because the filter, not the arithmetic, is what stops it going lower | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `stable_pair_full_rank_weight.bin` | an ordinary well-conditioned pose with both modes inside the unit disk, so the forward error against the binary128 reference resolves and passes with the whole solve exercised | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `unstable_mode_well_conditioned.bin` | a mode outside the unit disk on a well-conditioned pair, so the invariant-subspace reordering has a genuine separation to make rather than a trivial one | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `refusal_rank_deficient_weight.bin` | a singular weight factor, refused by the invertibility guard this target still carries and its continuous twin no longer does. The seed pins that asymmetry: the same bytes are a refusal here and an admitted regulator pose there | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `refusal_ill_conditioned_state_matrix.bin` | a state matrix whose singular-value ratio exceeds this target's own bound, refused before the weight factor is read at all. The guard exists because the solve forms the inverse transpose of that matrix as an intermediate, and it has no counterpart in the continuous target | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_ekf` | `rejected_step_poisoned_covariance.bin` | an initial state one finite-difference step short of the largest representable value, so the central-difference stencil's lower sample overflows, one Jacobian entry becomes infinite and the propagated covariance carries a NaN. `predict` is infallible by contract and is allowed to do this; the update that follows then rejects the step naming the carried covariance and leaves it bitwise untouched, which is the contract. The seed pins how that is CHECKED: an elementwise `operator!=` reports a NaN-carrying matrix as different from a byte-for-byte copy of itself, so the target read a met contract as a violation | `0148c3c` | 10 little-endian binary64: `x0(2), z(2), Q_diag(2), R_diag(2), a00, a11`; minimum 80 bytes |

## Decoded field values

Each field is the binary64 nearest the decimal shown. Every value is inside its decoder clamp, so
the decoded configuration is the configuration below with nothing altered; the seed files are exactly
the minimum decoded size, well inside the length cap the smoke run imposes.

`fuzz_trapezoidal`

| seed file | q0 | q1 | v_max | a_max | v0 | v1 | stretch |
|---|---|---|---|---|---|---|---|
| `plateau_cruise_above_both.bin` | 0 | 10 | 3 | 1 | 1.5 | 0.5 | 1.0526315789473684 |
| `ramp_through_cruise_between.bin` | 0 | 10 | 3 | 1 | 1.5 | 0.5 | 2.1052631578947367 |
| `valley_cruise_below_both.bin` | 0 | 10 | 3 | 1 | 1.5 | 0.5 | 5.2631578947368425 |
| `valley_equal_boundary_velocities.bin` | 0 | 1 | 1 | 1 | 0.9 | 0.9 | 9.9009900990099009 |
| `rejection_displacement_below_kinetic.bin` | 0 | 0.9 | 2 | 1 | 1 | 1 | 6.6066802089139518 |
| `zero_displacement_nonzero_boundary.bin` | 2 | 2 | 2 | 1 | 0.5 | -0.5 | 5 |
| `valley_reparametrized_boundary_window.bin` | 0 | 0.000244140625 | 0.9999999999999996 | 1e-6 | 0.9999999999999994 | 0.9999999999999942 | 1.0000000002328306 |
| `valley_equal_boundary_long_realization.bin` | 0 | 0.00390625 | 0.9999999999999821 | 1e-6 | 0.999999999999982 | 0.999999999999982 | 1.0000000009095642 |
| `valley_ramps_sweep_displacement_exactly.bin` | 0 | 1 | 10 | 1 | 1 | 1 | 2 |
| `plateau_rise_vanishing_boundary_velocity.bin` | 0 | 1e143 | 1e6 | 1e6 | -999999 | 1e-11 | 1e6 |

`fuzz_double_s`

| seed file | q0 | q1 | v_max | a_max | j_max | v0 | v1 | stretch |
|---|---|---|---|---|---|---|---|---|
| `no_cruise_unequal_boundary.bin` | 5.911630 | 6.893902 | 5.905937 | 18.316431 | 66.404132 | 2.447245 | 0.193441 | 1 |
| `no_cruise_unequal_boundary_stretched.bin` | 5.911630 | 6.893902 | 5.905937 | 18.316431 | 66.404132 | 2.447245 | 0.193441 | 1.1836876762187425 |
| `no_cruise_negative_boundary.bin` | 1.744413 | -2.633139 | 5.808603 | 13.340846 | 3.971480 | -1.630792 | -3.772260 | 1 |
| `no_cruise_negative_boundary_stretched.bin` | 1.744413 | -2.633139 | 5.808603 | 13.340846 | 3.971480 | -1.630792 | -3.772260 | 1.0140463037925442 |
| `rest_to_rest_closed_form.bin` | 0 | 10 | 3 | 2 | 5 | 0 | 0 | 3 |
| `equal_boundary_velocities.bin` | 0 | 5 | 3 | 2 | 10 | 1 | 1 | 1.5135135135135132 |
| `rejection_past_reachable_supremum.bin` | 0 | 5 | 3 | 2 | 10 | 1 | 1 | 4.0540540540540526 |

`fuzz_ekf`

| seed file | x0(0) | x0(1) | z(0) | z(1) | Q(0,0) | Q(1,1) | R(0,0) | R(1,1) | a00 | a11 |
|---|---|---|---|---|---|---|---|---|---|---|
| `rejected_step_poisoned_covariance.bin` | 0 | -1.7976931348623157e308 | 0 | 0 | 1 | 1 | 1 | 1 | 0.5 | -0.5 |

`x0(1)` is the negative of the largest finite binary64. The stencil's step there is
`cbrt(eps) * |x0(1)|`, so the lower sample leaves the representable range while every other field is
ordinary; nothing about the seed depends on the exact value beyond its being within one such step of
the range's edge.

A `stretch` of one is below the target's own `T_new > T` guard, so those two seeds exercise
construction and the dense scan without entering the time-scaling leg. That is deliberate: it keeps
the construction of each recorded counterexample pinned independently of the retiming.

`fuzz_care`

Decimal, because every stored field is an ordinary named quantity and each is the binary64 nearest
the decimal shown. The last column is derived rather than stored: it is the input weight the solve
actually receives, which is the stored eleventh field squared plus the floor the decoder derives from
the clamp bound.

| seed file | A(0,0) | A(0,1) | A(1,0) | A(1,1) | B(0) | B(1) | Q(0,0) | Q(0,1) | Q(1,0) | Q(1,1) | R stored | R solved with |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `marginal_mode_rank_deficient_weight.bin` | 0 | 1 | 0 | -0.5 | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 1.004 |
| `unstable_mode_rank_deficient_weight.bin` | 0.5 | 1 | 0 | -1 | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 1.004 |
| `stable_pair_full_rank_weight.bin` | -0.5 | 1 | 0 | -1.5 | 0 | 1 | 1 | 0 | 0.5 | 1 | 1 | 1.004 |
| `zero_weight_zero_solution_no_verdict.bin` | -0.5 | 1 | 0 | -1.5 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 1.004 |
| `refusal_unobservable_unstable_mode.bin` | 0.5 | 0 | 0 | -1 | 1 | 1 | 0 | 0 | 0 | 1 | 1 | 1.004 |

Every entry is inside the clamp bound and above the zero floor or exactly zero already, so the
decoded configuration is the configuration above with nothing altered.

`fuzz_dare`

Decimal, on the same convention.

| seed file | A(0,0) | A(0,1) | A(1,0) | A(1,1) | B(0) | B(1) | Q(0,0) | Q(0,1) | Q(1,0) | Q(1,1) | R stored | R solved with |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `unstable_mode_weakly_controllable.bin` | 2 | 0 | 0 | 0.5 | 0.085 | 1 | 1 | 0 | 0 | 1 | 0.99799799498397 | 0.99999999799202421 |
| `stable_pair_full_rank_weight.bin` | 0.5 | 0.2 | 0 | -0.4 | 1 | 0.5 | 1 | 0.3 | 0 | 0.8 | 0.5 | 0.254 |
| `unstable_mode_well_conditioned.bin` | 1.5 | 0.3 | 0 | 0.4 | 1 | 0.5 | 1 | 0.3 | 0 | 0.8 | 0.5 | 0.254 |
| `refusal_rank_deficient_weight.bin` | 0.5 | 0.2 | 0 | -0.4 | 1 | 0.5 | 1 | 0.5 | 2 | 1 | 0.5 | 0.254 |
| `refusal_ill_conditioned_state_matrix.bin` | 2 | 0 | 0 | 0.1 | 1 | 0.5 | 1 | 0.3 | 0 | 0.8 | 0.5 | 0.254 |

`unstable_mode_weakly_controllable.bin` stores an eleventh field chosen so the solve receives an
input weight of one to nine significant figures rather than exactly one, because the field is
squared and offset before use and no stored value reproduces one exactly. The family it belongs to
is the one the discrete error enumerator's reliance statement is measured on, `A = diag(2, 1/2)`,
`B = [d; 1]`, `Q = I`, `R = 1`, whose accepted range ends at `d = 3.6e-4`. **This target cannot
reach that boundary and the reason is its own filter, not its arithmetic.** Its controllability
conditioning bound refuses the family below a measured crossing of `0.0845 < d <= 0.085` -- located
identically under four builds, `clang++ 22.1.8` at `-O0`, `-O2` and `-O3` and `g++ 16.1.1` at `-O2`
-- and the seed sits at `0.085`, the first admitted rung, which is 2.4 decades above the boundary.
Written down here so nobody reads this corpus as pinning the boundary itself.

## Recorded reproducing inputs, which are deliberately NOT seeds

Four inputs that abort this tree's continuous target are recorded here in full and are **not** in
`fuzz_care/`, for one reason: the replay leg is a deterministic regression gate over the curated
directory, so a curated seed that aborts turns that gate red by construction and would replace a
regression signal with a permanent failure. An input that aborts is a defect to triage or an
entitlement question to answer; it is not a corpus entry. They are recorded rather than dropped
because the first of them is the one input a whole oracle rewrite was diagnosed from, and until now
it existed only as an opaque blob.

**They do not all abort at the same place, and the difference is the point.** Three reach the
pre-existing positive semi-definiteness pivot check; each carries a solution that is rank one in
exact arithmetic, so its smallest pivot is pure cancellation, and on all three that pivot is exactly
THREE units in the last place of the solution's largest entry against a floor standing at just over
two. The fourth aborts at the **accuracy** gate instead, on a pose whose closed loop sits `2.05e-9`
from the imaginary axis: there the conditioning exceeds `1/sqrt(eps)`, so a half-significand answer
is not attainable at binary64 by any backward-stable algorithm, and the target has no entitlement
test that would decline to ask for one. Two different questions, and neither is a seed.

Hexadecimal, because these are bit-specific: a decimal rendering of either does not survive the
round trip, and the exact bits are what the diagnosis rested on. The values are the pose AFTER the
decoder's clamp and zero-flush, which is what the target solves; the raw stored fields behind the
clamped entries are not recoverable beyond the bound they exceeded.

| input | A(0,0) | A(0,1) | A(1,0) | A(1,1) | B(0) | B(1) | Q(0,0) | Q(0,1) | Q(1,0) | Q(1,1) | R solved with |
|---|---|---|---|---|---|---|---|---|---|---|---|
| the oracle-rewrite diagnosis input | `-0x1p+1` | `0x0p+0` | `0x1p+1` | `0x0p+0` | `-0x1.fbfbfbf01p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1.fbfbfbfbf60ffp-4` | `-0x1.fbfbfbfbfbfbfp-4` | `0x0p+0` | `0x0p+0` | `0x1.0624dd2f1a9fcp-8` |
| the campaign input | `0x0p+0` | `0x1p+1` | `0x1p+1` | `0x0p+0` | `0x0p+0` | `-0x1.fbfbf02p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1.fbfbfbf002fbfp-4` | `-0x1.fbfbfbfp-4` | `0x1.0624dd2f1a9fcp-8` |
| the accuracy-gate input | `0x0p+0` | `0x0p+0` | `-0x1p+1` | `-0x1p+1` | `-0x1.fbdbfbfbebfbfp-4` | `0x0p+0` | `-0x1.fbfbf31bfbf00p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `0x0p+0` | `0x0p+0` | `0x1.004189374bc6ap+2` |
| the second definiteness input | `0x0p+0` | `0x1p+1` | `0x0p+0` | `-0x1p+1` | `0x0p+0` | `-0x1p+1` | `-0x1.fbfbfbfe3e3e3p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1p+1` | `-0x1p+1` | `0x1.0624dd2f1a9fcp-8` |

The same eleven values per input, one to a line, as the encoder invocation that regenerates each
pose. This is the useful form: it is a command, not a transcription. It reproduces an input the
target decodes to the pose above, which is not the same as reproducing the recorded input's own
bytes -- a clamped entry hides which value above the bound produced it, and the recorded inputs
carry such entries. Both invocations were run and both abort at the same check with the same pivot
as the recorded inputs do.

The oracle-rewrite diagnosis input:

```sh
tools/fuzz_corpus_encode.py --out oracle_rewrite_diagnosis.bin \
  --a00 -0x1p+1 \
  --a01 0x0p+0 \
  --a10 0x1p+1 \
  --a11 0x0p+0 \
  --b0 -0x1.fbfbfbf01p-4 \
  --b1 -0x1.fbfbfbfbfbfbfp-4 \
  --q00 -0x1.fbfbfbfbf60ffp-4 \
  --q01 -0x1.fbfbfbfbfbfbfp-4 \
  --q10 0x0p+0 \
  --q11 0x0p+0 \
  --r 0x0p+0
```

Its solution's four entries are all `0.031622776623…`, its smallest pivot is
`-2.0816681711721685e-17` against a floor of `-1.4043333884132805e-17`, a ratio of `1.4823`, and the
binary128 reference converges in two steps to a relative forward error squared of `6.82e-31`, which
is fifteen decades inside the accuracy criterion. **The answer is accurate and the target aborts on
it.**

The campaign input:

```sh
tools/fuzz_corpus_encode.py --out campaign_input.bin \
  --a00 0x0p+0 \
  --a01 0x1p+1 \
  --a10 0x1p+1 \
  --a11 0x0p+0 \
  --b0 0x0p+0 \
  --b1 -0x1.fbfbf02p-4 \
  --q00 -0x1.fbfbfbfbfbfbfp-4 \
  --q01 -0x1.fbfbfbfbfbfbfp-4 \
  --q10 -0x1.fbfbfbf002fbfp-4 \
  --q11 -0x1.fbfbfbfp-4 \
  --r 0x0p+0
```

Its smallest pivot is `-6.6613381477509392e-16` against a floor of `-4.653561350979949e-16`, a ratio
of `1.4315`. Here the binary128 reference **withdraws its own verdict**, because its own determinant
is negative too, so whether that answer is right is not settled by anything in this tree.

The accuracy-gate input, which is the only one of the four that does not reach the definiteness
check:

```sh
tools/fuzz_corpus_encode.py --out accuracy_gate.bin \
  --a00 0x0p+0 \
  --a01 0x0p+0 \
  --a10 -0x1p+1 \
  --a11 -0x1p+1 \
  --b0 -0x1.fbdbfbfbebfbfp-4 \
  --b1 0x0p+0 \
  --q00 -0x1.fbfbf31bfbf00p-4 \
  --q01 -0x1.fbfbfbfbfbfbfp-4 \
  --q10 0x0p+0 \
  --q11 0x0p+0 \
  --r -0x1p+1
```

Its state matrix has a marginal mode at zero whose eigenvector is `[1, -1]/sqrt(2)`, so the mode's
visibility through the weight factor is `|Q(0,0) - Q(0,1)| = 3.3061954049506959e-08`. The returned
solution places the closed loop at an abscissa of `-2.0478149176383909e-09`, giving a Lyapunov
separation whose reciprocal is `2.44e8` against the `6.71e7` past which a half-significand answer is
unattainable at binary64. The relative forward error squared is `1.2006016725322539e-15`, `2.33`
times the criterion, so the accuracy gate aborts -- while the definiteness pivot is strictly positive
and the closed loop is resolvably stable. The measured forward error `3.5e-8` sits just inside the
backward-stable prediction `kappa * eps = 5.4e-8`, so **the answer is at the arithmetic's own bound
and the criterion is below it.** Walking `Q(0,0)` across 41 units in the last place at that same
visibility, the solver declines on 37 rungs and answers on 4, and all four abort here.

The second definiteness input, which is the same class as the first two and was found independently:

```sh
tools/fuzz_corpus_encode.py --out second_definiteness.bin \
  --a00 0x0p+0 \
  --a01 0x1p+1 \
  --a10 0x0p+0 \
  --a11 -0x1p+1 \
  --b0 0x0p+0 \
  --b1 -0x1p+1 \
  --q00 -0x1.fbfbfbfe3e3e3p-4 \
  --q01 -0x1.fbfbfbfbfbfbfp-4 \
  --q10 -0x1p+1 \
  --q11 -0x1p+1 \
  --r 0x0p+0
```

Its smallest pivot is `-4.163336342344337e-17` against a floor of `-2.8140615587208564e-17`, a ratio
of `1.4795`. Its solution's largest entry is `0.06336703293626443`, whose unit in the last place is
`1.3877787807814457e-17`, so the pivot is **exactly three** of those units and the floor stands at
`2.0277` of them. The binary128 reference withdraws its own verdict here, as it does on the campaign
input, and the forward error computed for triage is `1.03e-30`, 14.3 decades inside the criterion.

**All four invocations were run and each reproduces its recorded pose, pivot, floor and forward error
exactly.** The four inputs the bytes were recovered from are 88, 104, 102 and 95 bytes. The decoder
copies exactly 88 and ignores anything past them, so the trailing bytes of the longer ones carry
nothing; a mutation length is not a field.

## Adding a seed

Author the bytes from named configuration values with `tools/fuzz_corpus_encode.py` rather than
assembling them by hand: a hand-assembled seed is unverifiable and undebuggable. The tool takes
values by NAME so two fields cannot be silently transposed, accepts hexadecimal float input as well
as decimal because a recorded input is bit-specific, and reads the file back after writing to report
the decoded configuration against the target's own clamps, so checking the emitted values is a step
it performs rather than one a reader may skip. It covers the two algebraic Riccati targets, which
share one field order; a target with a different layout needs the same treatment and not a throwaway
script.

Keep the file at or above the target's minimum decoded size and below the smoke run's length cap,
give it a descriptive name so a failing replay reads clearly, and add its row above.

Two more things, neither of which the replay log will tell you. **Verify the seed actually reaches
the branch it was written for**, because a filtered-out seed replays clean while covering nothing.
**Do not author a seed that is expected to abort**, because the replay leg is a deterministic
regression gate and every curated seed must replay without a finding.
