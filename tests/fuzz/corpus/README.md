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
  fuzz_online_planner_3rd/
  fuzz_online_planner_3rd_float/
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
| `fuzz_trapezoidal` | `boundary_velocity_at_limit.bin` | a boundary velocity equal to the limit, which builds a profile that cruises at the limit | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_double_s` | `no_cruise_unequal_boundary.bin` | the recorded counterexample, construction only: a shape with no cruise segment and unequal boundary velocities, the family whose earlier back-off loop could not terminate | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `no_cruise_unequal_boundary_stretched.bin` | the same configuration retimed halfway to its reachable supremum, so the seed also carries the time-scaling leg | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `no_cruise_negative_boundary.bin` | the recorded counterexample whose position was wrong by three orders of magnitude, construction only: no cruise segment, both boundary velocities negative against a negative displacement | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `no_cruise_negative_boundary_stretched.bin` | the same configuration retimed halfway to its reachable supremum | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `rest_to_rest_closed_form.bin` | rest-to-rest retiming, the closed-form path where the scale is a single quotient of two durations | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `equal_boundary_velocities.bin` | equal nonzero boundary velocities, the family the earlier construction happened to survive and which must therefore keep passing | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `rejection_past_reachable_supremum.bin` | a request at twice the reachable supremum: a typed rejection that leaves the profile untouched | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `boundary_velocity_at_limit.bin` | a boundary velocity equal to the limit, which builds a profile that cruises at the limit | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_care` | `marginal_mode_rank_deficient_weight.bin` | the ordinary output-map regulator: a mode exactly on the imaginary axis observed through a rank-one weight factor. The retired invertibility guard rejected every such factor outright, so this whole family was unexplored; the detectability test admits it, and this seed pins that the entitlement test runs at a mode of zero real part and admits there | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `unstable_mode_rank_deficient_weight.bin` | the same rank-one weight factor against a mode of strictly positive real part, so the entitlement test is exercised in the open right half-plane and not only on its boundary | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `stable_pair_full_rank_weight.bin` | every mode strictly stable and the factor full rank, so the entitlement loop tests no mode at all. The complementary branch to the two above: a pose whose filter is vacuous by construction and which must still be judged by the forward error | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `zero_weight_zero_solution_no_verdict.bin` | a weight factor flushed to exact zero against a stable pair, so the returned solution is the exact zero matrix and a RELATIVE forward error has no scale to be relative to. The reference converges and the oracle then declines, which is the largest single abstention class this target's population carries and is neither ill-conditioning nor an exhausted budget | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `refusal_unobservable_unstable_mode.bin` | a mode of positive real part invisible through the weight factor, so the pair is undetectable and no stabilizing solution is unique. The target refuses before it is entitled to demand an answer, which is a typed refusal and not an oracle verdict | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `definiteness_rank_one_diagnosis_pose.bin` | the diagnosis input below: a solution that is rank one in exact arithmetic, so the definiteness check's smallest pivot is pure cancellation, on a pose the binary128 reference converges on. The check must accept a matrix within one rounding per entry of a positive semi-definite one, and the solution must reach the target as an exact power-of-two multiple of the one the solver's check accepted | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `definiteness_rank_one_campaign_pose.bin` | the campaign input below: the same rank-one cancellation on a pose where the binary128 reference does not converge, so there is no accuracy verdict and the definiteness check and the closed-loop spectrum are the checks the pose reaches | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_care` | `definiteness_rank_one_seedless_pose.bin` | the second definiteness input below, found independently of the other two: the same rank-one cancellation on a different state and input matrix, again with no converged reference | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `unstable_mode_weakly_controllable.bin` | the one-parameter family the error enumerator's own reliance statement is measured on, at the weakest input coupling this target's controllability-conditioning filter admits. It pins a pose reaching toward the resolution boundary rather than only well-conditioned ones, and it is 2.4 decades above that boundary because the filter, not the arithmetic, is what stops it going lower | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `stable_pair_full_rank_weight.bin` | an ordinary well-conditioned pose with both modes inside the unit disk, so the forward error against the binary128 reference resolves and passes with the whole solve exercised | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `unstable_mode_well_conditioned.bin` | a mode outside the unit disk on a well-conditioned pair, so the invariant-subspace reordering has a genuine separation to make rather than a trivial one | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `refusal_rank_deficient_weight.bin` | a singular weight factor, refused by the invertibility guard this target still carries and its continuous twin no longer does. The seed pins that asymmetry: the same bytes are a refusal here and an admitted regulator pose there | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_dare` | `refusal_ill_conditioned_state_matrix.bin` | a state matrix whose singular-value ratio exceeds this target's own bound, refused before the weight factor is read at all. The guard exists because the solve forms the inverse transpose of that matrix as an intermediate, and it has no counterpart in the continuous target | `ef911e5` | 11 little-endian binary64: `A(0,0), A(0,1), A(1,0), A(1,1), B(0), B(1), Q(0,0), Q(0,1), Q(1,0), Q(1,1), R`; minimum 88 bytes |
| `fuzz_ekf` | `rejected_step_poisoned_covariance.bin` | an initial state one finite-difference step short of the largest representable value, so the central-difference stencil's lower sample overflows, one Jacobian entry becomes infinite and the propagated covariance carries a NaN. `predict` is infallible by contract and is allowed to do this; the update that follows then rejects the step naming the carried covariance and leaves it bitwise untouched, which is the contract. The seed pins how that is CHECKED: an elementwise `operator!=` reports a NaN-carrying matrix as different from a byte-for-byte copy of itself, so the target read a met contract as a violation | `0148c3c` | 10 little-endian binary64: `x0(2), z(2), Q_diag(2), R_diag(2), a00, a11`; minimum 80 bytes |
| `fuzz_online_planner_3rd` | `retarget_cruise_unit_limits.bin` | the recorded counterexample: a retarget to a farther target while cruising at the velocity limit on unit limits, where the earlier replan kept accelerating from the cruise state and exceeded the velocity limit on its first sample | `d8b2b04` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes |
| `fuzz_online_planner_3rd` | `retarget_cruise_review_limits.bin` | the same cruise retarget on non-unit limits whose acceleration phase reaches the acceleration limit, so the cruise is entered through a constant-acceleration phase | `d8b2b04` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes |
| `fuzz_online_planner_3rd` | `retarget_ramp_unit_limits.bin` | a retarget at half the acceleration limit while the jerk is still raising it, so the replan starts from a nonzero acceleration and must null it first | `d8b2b04` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes |
| `fuzz_online_planner_3rd` | `retarget_nulling_unit_limits.bin` | a retarget at half the acceleration limit while the jerk is bringing it back to zero, the phase that ends at the velocity limit | `d8b2b04` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes |
| `fuzz_online_planner_3rd` | `cr04_constant_phase_end.bin` | the recorded counterexample: a retarget from the acceleration ramp on a jerk limit large enough that the replan's closing jerk phase is shorter than the resolution of the sample times, so the last sample in the constant-acceleration phase before it read one unit in the last place above the velocity limit. It aborted at `1f2ea88`, before that phase was evaluated backward from its analytic end | `c45c886` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes |
| `fuzz_online_planner_3rd_float` | `cr04_constant_phase_end.bin` | the same counterexample in `float` on ordinary limits, where the closing jerk phase lasts less than one sample-time resolution near t = 1809 s. It aborted at `1f2ea88` | `c45c886` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes, each field narrowed to `float` after the clamps |
| `fuzz_online_planner_3rd_float` | `retarget_in_place_large_limits.bin` | a retarget after the first sample to the position the axis still holds, on limits of `1e6` and the jerk floor `1e-3`, while the first jerk phase has raised the acceleration to 2001 times `j_max dt`. The replan was taken as a no-op and the acceleration dropped to zero in one sample. It aborted at `361f725` | `c45c886` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes, each field narrowed to `float` after the clamps |
| `fuzz_online_planner_3rd_float` | `retarget_in_place_negative_position.bin` | the same retarget in place at a negative position with the opposite acceleration sign. It aborted at `361f725` | `c45c886` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes, each field narrowed to `float` after the clamps |
| `fuzz_online_planner_3rd_float` | `retarget_in_place_small_acceleration_limit.bin` | the same retarget in place on `a_max = 29.2`, where the dropped acceleration was 478 times `j_max dt`. It aborted at `361f725` | `c45c886` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes, each field narrowed to `float` after the clamps |
| `fuzz_online_planner_3rd_float` | `retarget_in_place_unit_scale_position.bin` | the same retarget in place at position 2.125. It aborted at `361f725` | `c45c886` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes, each field narrowed to `float` after the clamps |
| `fuzz_online_planner_3rd_float` | `retarget_in_place_coarse_sample_period.bin` | the same retarget in place at a sample period of 5.3e-4 s, near position 2.26e6. It aborted at `361f725` | `c45c886` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes, each field narrowed to `float` after the clamps |
| `fuzz_online_planner_3rd_float` | `retarget_in_place_moderate_limits.bin` | the same retarget in place on `v_max = 130` and `a_max = 32.5`. It aborted at `361f725` | `c45c886` | 9 little-endian binary64: `pos, target, v_max, a_max, j_max, dt, t_start, retarget_fraction, retarget_target`; minimum 72 bytes, each field narrowed to `float` after the clamps |

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
| `boundary_velocity_at_limit.bin` | 0 | 10 | 1 | 1 | 1 | 0 | 1 |

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
| `boundary_velocity_at_limit.bin` | 0 | 10 | 1 | 1 | 1 | 1 | 0 | 1 |

`fuzz_ekf`

| seed file | x0(0) | x0(1) | z(0) | z(1) | Q(0,0) | Q(1,1) | R(0,0) | R(1,1) | a00 | a11 |
|---|---|---|---|---|---|---|---|---|---|---|
| `rejected_step_poisoned_covariance.bin` | 0 | -1.7976931348623157e308 | 0 | 0 | 1 | 1 | 1 | 1 | 0.5 | -0.5 |

`fuzz_online_planner_3rd`

| seed file | pos | target | v_max | a_max | j_max | dt | t_start | retarget_fraction | retarget_target |
|---|---|---|---|---|---|---|---|---|---|
| `retarget_cruise_unit_limits.bin` | 0 | 1000 | 1 | 1 | 1 | 0.01 | 0 | 0.5 | 2000 |
| `retarget_cruise_review_limits.bin` | 0 | 700 | 0.7 | 1.3 | 2.9 | 0.005 | 0 | 0.5 | 1400 |
| `retarget_ramp_unit_limits.bin` | 0 | 1000 | 1 | 1 | 1 | 0.01 | 0 | 0.025 | 2000 |
| `retarget_nulling_unit_limits.bin` | 0 | 1000 | 1 | 1 | 1 | 0.01 | 0 | 0.075 | 2000 |
| `cr04_constant_phase_end.bin` | 0 | 0x1.095967f30a7fep+17 | 0x1.2c6932p+4 | 0x1.541afep-7 | 0x1.541afep+36 | 0x1.ea91a70c48f32p-1 | 0 | 0x1.808b5acbe418fp-3 | 0x1.095967f30a7fep+18 |

`fuzz_online_planner_3rd_float`

| seed file | pos | target | v_max | a_max | j_max | dt | t_start | retarget_fraction | retarget_target |
|---|---|---|---|---|---|---|---|---|---|
| `cr04_constant_phase_end.bin` | 0 | 0x1.095968ae6013p+17 | 0x1.2c6932p+4 | 0x1.541afep-7 | 0x1.10b2b4p+7 | 0x1.d25bacp-1 | 0 | 0x1.88bb00858be25p-6 | 0x1.095968ae6013p+18 |
| `retarget_in_place_large_limits.bin` | 0x1.fffff0affffffp+16 | 0x0.000007c00470ap-1022 | -0x1.0000000000000p+1009 | -0x1.000f5000027ffp+865 | 0x1.f0aff40b9b8ffp-1008 | 0x0.0c2bf4647000ap-1022 | 0x1.4e00129000000p+105 | 0x1.0ffffb5be0000p-959 | 0x1.fffff00000000p+16 |
| `retarget_in_place_negative_position.bin` | -0x1.d26e7e7215b00p-13 | -0x1.6fff9ee1918f5p-79 | -0x1.7e70ae7000000p+639 | -0x1.ae7c7e7e7e7e7p+625 | 0x0.140000001e7e7p-1022 | 0x1.5b80000000000p-608 | 0x1.0000000000000p+129 | 0x0.00001e6e7e70ap-1022 | -0x1.d26e7e7000000p-13 |
| `retarget_in_place_small_acceleration_limit.bin` | 0x1.f3f3f3f373f0ep-12 | 0x1.03f3ffb3f3f3fp-11 | -0x1.0c0c0c0c0c13fp+13 | 0x1.d3f3f3f3f3faep+4 | 0x1.f3f3f3f3f3f3fp-12 | 0x0.00000003f003fp-1022 | 0x1.f3f3f3f3f0000p-12 | 0x1.f3f3f3f4a3f3fp-12 | 0x1.f3f3f3f3f403dp-12 |
| `retarget_in_place_unit_scale_position.bin` | 0x1.100000000000ap+1 | 0x0.0000001000000p-1022 | 0x1.0010000040000p+241 | 0x1.5555555555555p+342 | 0x0.0000000000055p-1022 | 0x0.000ff00ffff00p-1022 | 0x1.777774a016700p+888 | 0x0.0000000007777p-1022 | 0x1.100000000000ap+1 |
| `retarget_in_place_coarse_sample_period.bin` | 0x1.1414141404341p+21 | 0x1.14541c0414141p+21 | 0x1.1414141414141p+29 | 0x1.100410a07cf00p+21 | 0x0.0000000000001p-1022 | -0x1.1415b41410041p-11 | 0x1.1414141bebebep+21 | -0x1.ebeb641414141p-20 | 0x1.1414141414141p+21 |
| `retarget_in_place_moderate_limits.bin` | 0x1.040404040402dp+5 | 0x1.0404040404040p-299 | 0x1.0402c4040402cp+7 | 0x1.0404040403940p+5 | -0x1.4a34040404040p-871 | 0x1.dfd00007b2192p-317 | 0x1.0404240404040p+5 | 0x1.14040404040c1p-493 | 0x1.0404040404021p+5 |

The two `cr04_constant_phase_end.bin` seeds are recorded inputs, so their fields are given in hexadecimal, exactly
as stored. The `float` target narrows each field after the clamps, so it plans with target `0x1.095968p+17`,
retarget target `0x1.095968p+18` and the other fields unchanged. The retarget follows sample 47 in `float` and
sample 375 in `double`, on the ramp, and both jerk limits are inside their targets' clamps.

The six `retarget_in_place_*.bin` seeds are recorded inputs as well, stored byte for byte with the trailing bytes the
fuzzer left. Each clamps `j_max` to `1e-3`. Five clamp `t_start` to `2000 dt`, so their first sample lies 2001 sample
periods into the first jerk phase; `retarget_in_place_small_acceleration_limit.bin` starts at 477 sample periods. Each retargets after that sample (`retarget_fraction` truncates to index 0) to its own
starting position, which the axis still holds to the position's resolution.

The retarget is applied after the scan sample whose index is `retarget_fraction * 1999`, truncated, so
the four seeds replan at t = 10, 5, 0.5 and 1.5 s. Sampled there, the planner reports v = 1 and a = 0,
v = 0.7 and a = 0, v = 0.125 and a = 0.5, and v = 0.875 and a = 0.5: cruise, cruise, the rising
acceleration ramp and the nulling phase. `tools/fuzz_corpus_encode.py` covers only the Riccati
layout, so these were written from the named values with a one-line `struct.pack`, for example:

```sh
python3 -c "import struct; open('retarget_cruise_unit_limits.bin', 'wb').write(struct.pack('<9d', 0, 1000, 1, 1, 1, 0.01, 0, 0.5, 2000))"
```

`x0(1)` is the negative of the largest finite binary64. The stencil's step there is
`cbrt(eps) * |x0(1)|`, so the lower sample leaves the representable range while every other field is
ordinary; nothing about the seed depends on the exact value beyond its being within one such step of
the range's edge.

A `stretch` of one is below the target's own `T_new > T` guard, so every seed that stores it exercises
construction and the dense scan without entering the time-scaling leg. That is deliberate: it keeps
the construction of each recorded counterexample, and of a boundary velocity at the limit, pinned
independently of the retiming.

`fuzz_care`

Decimal for the first five, because every stored field there is an ordinary named quantity and each
is the binary64 nearest the decimal shown. Hexadecimal for the three `definiteness_rank_one_*.bin`
seeds, because they are recorded inputs and bit-specific, exactly as stored. The last column is
derived rather than stored: it is the input weight the solve actually receives, which is the stored
eleventh field squared plus the floor the decoder derives from the clamp bound.

| seed file | A(0,0) | A(0,1) | A(1,0) | A(1,1) | B(0) | B(1) | Q(0,0) | Q(0,1) | Q(1,0) | Q(1,1) | R stored | R solved with |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `marginal_mode_rank_deficient_weight.bin` | 0 | 1 | 0 | -0.5 | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 1.004 |
| `unstable_mode_rank_deficient_weight.bin` | 0.5 | 1 | 0 | -1 | 0 | 1 | 1 | 0 | 0 | 0 | 1 | 1.004 |
| `stable_pair_full_rank_weight.bin` | -0.5 | 1 | 0 | -1.5 | 0 | 1 | 1 | 0 | 0.5 | 1 | 1 | 1.004 |
| `zero_weight_zero_solution_no_verdict.bin` | -0.5 | 1 | 0 | -1.5 | 0 | 1 | 0 | 0 | 0 | 0 | 1 | 1.004 |
| `refusal_unobservable_unstable_mode.bin` | 0.5 | 0 | 0 | -1 | 1 | 1 | 0 | 0 | 0 | 1 | 1 | 1.004 |
| `definiteness_rank_one_diagnosis_pose.bin` | `-0x1p+1` | `0x0p+0` | `0x1p+1` | `0x0p+0` | `-0x1.fbfbfbf01p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1.fbfbfbfbf60ffp-4` | `-0x1.fbfbfbfbfbfbfp-4` | `0x0p+0` | `0x0p+0` | `0x0p+0` | `0x1.0624dd2f1a9fcp-8` |
| `definiteness_rank_one_campaign_pose.bin` | `0x0p+0` | `0x1p+1` | `0x1p+1` | `0x0p+0` | `0x0p+0` | `-0x1.fbfbf02p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1.fbfbfbf002fbfp-4` | `-0x1.fbfbfbfp-4` | `0x0p+0` | `0x1.0624dd2f1a9fcp-8` |
| `definiteness_rank_one_seedless_pose.bin` | `0x0p+0` | `0x1p+1` | `0x0p+0` | `-0x1p+1` | `0x0p+0` | `-0x1p+1` | `-0x1.fbfbfbfe3e3e3p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1p+1` | `-0x1p+1` | `0x0p+0` | `0x1.0624dd2f1a9fcp-8` |

Every entry is within the clamp bound, the bound itself included, and above the zero floor or exactly
zero already, so the decoded configuration is the configuration above with nothing altered.

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

## Recorded inputs

Four inputs found by exploring the continuous target are recorded here in full. Three are curated
seeds; the fourth is recorded and is deliberately not a seed.

**The three definiteness inputs are seeds.** Each carries a solution that is rank one in exact
arithmetic, so the smallest pivot of a definiteness check on it is pure cancellation. The solver
divides the weights by an exact power of two, decides definiteness on the equilibrated solution with a
check that accepts every matrix within one rounding per entry of a positive semi-definite one, and
returns an exact power-of-two multiple of the solution that check accepted. The target applies the
same check to the returned solution, and an exact power-of-two scaling does not change its verdict,
so the two agree. All three reach that check and replay clean, which makes them its regression
inputs: they are the `definiteness_rank_one_*.bin` rows of the tables above.

**The accuracy-gate input is not a seed.** The library declines the pose, and the target's
conditioning entitlement excludes it independently: the condition number `care_condition_number`
estimates for it (SLICOT SB02QD), times eps, exceeds sqrt(eps), so a half-significand forward error is not attainable at binary64 by a
backward-stable algorithm and the target does not ask for one. `tests/unit/care_entitlement_test.cpp`
asserts that classification from the binary128 reference alone. As a seed it would replay the
solver's decline and nothing else, so it stays a recorded input that the unit test pins.

Hexadecimal, because these are bit-specific: a decimal rendering does not survive the round trip.
The values are the pose AFTER the decoder's clamp and zero-flush, which is what the target solves;
the raw stored fields behind the clamped entries are not recoverable beyond the bound they exceeded.

| input | A(0,0) | A(0,1) | A(1,0) | A(1,1) | B(0) | B(1) | Q(0,0) | Q(0,1) | Q(1,0) | Q(1,1) | R solved with |
|---|---|---|---|---|---|---|---|---|---|---|---|
| the diagnosis input | `-0x1p+1` | `0x0p+0` | `0x1p+1` | `0x0p+0` | `-0x1.fbfbfbf01p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1.fbfbfbfbf60ffp-4` | `-0x1.fbfbfbfbfbfbfp-4` | `0x0p+0` | `0x0p+0` | `0x1.0624dd2f1a9fcp-8` |
| the campaign input | `0x0p+0` | `0x1p+1` | `0x1p+1` | `0x0p+0` | `0x0p+0` | `-0x1.fbfbf02p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1.fbfbfbf002fbfp-4` | `-0x1.fbfbfbfp-4` | `0x1.0624dd2f1a9fcp-8` |
| the accuracy-gate input | `0x0p+0` | `0x0p+0` | `-0x1p+1` | `-0x1p+1` | `-0x1.fbdbfbfbebfbfp-4` | `0x0p+0` | `-0x1.fbfbf31bfbf00p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `0x0p+0` | `0x0p+0` | `0x1.004189374bc6ap+2` |
| the second definiteness input | `0x0p+0` | `0x1p+1` | `0x0p+0` | `-0x1p+1` | `0x0p+0` | `-0x1p+1` | `-0x1.fbfbfbfe3e3e3p-4` | `-0x1.fbfbfbfbfbfbfp-4` | `-0x1p+1` | `-0x1p+1` | `0x1.0624dd2f1a9fcp-8` |

The same eleven values per input, one to a line, as the encoder invocation that regenerates each
pose. It reproduces an input the target decodes to the pose above, which is not the same as
reproducing the recorded input's own bytes: a clamped entry hides which value above the bound
produced it. The three definiteness invocations write the curated seed files themselves.

The diagnosis input:

```sh
tools/fuzz_corpus_encode.py --out tests/fuzz/corpus/fuzz_care/definiteness_rank_one_diagnosis_pose.bin \
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

Its solution's four entries are all `0.031622776623…`. The binary128 reference converges in two
steps, so this pose carries an accuracy verdict as well as the definiteness check.

The campaign input:

```sh
tools/fuzz_corpus_encode.py --out tests/fuzz/corpus/fuzz_care/definiteness_rank_one_campaign_pose.bin \
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

Here the binary128 reference does not converge, so the target reaches no accuracy verdict and
whether this answer is accurate is not settled by anything in this tree.

The accuracy-gate input, the one of the four that does not reach the definiteness check:

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
visibility through the weight factor is `|Q(0,0) - Q(0,1)| = 3.3061954049506959e-08`. That weak
visibility is what puts the stabilizing closed loop next to the imaginary axis and the conditioning
past the entitlement.

The second definiteness input, the same class as the diagnosis and campaign inputs:

```sh
tools/fuzz_corpus_encode.py --out tests/fuzz/corpus/fuzz_care/definiteness_rank_one_seedless_pose.bin \
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

Its solution's largest entry is `0.06336703293626454`. The binary128 reference does not converge
here either.

The inputs these poses were recorded from are 88, 104, 102 and 95 bytes. The decoder copies exactly
88 and ignores anything past them, so the trailing bytes of the longer ones carry nothing and each
88-byte encoding above decodes to the same pose; a mutation length is not a field.

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
