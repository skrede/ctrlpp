# Curated fuzz seed corpus

A small, hand-authored regression set for the trajectory fuzz targets. One directory per target,
named to match the target binary exactly, which is libFuzzer's own corpus-directory convention:
the directory is handed to the binary with no translation layer.

```
tests/fuzz/corpus/
  fuzz_trapezoidal/
  fuzz_double_s/
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

## Provenance

| target | seed file | what it reproduces | source revision | decoder contract assumed |
|---|---|---|---|---|
| `fuzz_trapezoidal` | `plateau_cruise_above_both.bin` | retiming that lands on a cruise velocity at or above both boundary velocities | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `ramp_through_cruise_between.bin` | retiming that lands on a cruise velocity strictly between the two boundary velocities, where the branch is a hyperbola and not a quadratic | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `valley_cruise_below_both.bin` | retiming that lands on a cruise velocity at or below both boundary velocities, so the axis dips and climbs back | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `valley_equal_boundary_velocities.bin` | the recorded counterexample: equal boundary velocities where the plateau branch's smaller root falls outside its own validity interval and produces a negative acceleration phase, and the valley root is the admissible one | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `rejection_displacement_below_kinetic.bin` | displacement below the distance the two boundary velocities already sweep, so the reachable durations have a finite supremum and a request past it must be a typed rejection that leaves the profile untouched | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_trapezoidal` | `zero_displacement_nonzero_boundary.bin` | nothing to traverse and opposing boundary velocities, asked for a longer duration: a typed rejection, not a duration the profile does not realize | `b03a1a2` | 7 little-endian binary64: `q0, q1, v_max, a_max, v0, v1, stretch`; minimum 56 bytes |
| `fuzz_double_s` | `no_cruise_unequal_boundary.bin` | the recorded counterexample, construction only: a shape with no cruise segment and unequal boundary velocities, the family whose earlier back-off loop could not terminate | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `no_cruise_unequal_boundary_stretched.bin` | the same configuration retimed halfway to its reachable supremum, so the seed also carries the time-scaling leg | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `no_cruise_negative_boundary.bin` | the recorded counterexample whose position was wrong by three orders of magnitude, construction only: no cruise segment, both boundary velocities negative against a negative displacement | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `no_cruise_negative_boundary_stretched.bin` | the same configuration retimed halfway to its reachable supremum | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `rest_to_rest_closed_form.bin` | rest-to-rest retiming, the closed-form path where the scale is a single quotient of two durations | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `equal_boundary_velocities.bin` | equal nonzero boundary velocities, the family the earlier construction happened to survive and which must therefore keep passing | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |
| `fuzz_double_s` | `rejection_past_reachable_supremum.bin` | a request at twice the reachable supremum: a typed rejection that leaves the profile untouched | `b03a1a2` | 8 little-endian binary64: `q0, q1, v_max, a_max, j_max, v0, v1, stretch`; minimum 64 bytes |

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

A `stretch` of one is below the target's own `T_new > T` guard, so those two seeds exercise
construction and the dense scan without entering the time-scaling leg. That is deliberate: it keeps
the construction of each recorded counterexample pinned independently of the retiming.

## Adding a seed

Author the bytes from named configuration values with a throwaway encoder rather than assembling
them by hand: a hand-assembled seed is unverifiable and undebuggable. Check the emitted values
against the target's own clamps so the decoded configuration is the one intended, keep the file at or
above the target's minimum decoded size and below the smoke run's length cap, give it a descriptive
name so a failing replay reads clearly, and add its row above.
