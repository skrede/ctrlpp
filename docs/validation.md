# Validation Status

How each ctrlpp component is tested and what level of independent validation
it has received.

## Validation levels

| Level | Label | Meaning |
|-------|-------|---------|
| 1 | Unit tested | Catch2 unit tests, property tests, or fuzz targets |
| 2 | Cross-validated | Output compared against an independent reference implementation (GNU Octave toolbox) |
| 3 | Simulated system | Validated in a closed-loop simulation with realistic plant dynamics |
| 4 | Physical system | Validated on real hardware |

Higher levels imply all lower levels.  Every component in ctrlpp has at least
level 1 coverage.

## Cross-validation against GNU Octave

The `validation/` directory registers 18 test cases that run identical
scenarios in both ctrlpp (C++20) and GNU Octave, then compare outputs
column-by-column.  The comparison uses a combined tolerance criterion: pass iff
`|ref - cand| <= atol + rtol * |ref|` for every element (default
`atol = 1e-10`, `rtol = 1e-8`).

The reference scripts collectively require four Octave packages: **control**,
**signal**, **splines**, and **quaternion**.  `control` is loaded by 14 cases,
`signal` by 2, and `splines` and `quaternion` by 1 each; one case
(`fir_filter`) loads no package at all.

`validation/cases/` holds one directory beyond the 18 registered cases:
`ekf_numerical_jacobian` is not a cross-validation case but a pair of mutually
independent Octave and Python oracles for the linear Kalman recursion exercised
by an EKF unit test, so it has neither a reference script named after the case
nor a build file, and the harness skips it.

Run the suite:

```bash
cd validation && CTRLPP_VALIDATE_JOBS=6 ./validate.sh
```

`CTRLPP_VALIDATE_JOBS` sets the build parallelism (default 2) and
`CTRLPP_VALIDATE_GENERATOR` sets the CMake generator (default `Unix Makefiles`).

**The reference environment is not pinned.** The results below were produced by
a local run on GNU Octave 11.3.0 with Control 4.2.2 and Quaternion 2.4.2; the
`signal` and `splines` packages were not installed on that host, so the three
cases that need them did not run and carry no verdict.  Hermetic provisioning of
a fixed reference environment, and promotion of this suite into continuous
integration, are tracked separately and are not claimed here.

### Results

Regenerated from the run described above.  "Digits" is the minimum digits of
agreement across the signals of a case, `-log10(max relative error)`, and "max
abs error" is the absolute error of that worst signal.

| Component | Octave function | Worst-case digits | Max abs error | Verdict |
|-----------|----------------|:-----------------:|:-------------:|:-------:|
| `dare` | `dare()` | 12.9 | 5.66e-12 | PASS |
| `care` | `care()` | 14.7 | 4.00e-15 | PASS |
| `lqr` (infinite horizon) | `dlqr()` | 10.0 | 3.51e-14 | PASS |
| `lqr` (finite horizon) | backward Riccati | 12.5 | 7.99e-15 | PASS |
| `lqi` | `dlqr()` augmented | 2.8 | 6.39e-12 | PASS |
| `pid` (linear PI) | `lsim()` | 10.8 | 1.05e-14 | PASS |
| `place` | `place()` | 13.9 | 1.14e-15 | PASS |
| `discretise` (ZOH) | `c2d()` | 15.6 | 1.39e-17 | PASS |
| `analysis` (poles) | `pole()`, `ctrb()`, `obsv()` | 15.4 | 2.08e-17 | PASS |
| `tf2ss` / `ss2tf` | `tf2ss()`, `ss2tf()` | not run | not run | NOT RUN, needs `signal` |
| `kalman_filter` | time-varying KF | 10.4 | 2.08e-14 | PASS |
| `luenberger_observer` | `place()` on dual | 13.1 | 4.61e-15 | PASS |
| `butterworth` (4th order) | `butter()`, `filter()` | not run | not run | NOT RUN, needs `signal` |
| `fir` | `filter()` | 15.0 | 1.11e-15 | PASS |
| `cubic_spline` (natural) | `csape()`, `ppder()` | not run | not run | NOT RUN, needs `splines` |
| `so3::exp` / `so3::log` | `rot2q()`, `q2rot()` | 15.6 | 1.11e-16 | PASS |
| `batch_arx` | `arx()` | 13.5 | 9.66e-15 | PASS |
| `moesp` (cross-algorithm) | `n4sid()` | 10.7 | 2.59e-13 | PASS |

Census for that run: 15 cases compared and passed, 0 compared and failed, 3 not
run because their packages were absent, and 1 directory skipped because it is
not a case.  The harness exits nonzero because it counts the three unrunnable
cases as failures; that exit status reflects the missing packages, not a
numerical disagreement.

The `moesp` row is a **cross-algorithm** comparison.  Octave's `n4sid()` is a
different subspace identification algorithm than the routine under test, so the
row shows that two independent subspace methods recover the same input-output
behavior from the same record.  It is not a MOESP-specific reference, and it
does not validate any MOESP internal quantity.

The `batch_arx` row exercises the fallible entry point: `batch_arx` returns
`ctrlpp::expected<arx_result<...>, ctrlpp::sysid_error>`, and the case treats a
rejected record as a case failure rather than emitting an empty comparison file.

### What this page does and does not claim

These numbers come from one local run on one host with the package set named
above.  They are not from a pinned environment and not from continuous
integration, so a reader cannot assume they reproduce byte-for-byte elsewhere.
Three rows have no verdict at all because the packages they need were absent on
that host; they are recorded as not run rather than carried forward from an
earlier run.  Making the reference environment reproducible and running this
suite automatically are separate pieces of work.

## Full component matrix

| Component | Level | Notes |
|-----------|:-----:|-------|
| **Control** | | |
| `pid` (linear core) | 2 | Cross-validated against Octave `lsim()` (backward Euler PI) |
| `pid` policies (saturation, anti-windup, rate limiting, gain scheduling) | 1 | Nonlinear; Octave `pid()` is linear-only |
| `pid_performance` | 1 | IAE/ISE/ITAE metrics |
| `lqr` (infinite horizon) | 2 | Cross-validated against Octave `dlqr()` |
| `lqr` (finite horizon) | 2 | Cross-validated against Octave backward Riccati |
| `lqi` | 2 | Cross-validated against Octave `dlqr()` on augmented system |
| `dare` | 2 | Cross-validated against Octave `dare()` |
| `care` | 2 | Cross-validated against Octave `care()` |
| `place` | 2 | Cross-validated against Octave `place()` |
| `mrac` | 1 | No MRAC in any Octave package |
| `l1` | 1 | No L1 adaptive in any Octave package |
| **Estimation** | | |
| `kalman_filter` | 2 | Cross-validated against time-varying Kalman equations |
| `luenberger_observer` | 2 | Cross-validated against Octave `place()` on dual system |
| `ekf` | 1 | Nonlinear; no Octave equivalent |
| `ukf` | 1 | Nonlinear; no Octave equivalent |
| `particle_filter` | 1 | Stochastic; no Octave equivalent |
| `mekf` | 1 | SO(3) error-state; no Octave equivalent. Monte Carlo NEES validated |
| `manifold_ukf` | 1 | SO(3) manifold; no Octave equivalent. Monte Carlo NEES validated |
| `complementary_filter` | 1 | Mahony filter; no Octave equivalent |
| **Model** | | |
| `discretise` (ZOH) | 2 | Cross-validated against Octave `c2d()` |
| `state_space` | 2 | Implicitly validated via all state-space cross-validation cases |
| `transfer_function` | 2 | Cross-validated via `tf2ss` / `ss2tf` round-trip; same case, same `signal` package requirement |
| `analysis` (poles, stability, controllability, observability) | 2 | Cross-validated against Octave `pole()`, `ctrb()`, `obsv()` |
| `conversion` (`tf2ss`, `ss2tf`) | 2 | Cross-validated against Octave `tf2ss()`, `ss2tf()`; not reproduced in the reference run, the case needs the `signal` package |
| `propagate` | 2 | Implicitly validated via LQR and Kalman time-series cases |
| **MPC / MHE** | | |
| `mpc` | 1 | Optimization-based; requires solver integration |
| `nmpc` | 1 | Optimization-based; requires solver integration |
| `mhe` | 1 | Optimization-based estimation |
| `nmhe` | 1 | Optimization-based estimation |
| **Signal Processing** | | |
| `butterworth` (cascaded biquad) | 2 | Cross-validated against Octave `butter()` + `filter()`; not reproduced in the reference run, the case needs the `signal` package |
| `fir` | 2 | Cross-validated against Octave `filter()` |
| `biquad` (low-pass, notch, dirty derivative) | 1 | RBJ cookbook coefficients; no direct Octave equivalent |
| `vector_biquad` | 1 | Multi-channel wrapper over biquad |
| **System Identification** | | |
| `batch_arx` | 2 | Cross-validated against Octave `arx()`; returns `expected<arx_result<...>, sysid_error>` and the case fails on a rejected record |
| `moesp` | 2 | Cross-algorithm comparison against Octave `n4sid()`, a different subspace method; agreement is on input-output behavior, not on MOESP internals |
| `rls` | 1 | No RLS in any Octave package |
| `recursive_arx` | 1 | No recursive ARX in any Octave package |
| **Lie Groups** | | |
| `so3` (exp, log, quaternion) | 2 | Cross-validated against Octave quaternion package |
| **Trajectory** | | |
| `cubic_spline` (natural) | 2 | Cross-validated against Octave `csape()` + `ppder()`; not reproduced in the reference run, the case needs the `splines` package |
| `cubic_path` / `cubic_trajectory` | 1 | Polynomial evaluation; unit tested |
| `quintic_path` / `quintic_trajectory` | 1 | Polynomial evaluation; unit tested |
| `septic_path` / `septic_trajectory` | 1 | Polynomial evaluation; unit tested |
| `harmonic_path` / `cycloidal_path` | 1 | Trigonometric primitives; unit tested |
| `trapezoidal_trajectory` | 1 | Velocity profile; no Octave equivalent |
| `double_s_trajectory` | 1 | 7-segment S-curve; no Octave equivalent |
| `modified_sin_trajectory` / `modified_trap_trajectory` | 1 | Specialized profiles; no Octave equivalent |
| `bspline_trajectory` | 1 | B-spline evaluation; unit tested |
| `smoothing_spline` | 1 | Regularized spline; unit tested |
| `online_planner_2nd` / `online_planner_3rd` | 1 | Real-time planners; unit tested |
| `rescale_to` / `synchronize` | 1 | Trajectory time scaling; unit tested. `rescale_to` returns `expected<void, trajectory_error>` and `synchronize` returns the same over a variadic pack or a `std::span<Profile>`, rejecting all-or-nothing with nothing mutated |
