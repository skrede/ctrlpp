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

Separately from the cases, `validation/oracles/linear_kalman/` holds a pair of
mutually independent Octave and Python implementations of the linear Kalman
recursion; they are not a cross-validation case but the external oracle for the
golden state an EKF unit test asserts against.

Run the suite:

```bash
cd validation && CTRLPP_VALIDATE_JOBS=6 ./validate.sh
```

`CTRLPP_VALIDATE_JOBS` sets the build parallelism (default 2),
`CTRLPP_VALIDATE_GENERATOR` sets the CMake generator (default `Unix Makefiles`),
and `CTRLPP_VALIDATE_OCTAVE` selects the interpreter to invoke (default
`octave`).  Every case is registered as a test, so a case whose executable was
never built fails the run rather than being reported as skipped.

A case can be waived by placing a `waiver.cfg` in its directory under
`validation/cases/`, holding `WAIVER_REASON` and `WAIVER_DATE`.  Both keys are
mandatory; a missing or empty one stops the configure and names the case.  A
waived case is registered as a disabled test, so it does not run and does not
fail the suite, but it keeps its row in the table below with its reason
reproduced verbatim.  No case is waived today.

`tools/validation_results.py --gate` fails a run in which any case failed, or
did not run without a waiver, or is a disabled case carrying no waiver file; the
test runner alone lets the last two cost nothing.
`tools/validation_results.py --check docs/validation.md` fails when the table
below no longer matches the run, and prints the difference.

The results below come from the pinned reference environment in continuous
integration rather than from a developer's host: GNU Octave 8.4.0 with Control
4.2.3, Signal 1.4.7, Splines 1.3.5 and Quaternion 2.4.2, against the library
built in Release on the `ubuntu-24.04` runner image.  All four reference
packages are present there, so every registered case runs.

The reference environment itself is described by two files.
`validation/octave-packages.sha256` holds the content hash of each of the four
reference package tarballs, and the continuous-integration workflow verifies
every tarball against it before installing any of them, asserts the
interpreter's version at start, and then runs this suite and the freshness
comparison above.  To reproduce that environment locally, download the four
tarballs at the versions the workflow names, check them with
`sha256sum -c validation/octave-packages.sha256`, and install them with `control`
first, because `signal` declares a dependency on it.

That pinning is not a claim of hermetic reproduction.  The reference *packages*
are pinned by content hash and the interpreter's version is asserted; the
interpreter itself, its linear-algebra library, the compiler and the runner image
are not pinned, so reproducing a run from months ago is not achievable on this
route.

### Results

The table below is generated from the run described above by
`tools/validation_results.py`; it is regenerated with `--write` in a reviewed
commit, so a hand edit to it is a difference the check mode reports rather than
a correction.  "Digits" is the minimum digits of agreement across the signals of
a case, `-log10(max relative error)`, and "max abs error" is the absolute error
of that worst signal.  A row reading `n/a` is a case that failed before it
produced anything to compare.

<!-- results:begin: generated from a harness run; edit the harness, not this region -->
| Component | Octave function | Worst-case digits | Max abs error | Verdict |
|-----------|----------------|:-----------------:|:-------------:|:-------:|
| `dare` | `dare()` | 13.5 | 1.49e-12 | PASS |
| `care` | `care()` | 14.6 | 4.88e-15 | PASS |
| `lqr` (infinite horizon) | `dlqr()` | 10.6 | 7.88e-15 | PASS |
| `lqr` (finite horizon) | backward Riccati | 12.6 | 9.10e-15 | PASS |
| `lqi` | `dlqr()` augmented | 2.8 | 9.73e-14 | PASS |
| `pid` (linear PI) | `lsim()` | 10.8 | 1.05e-14 | PASS |
| `place` | `place()` | 13.7 | 2.21e-15 | PASS |
| `discretize` (ZOH) | `c2d()` | 15.6 | 1.39e-17 | PASS |
| `analysis` (poles) | `pole()`, `ctrb()`, `obsv()` | 15.4 | 2.08e-17 | PASS |
| `tf2ss` / `ss2tf` | `tf2ss()`, `ss2tf()` | 15.0 | 1.11e-15 | PASS |
| `kalman_filter` | time-varying KF | 10.6 | 5.55e-15 | PASS |
| `luenberger_observer` | `place()` on dual | 12.9 | 4.72e-15 | PASS |
| `butterworth` (4th order) | `butter()`, `filter()` | 13.2 | 5.94e-14 | PASS |
| `fir` | `filter()` | 15.0 | 1.11e-15 | PASS |
| `cubic_spline` (natural) | `csape()`, `ppder()` | 12.5 | 1.11e-15 | PASS |
| `so3::exp` / `so3::log` | `rot2q()`, `q2rot()` | 15.5 | 1.11e-16 | PASS |
| `batch_arx` | `arx()` | 13.8 | 3.08e-15 | PASS |
| `moesp` (cross-algorithm) | `n4sid()` | 10.7 | 2.74e-13 | PASS |
<!-- results:end -->

Census for that run: all 18 rows carry PASS, none carries FAIL, and none carries
WAIVED.  The three cases that depend on the `signal` and `splines` packages
(`tf2ss` / `ss2tf`, `butterworth`, and `cubic_spline`) run here because the
pinned environment installs those packages; the workflow asserts those three
verdicts by name rather than resting on the suite's exit status.

The `moesp` row is a **cross-algorithm** comparison.  Octave's `n4sid()` is a
different subspace identification algorithm than the routine under test, so the
row shows that two independent subspace methods recover the same input-output
behavior from the same record.  It is not a MOESP-specific reference, and it
does not validate any MOESP internal quantity.

The `batch_arx` row exercises the fallible entry point: `batch_arx` returns
`ctrlpp::expected<arx_result<...>, ctrlpp::sysid_error>`, and the case treats a
rejected record as a case failure rather than emitting an empty comparison file.

### What this page does and does not claim

What the run establishes is that all eighteen cases agree with their Octave
references to the digits shown, in the environment named above.  It does not
establish that those digits reproduce elsewhere: the interpreter, its
linear-algebra library, the compiler and the runner image are all unpinned, and
the freshness check that guards this table is byte-exact, so a run in a
different environment will fail it until the region is regenerated in a
reviewed commit.  The table is the published record of one run, not a standing
claim about every environment.

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
| `discretize` (ZOH) | 2 | Cross-validated against Octave `c2d()` |
| `state_space` | 2 | Implicitly validated via all state-space cross-validation cases |
| `transfer_function` | 2 | Cross-validated via `tf2ss` / `ss2tf` round-trip; same case |
| `analysis` (poles, stability, controllability, observability) | 2 | Cross-validated against Octave `pole()`, `ctrb()`, `obsv()` |
| `conversion` (`tf2ss`, `ss2tf`) | 2 | Cross-validated against Octave `tf2ss()`, `ss2tf()`; the case needs the `signal` package |
| `propagate` | 2 | Implicitly validated via LQR and Kalman time-series cases |
| **MPC / MHE** | | |
| `mpc` | 1 | Optimization-based; requires solver integration |
| `nmpc` | 1 | Optimization-based; requires solver integration |
| `mhe` | 1 | Optimization-based estimation |
| `nmhe` | 1 | Optimization-based estimation |
| **Signal Processing** | | |
| `butterworth` (cascaded biquad) | 2 | Cross-validated against Octave `butter()` + `filter()`; the case needs the `signal` package |
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
| `cubic_spline` (natural) | 2 | Cross-validated against Octave `csape()` + `ppder()`; the case needs the `splines` package |
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
