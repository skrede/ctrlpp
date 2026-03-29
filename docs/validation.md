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

The `validation/` directory contains 17 test cases that run identical
scenarios in both ctrlpp (C++23) and GNU Octave (with the Control and Signal
packages), then compare outputs column-by-column.  The comparison uses a
combined tolerance criterion: pass iff `|ref - cand| <= atol + rtol * |ref|`
for every element (default `atol = 1e-10`, `rtol = 1e-8`).

Reference environment: GNU Octave 11.1.0, Control 4.2.1, Signal 1.4.7.

Run the suite:

```bash
cd validation && ./validate.sh
```

### Results

| Component | Octave function | Worst-case digits | Max abs error | Verdict |
|-----------|----------------|:-----------------:|:-------------:|:-------:|
| `dare` | `dare()` | 14 | 1.26e-12 | PASS |
| `lqr` (infinite horizon) | `dlqr()` | 10 | 9.33e-15 | PASS |
| `lqr` (finite horizon) | backward Riccati | 12 | 7.99e-15 | PASS |
| `lqi` | `dlqr()` augmented | 3 | 7.39e-13 | PASS |
| `pid` (linear PI) | `lsim()` | 11 | 1.05e-14 | PASS |
| `place` | `place()` | 14 | 1.14e-15 | PASS |
| `discretise` (ZOH) | `c2d()` | 16 | 1.39e-17 | PASS |
| `analysis` (poles) | `pole()`, `ctrb()`, `obsv()` | 15 | 2.08e-17 | PASS |
| `tf2ss` / `ss2tf` | `tf2ss()`, `ss2tf()` | 15 | 1.11e-15 | PASS |
| `kalman_filter` | time-varying KF | 11 | 5.56e-15 | PASS |
| `luenberger_observer` | `place()` on dual | 13 | 4.61e-15 | PASS |
| `butterworth` (4th order) | `butter()`, `filter()` | 13 | 4.00e-14 | PASS |
| `fir` | `filter()` | 15 | 1.11e-15 | PASS |
| `cubic_spline` (natural) | `csape()`, `ppder()` | 12 | 1.11e-15 | PASS |
| `so3::exp` / `so3::log` | `rot2q()`, `q2rot()` | 16 | 1.11e-16 | PASS |
| `batch_arx` | `arx()` | 14 | 9.77e-15 | PASS |
| `n4sid` | `n4sid()` | 11 | 2.61e-13 | PASS |

"Digits" is the minimum digits of agreement across all signals in a case:
`-log10(max relative error)`.

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
| `transfer_function` | 2 | Cross-validated via `tf2ss` / `ss2tf` round-trip |
| `analysis` (poles, stability, controllability, observability) | 2 | Cross-validated against Octave `pole()`, `ctrb()`, `obsv()` |
| `conversion` (`tf2ss`, `ss2tf`) | 2 | Cross-validated against Octave `tf2ss()`, `ss2tf()` |
| `propagate` | 2 | Implicitly validated via LQR and Kalman time-series cases |
| **MPC / MHE** | | |
| `mpc` | 1 | Optimization-based; requires solver integration |
| `nmpc` | 1 | Optimization-based; requires solver integration |
| `mhe` | 1 | Optimization-based estimation |
| `nmhe` | 1 | Optimization-based estimation |
| **Signal Processing** | | |
| `butterworth` (cascaded biquad) | 2 | Cross-validated against Octave `butter()` + `filter()` |
| `fir` | 2 | Cross-validated against Octave `filter()` |
| `biquad` (low-pass, notch, dirty derivative) | 1 | RBJ cookbook coefficients; no direct Octave equivalent |
| `vector_biquad` | 1 | Multi-channel wrapper over biquad |
| **System Identification** | | |
| `batch_arx` | 2 | Cross-validated against Octave `arx()` |
| `n4sid` | 2 | Cross-validated against Octave `n4sid()` |
| `rls` | 1 | No RLS in any Octave package |
| `recursive_arx` | 1 | No recursive ARX in any Octave package |
| **Lie Groups** | | |
| `so3` (exp, log, quaternion) | 2 | Cross-validated against Octave quaternion package |
| **Trajectory** | | |
| `cubic_spline` (natural) | 2 | Cross-validated against Octave `csape()` + `ppder()` |
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
| `time_scaling` / `synchronize` | 1 | Trajectory utilities; unit tested |
