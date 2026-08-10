# Benchmarks

Performance measurements for ctrlpp components and comparison comparisons
against other C++ control libraries.

## Running benchmarks

The benchmark suite is a standalone CMake project under `benchmarks/`.

### Prerequisites

- C++20 compiler
- Eigen 3.4+
- OSQP 1.0+ and NLopt (for MPC/MHE benchmarks)
- nanobench (fetched automatically via FetchContent)

### Internal benchmarks

```bash
cd benchmarks
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./bench.sh
```

### Competitive benchmarks

Each competitor is gated behind its own CMake option, all OFF by default.
Enable the ones you have installed:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCTRLPP_BENCH_RUCKIG=ON \
  -DCTRLPP_BENCH_OSQP_EIGEN=ON \
  -DCTRLPP_BENCH_LIBMPC=ON \
  -DCTRLPP_BENCH_CT=ON \
  -DCTRLPP_BENCH_HPIPM=ON \
  -DCTRLPP_BENCH_DRAKE=ON
cmake --build build -j$(nproc)
./bench.sh --comparison
```

| Option | Library | Install (Arch) | Notes |
|--------|---------|----------------|-------|
| `CTRLPP_BENCH_RUCKIG` | ruckig | `yay -S ruckig` | Online trajectory generation |
| `CTRLPP_BENCH_OSQP_EIGEN` | osqp-eigen | `yay -S osqp-eigen` | Eigen wrapper for OSQP QP solver |
| `CTRLPP_BENCH_LIBMPC` | libmpc++ | FetchContent (automatic) | Linear MPC library |
| `CTRLPP_BENCH_CT` | ETH Control Toolbox | `yay -S control-toolbox-optcon` | LQR, iLQR, MPC |
| `CTRLPP_BENCH_HPIPM` | HPIPM | `yay -S hpipm-git blasfeo-git` | High-performance interior-point QP |
| `CTRLPP_BENCH_DRAKE` | Drake | Manual (Bazel) | `find_package` only, no FetchContent |

## Internal benchmark results

> **Stale environment, pending remeasurement.** Every figure in the table below
> was taken on a compiler version that is no longer installed on the measuring
> station, against a linear-algebra dependency recorded as a version range
> rather than a version, and with processor frequency scaling left enabled, so
> the clock was not held at any stated frequency.  None of the three can be
> reconstructed after the fact.  Read the figures as an order of magnitude and
> not as a current measurement.

All measurements taken with nanobench, `Release` build, full perf counters
enabled.  Object construction is outside the benchmark lambda — only the
hot-path call is measured.

| Component | Call | ns/op | ops/sec | Instructions |
|-----------|------|------:|--------:|-------------:|
| `pid` | `compute()` | 3.0 | 328M | 52 |
| `lqr` | `compute()` | 1.1 | 938M | 10 |
| `mrac` | `evaluate()` | 3.8 | 261M | 32 |
| `l1` | `evaluate()` | 9.5 | 106M | 42 |
| `biquad` | `process()` | 5.3 | 190M | 19 |
| `fir` (5-tap) | `process()` | 32.4 | 31M | 152 |
| `kalman_filter` | `predict()` | 17.4 | 58M | 320 |
| `kalman_filter` | `update()` | 297.2 | 3.4M | 3385 |
| `ekf` | `predict()` | 44.7 | 22M | 538 |
| `ekf` | `update()` | 279.7 | 3.6M | 3003 |
| `ukf` | `predict()` | 109.3 | 9.2M | 1412 |
| `ukf` | `update()` | 172.9 | 5.8M | 2360 |
| `so3::exp` | exp map | 17.9 | 56M | 192 |
| `so3::log` | log map | 16.8 | 60M | 191 |
| `cubic_spline` | `evaluate()` | 3.7 | 272M | 88 |
| `online_planner_3rd` | `sample()` | 2.8 | 352M | 31 |
| `dare` | solve | 17.8 us | 56K | 191K |
| `lqr_gain` | DARE + gain | 3.8 us | 266K | 36.5K |
| `batch_arx` | `identify()` | 3.5 us | 283K | 41.2K |
| `moesp` | `identify()` | 578 us | 1.7K | 6.9M |

System dimensions: NX=2 for estimators and sysid, NX=4/NU=2 for LQR/DARE,
5-tap for FIR, 5-knot natural cubic spline.

## Discrete Riccati: the acceptance check and the equilibrated path

The figures in this section were taken on their own machine settings and their
own toolchain, both stated below, and are not read across to the table above,
which was taken on the compiler and optimization level named under
[Environment](#environment).

**What each figure is per.** Every number below is one `ctrlpp::dare` call at the
stated `NX` / `NU`, on the discrete damped chain the benchmark builds -- a
forward-Euler step of the same continuous chain the CARE bakeoff sweeps. The
acceptance-check and forward-error rows are one call of that component alone, on
operands the solve already produced, formed outside the measured region so the
row is the component's own cost and not the cost of preparing its inputs.

**Machine.** AMD Ryzen 7 5800X3D, governor `performance`, boost on, pinned to one
logical CPU with `taskset`. The measurement held a CPU-exclusive station window
and nothing else ran in it. The harness reported a frequency range of 1,755 to
4,553 MHz, which is the boost range rather than an idle-state range, so the clock
was not held at one stated frequency; the run-to-run control at the end of this
section is reported for that reason rather than assumed away.

**Build.** `g++ (GNU) 16.1.1 20260728`, `-O3 -DNDEBUG -std=c++23`, Eigen 3.4.1
(the system release under `/usr/include/eigen3`, resolved by `find_package`, and
not the 3.4.0 the test tree fetches -- the two differ on this solver's stack
figures elsewhere), nanobench 4.3.11, scalar `double`, x86-64 Linux.

**Dispersion.** The harness runs 10 warmup iterations and 51 epochs of at least
1 ms each. 51 is odd, so every reported time is an observed median rather than an
interpolation, and the `+-` beside it is the median absolute percent error over
those same 51 epochs. It is a within-run dispersion and it is not a confidence
interval.

### At a weight scale of one

`ctrlpp::dare` equilibrates before it solves. When the caller's weights already
have a scale of exactly one, equilibration is inactive and the entry point runs
its acceptance check once.

| `NX` / `NU` | `dare` whole entry point | acceptance check alone | forward-error estimate alone |
|---|---:|---:|---:|
| 2 / 1 | 2,328.6 ns +-0.3% | 464.3 ns +-0.4% | 304.0 ns +-0.7% |
| 4 / 2 | 11,231 ns +-0.2% | 3,343.4 ns +-0.3% | 1,788.4 ns +-0.7% |
| 6 / 2 | 33,864 ns +-0.1% | 8,115.1 ns +-0.2% | 5,205.8 ns +-0.4% |
| 8 / 2 | 66,259 ns +-0.3% | 21,319 ns +-0.1% | 16,244 ns +-0.3% |
| 12 / 3 | 228,222 ns +-0.2% | 101,327 ns +-0.2% | 88,239 ns +-0.2% |
| 15 / 3 | 470,065 ns +-0.3% | 278,215 ns +-0.3% | 257,283 ns +-0.4% |

`NX = 15` is a compile-time ceiling and not a choice. The forward-error estimate
holds an `M x M` operator with `M = NX(NX+1)/2` as a fixed-size Eigen object, so
Eigen's 128 KiB stack-allocation limit admits `M <= 128`, that is `NX <= 15`. At
`NX = 16` the acceptance check cannot be instantiated at all.

### The equilibrated path, swept over the weight scale

When the weight scale is anything other than one the entry point takes a longer
route: the acceptance check at the equilibrated scale, a rescale of the answer, a
definiteness factorization of the rescaled answer, the acceptance check **again**
at the caller's scale, and a gain agreement between the two scales. Both weights
are multiplied by `s` rather than `Q` alone, so `Q/s` is the identity exactly and
the Schur solve underneath is the same work at every rung; scaling `Q` alone would
move the weight ratio with the scale and confound this path's cost with a change
in the solve.

The scale was swept over a geometric ladder of **exactly 7 rungs** -- three
decades below one, three above, and one rung immediately beside it. The near-one
rung is `1 + 2^-16`, which is exactly representable, so its weight scale is that
value and not a rounding of it; it separates "the gate is open" from "the pose is
far from equilibrated".

**The equilibrated path's cost is the difference between the two whole-entry-point
rows at the same size.** It is not timed separately and it is not asserted
anywhere. Every difference below is the scaled row in that same row of this table
minus the whole-entry-point row for that `NX` in the table above, so a reader
recomputes all 42 of them from published figures.

| `NX` / `NU` | weight scale `s` | whole entry point at `s` | equilibrated path = difference | as a share of the scale-one row |
|---|---|---:|---:|---:|
| 2 / 1 | `1e-06` | 2,773.6 ns +-0.3% | 445.0 ns | 19.1% |
| 2 / 1 | `1e-04` | 2,777.9 ns +-0.2% | 449.4 ns | 19.3% |
| 2 / 1 | `1e-02` | 2,776.5 ns +-0.2% | 447.9 ns | 19.2% |
| 2 / 1 | `1 + 2^-16` | 2,778.3 ns +-0.2% | 449.7 ns | 19.3% |
| 2 / 1 | `1e+02` | 2,778.5 ns +-0.4% | 449.9 ns | 19.3% |
| 2 / 1 | `1e+04` | 2,777.3 ns +-0.2% | 448.8 ns | 19.3% |
| 2 / 1 | `1e+06` | 2,776.8 ns +-0.2% | 448.3 ns | 19.3% |
| 4 / 2 | `1e-06` | 14,658 ns +-0.1% | 3,427.1 ns | 30.5% |
| 4 / 2 | `1e-04` | 14,654 ns +-0.2% | 3,423.5 ns | 30.5% |
| 4 / 2 | `1e-02` | 14,699 ns +-0.4% | 3,468.4 ns | 30.9% |
| 4 / 2 | `1 + 2^-16` | 14,664 ns +-0.2% | 3,432.6 ns | 30.6% |
| 4 / 2 | `1e+02` | 14,668 ns +-0.2% | 3,436.6 ns | 30.6% |
| 4 / 2 | `1e+04` | 14,678 ns +-0.2% | 3,446.8 ns | 30.7% |
| 4 / 2 | `1e+06` | 14,645 ns +-0.2% | 3,414.3 ns | 30.4% |
| 6 / 2 | `1e-06` | 42,303 ns +-0.2% | 8,438.6 ns | 24.9% |
| 6 / 2 | `1e-04` | 42,281 ns +-0.1% | 8,416.9 ns | 24.9% |
| 6 / 2 | `1e-02` | 42,274 ns +-0.1% | 8,409.6 ns | 24.8% |
| 6 / 2 | `1 + 2^-16` | 42,312 ns +-0.2% | 8,447.3 ns | 24.9% |
| 6 / 2 | `1e+02` | 42,282 ns +-0.1% | 8,418.3 ns | 24.9% |
| 6 / 2 | `1e+04` | 42,312 ns +-0.1% | 8,448.2 ns | 24.9% |
| 6 / 2 | `1e+06` | 42,288 ns +-0.1% | 8,424.1 ns | 24.9% |
| 8 / 2 | `1e-06` | 88,065 ns +-0.2% | 21,807 ns | 32.9% |
| 8 / 2 | `1e-04` | 88,038 ns +-0.3% | 21,780 ns | 32.9% |
| 8 / 2 | `1e-02` | 88,118 ns +-0.2% | 21,859 ns | 33.0% |
| 8 / 2 | `1 + 2^-16` | 88,032 ns +-0.3% | 21,773 ns | 32.9% |
| 8 / 2 | `1e+02` | 88,122 ns +-0.3% | 21,863 ns | 33.0% |
| 8 / 2 | `1e+04` | 87,946 ns +-0.3% | 21,687 ns | 32.7% |
| 8 / 2 | `1e+06` | 88,008 ns +-0.3% | 21,750 ns | 32.8% |
| 12 / 3 | `1e-06` | 332,348 ns +-0.3% | 104,125 ns | 45.6% |
| 12 / 3 | `1e-04` | 331,393 ns +-0.2% | 103,170 ns | 45.2% |
| 12 / 3 | `1e-02` | 331,300 ns +-0.2% | 103,078 ns | 45.2% |
| 12 / 3 | `1 + 2^-16` | 331,608 ns +-0.1% | 103,385 ns | 45.3% |
| 12 / 3 | `1e+02` | 331,997 ns +-0.3% | 103,774 ns | 45.5% |
| 12 / 3 | `1e+04` | 331,123 ns +-0.1% | 102,900 ns | 45.1% |
| 12 / 3 | `1e+06` | 334,263 ns +-0.3% | 106,041 ns | 46.5% |
| 15 / 3 | `1e-06` | 754,410 ns +-0.4% | 284,345 ns | 60.5% |
| 15 / 3 | `1e-04` | 751,205 ns +-0.3% | 281,140 ns | 59.8% |
| 15 / 3 | `1e-02` | 752,670 ns +-0.5% | 282,605 ns | 60.1% |
| 15 / 3 | `1 + 2^-16` | 750,120 ns +-0.3% | 280,054 ns | 59.6% |
| 15 / 3 | `1e+02` | 751,260 ns +-0.3% | 281,195 ns | 59.8% |
| 15 / 3 | `1e+04` | 749,010 ns +-0.2% | 278,945 ns | 59.3% |
| 15 / 3 | `1e+06` | 750,420 ns +-0.2% | 280,355 ns | 59.6% |

### Does the overhead depend on the weight scale?

**No. It is flat in the weight scale over the twelve decades from `1e-06` to
`1e+06`, and it depends on the state dimension instead.** Budget it against `NX`
and ignore how far from equilibrated the pose is.

The retired-instruction counters settle this more sharply than the times can. The
scaled call retires the same instruction count at every rung:

| `NX` | scale-one call | scaled call, all 7 rungs | equilibrated path |
|---|---:|---:|---:|
| 2 | 27,542 | 32,254 | 4,712 |
| 4 | 152,844 | 190,977 | 38,133 |
| 6 | 511,070 | 624,745 | 113,675 |
| 8 | 1,088,895 | 1,415,862 to 1,416,318 | 326,967 to 327,423 |
| 12 | 4,030,252 | 5,797,741 | 1,767,489 |
| 15 | 8,734,242 | 13,923,238 | 5,188,996 |

At five of the six sizes the count is one value across all 7 rungs. Only `NX = 8`
resolves into three distinct values, spanning 456 instructions out of 1,415,862,
which is 0.032% of the row. The equilibrated path executes the same work at
`1e-06` as at `1e+06`.

The per-rung spread in the timed share is 0.11 to 1.45 percentage points across
the six sizes. That is the size of the difference of two rows each carrying 0.1
to 0.7% dispersion, and it does not order itself along the ladder, so it is noise
and not a trend.

What the overhead does depend on is `NX`: 19.3% at `NX = 2`, 30.6% at `NX = 4`,
24.9% at `NX = 6`, 32.9% at `NX = 8`, 45.3% at `NX = 12` and 59.6% at `NX = 15`,
taking the near-one rung as the representative. The rise is not monotone -- `NX = 6`
sits below `NX = 4` in both runs and in the instruction counts, so that dip is a
property of the work and not of the timing.

### Counted operations and measured time are different quantities

The counted-operation figures stand and are not replaced by anything above. The
forward-error estimate assembles its operator in exactly `NX^4` operations and
solves it in `(2/3) M^3` with `M = NX(NX+1)/2`; against a seven-iteration Schur
solve that leading-term ratio is 3.6% at `NX = 2`, 10.7% at `NX = 4`, 24.5% at
`NX = 6` and 47.7% at `NX = 8`. `docs/api/control/dare.md` carries the derivation.

Measured on this machine, the same ratio -- the forward-error row over the solve
with no acceptance arithmetic in it, which is the whole entry point minus the
acceptance check -- is 16.3%, 22.7%, 20.2% and 36.1% at those four dimensions.

**The two do not agree, and neither is wrong.** A leading-term operation count and
a wall time answer different questions: the count is what an implementation is
checked against and carries to any machine, while the time includes everything
this machine does and carries to no other. The measured share exceeds the counted
one by more than four times at `NX = 2` and falls below it at `NX = 6` and
`NX = 8`, so the counted ratio is not a wall-time predictor in either direction.
**Why** they diverge is not established by this measurement; a mechanism that
would fit -- fixed per-call cost dominating the smallest problem, and the
estimator's real-arithmetic kernels retiring more work per cycle than the Schur
solve's complex arithmetic -- is a conjecture and was not measured.

### Run-to-run reproducibility

The set was run twice in the same window on the same binary. The second run is a
reproducibility control: it is reported beside the first, never averaged with it
and never differenced against it. Each run's overhead is a difference of two rows
taken inside that same run, so the two give independently computed answers to the
same question.

| `NX` | overhead band over the 7 rungs, published run | overhead band over the 7 rungs, control run |
|---|---:|---:|
| 2 | 19.1% to 19.3% | 20.2% to 21.7% |
| 4 | 30.4% to 30.9% | 30.5% to 30.8% |
| 6 | 24.8% to 24.9% | 24.8% to 26.2% |
| 8 | 32.7% to 33.0% | 33.3% to 33.6% |
| 12 | 45.1% to 46.5% | 43.4% to 43.6% |
| 15 | 59.3% to 60.5% | 59.4% to 60.6% |

The bands overlap at `NX = 4`, `NX = 6` and `NX = 15` and do not overlap at
`NX = 2`, `NX = 8` and `NX = 12`. **Run-to-run variation in these shares is
therefore wider than the within-run dispersion, and wider than the per-rung
spread the ladder shows.** Read the shares to about two percentage points and no
finer. The instruction counts, which are identical in both runs, carry the
scale-independence conclusion; the times do not have to.

## Competitive comparisons

> **Stale environment, pending remeasurement.** Every figure and every ratio in
> this section was taken on a compiler version that is no longer installed on
> the measuring station, against a linear-algebra dependency recorded as a
> version range rather than a version, and with processor frequency scaling
> left enabled, so the clock was not held at any stated frequency.  None of the
> three can be reconstructed after the fact.  Read the ratios as an order of
> magnitude and not as a current measurement.

nanobench relative mode: ctrlpp is the baseline (100%).  Higher percentage
means the competitor is faster; lower means ctrlpp wins.

### Trajectory generation: ctrlpp vs ruckig

| Benchmark | ctrlpp | ruckig | Ratio |
|-----------|-------:|-------:|------:|
| 3rd-order online planner | 2.8 ns (32 ins) | 57.6 ns (414 ins) | **ctrlpp 20x faster** |

Both solve the same problem: jerk-limited online trajectory generation with
mid-motion retargeting.  ctrlpp's `online_planner_3rd` evaluates a
precomputed polynomial, while ruckig's `Ruckig::update` recomputes the full
time-optimal profile each call.

### MPC: ctrlpp vs libmpc++

| Benchmark | ctrlpp | libmpc++ | Ratio |
|-----------|-------:|---------:|------:|
| Linear MPC solve (NX=4, NU=2, N=10) | 62 us (876K ins) | 305 us (4.7M ins) | **ctrlpp 5x faster** |

Both use OSQP as the QP backend.  The difference is in the QP formulation
layer: ctrlpp builds the sparse QP directly from the state-space model,
while libmpc++ has a heavier abstraction with logging and intermediate
allocations.

### QP solve: ctrlpp vs osqp-eigen

| Benchmark | ctrlpp | osqp-eigen | Ratio |
|-----------|-------:|-----------:|------:|
| QP solve (MPC-sized) | 5.8 us (62K ins) | 5.7 us (61K ins) | ~parity |

Both are thin Eigen wrappers around the same OSQP C solver.  The near-parity
confirms that ctrlpp's QP layer adds negligible overhead.

### LQR gain: ctrlpp vs ETH Control Toolbox

| Benchmark | ctrlpp | ct_optcon | Ratio |
|-----------|-------:|----------:|------:|
| LQR gain (NX=4, NU=2) | 9.0 us (102K ins) | 3.1 us (47K ins) | **ct 2.8x faster** |

The two sides do not solve the same equation, and they do not do the same
amount of work around the solve.  What each one does:

- **Arithmetic**: both factor a `2n x 2n` matrix in real arithmetic.  ct
  applies `Eigen::RealSchur` to the Hamiltonian; `ctrlpp::dare` applies
  `Eigen::RealSchur` to the symplectic matrix
  (`lib/ctrlpp/include/ctrlpp/control/dare.h`).
- **Schur reordering**: ct hands its real Schur factor to the system LAPACK
  library's reordering routine, taking that path only when the toolbox is
  built against LAPACK and falling back to a matrix-sign iteration when it is
  not.  ctrlpp reorders in the repository, in
  `lib/ctrlpp/include/ctrlpp/detail/schur_reorder.h`: the Bai and Demmel 1993
  predicated swap kernel over all four adjacent-block configurations, with
  Murnaghan standardization of every `2x2` block a swap touches.  The
  eigenvalues it moves to the leading position are the ones the discrete
  solver's predicate accepts, `|lambda| < 1 - 2n eps ||T||_max`: inside the
  unit disk by more than the Schur factor's own backward error.
- **Problem formulation**: ct solves the continuous-time ARE (CARE) while
  `ctrlpp::lqr_gain` solves the discrete-time ARE (DARE).  The symplectic
  matrix requires an inverse of `A^T` during setup that the Hamiltonian does
  not; both sides invert `R`.
- **Acceptance**: `ctrlpp::dare` verifies the matrix it is about to return
  before returning it -- a definiteness factorization, a forward-error
  estimate against a half-significand criterion, and a closed-loop eigenvalue
  check -- and repeats that verification at the caller's weight scale whenever
  the scale is not one.  ct returns the extracted solution unverified.  That
  verification is a measurable share of the discrete entry point's cost;
  [Discrete Riccati: the acceptance check and the equilibrated
  path](#discrete-riccati-the-acceptance-check-and-the-equilibrated-path)
  measures it on its own.

## Environment

Results above were collected on:

- CPU: AMD (Zen-class), processor frequency scaling left enabled.  The clock
  was not held at any stated frequency and no governor, boost range or core
  pinning was recorded, so these figures cannot be reproduced even on this
  machine.
- OS: EndeavourOS (Arch Linux), kernel 6.18
- Compiler: GCC 15.2.1, `-O2` (Release)
- Eigen: not recorded.  `3.4+` names a release family, not the build these
  figures were measured against.
- OSQP: 1.0.0
- nanobench: 4.3.11

For stable results, disable CPU frequency scaling before benchmarking:

```bash
sudo cpupower frequency-set -g performance
```

Or use `pyperf system tune` as nanobench suggests.
