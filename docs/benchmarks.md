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

Taken with nanobench in a `Release` build with performance counters enabled,
under the conditions the [Environment](#environment) section states in full:
one processor core, pinned, on an otherwise idle machine, three serial passes,
medians reported.  Object construction is outside the benchmark lambda, so only
the hot-path call is measured.

**Instructions are the archival-grade column.**  They are frequency invariant,
and they reproduced to 0.00% across all three passes on every row here except
`ukf::update` (1.34%).  Wall time is directional: processor boost is on, so the
clock is not held at a stated frequency and `ns/op` carries boost variance.
Where the two columns disagree about a ratio, trust the instruction count.

The **spread** column is the cross-pass spread on wall time,
`(max - min) / median` over the three passes.

| Component | Call | ns/op | ops/sec | Instructions | Compiler | Spread |
|-----------|------|------:|--------:|-------------:|----------|-------:|
| `pid` | `compute()` | 8.12 | 123M | 162 | GCC 16.1.1 | 1.06% |
| `lqr` | `compute()` | 0.68 | 1477M | 10 | GCC 16.1.1 | 1.00% |
| `mrac` | `evaluate()` | 5.15 | 194M | 77 | GCC 16.1.1 | 0.00% |
| `l1` | `evaluate()` | 10.24 | 98M | 88 | GCC 16.1.1 | 0.00% |
| `biquad` | `process()` | 5.13 | 195M | 19 | GCC 16.1.1 | 0.39% |
| `fir` (5-tap) | `process()` | 30.09 | 33M | 170 | GCC 16.1.1 | 4.17% |
| `kalman_filter` | `predict()` | 17.04 | 59M | 311 | GCC 16.1.1 | 1.81% |
| `kalman_filter` | `update()` | 313.5 | 3.2M | 3572 | GCC 16.1.1 | 2.51% |
| `ekf` | `predict()` | 32.95 | 30M | 520 | GCC 16.1.1 | 1.10% |
| `ekf` | `update()` | 223.6 | 4.5M | 2659 | GCC 16.1.1 | 1.79% |
| `ukf` | `predict()` | 86.28 | 12M | 1183 | GCC 16.1.1 | 0.64% |
| `ukf` | `update()` | 306.3 | 3.3M | 3342 | GCC 16.1.1 | 2.70% |
| `so3::exp` | exp map | 41.27 | 24M | 217 | GCC 16.1.1 | 0.94% |
| `so3::log` | log map | 15.96 | 63M | 190 | GCC 16.1.1 | 0.76% |
| `cubic_spline` | `evaluate()` | 3.80 | 263M | 92 | GCC 16.1.1 | 1.22% |
| `online_planner_3rd` | `sample()` | 8.25 | 121M | 128 | GCC 16.1.1 | 0.23% |
| `dare` | solve + acceptance check | 10.78 us | 92.8K | 153K | GCC 16.1.1 | 2.08% |
| `lqr_gain` | DARE + gain | 3.01 us | 333K | 36.4K | GCC 16.1.1 | 0.34% |
| `batch_arx` | `identify()` | 4.38 us | 228K | 58.1K | GCC 16.1.1 | 0.51% |
| `moesp` | `identify()` | 567 us | 1.8K | 6.88M | GCC 16.1.1 | 0.68% |

System dimensions: NX=2 for estimators and sysid, NX=4/NU=2 for LQR/DARE,
5-tap for FIR, 5-knot natural cubic spline.

### What each row checks about its own answer

A speed number says nothing about whether the answer is right, so every row
either carries a figure that does say so or states plainly that no such figure
exists.  These are computed outside every timed region and are independent of
processor contention.

| Row | What is checked | Value |
|---|---|---:|
| `biquad::process` | relative deviation of the settled constant response from the zero-frequency gain of its own coefficients | 8.3e-15 |
| `fir::process` | max relative deviation of its own impulse response from its tap vector | 0.0 |
| `so3::exp` | quaternion distance between `exp(2w)` and `exp(w)` composed with itself | 0.0 |
| `so3::log` | relative deviation of `log(exp(w))` from the rotation vector it was built from | 1.2e-16 |
| `cubic_spline::evaluate` | max relative deviation from the knot positions it interpolates | 1.1e-16 |
| `online_planner_3rd::sample` | max relative violation of the point symmetry of its own profile about its own midpoint | 8.9e-16 |
| `dare` | relative residual of its own discrete Riccati solution | 2.7e-15 |
| `lqr_gain` | discrete closed-loop spectral radius of its own gain (below one certifies stability) | 0.9914 |
| `batch_arx::identify` | max relative deviation of its Markov parameters from the generating system's | 2.3e-15 |
| `moesp::identify` | same metric | 6.4e-11 |

`pid`, `lqr::compute`, `mrac`, `l1` and the six estimator rows carry **no such
figure, and the omission is deliberate**.  Each times one step of a recursion,
so any residual formed from the same difference equation would re-derive the
step from itself and report zero whether or not the step is correct.  That is a
tautology, not a check.

`lqr_gain` also admits a residual of the Riccati solution reconstructed from its
own gain, measured at 1.4e-17.  **Read that as an upper bound, not as a
measurement of suboptimality.**  The reconstruction is quadratic in the gain
error, so at a gain error near machine precision the true suboptimality term
falls below anything double precision represents, and the figure is the
reconstruction's own numerical floor.  The spectral radius above is the sharper
of the two statements.

### Reproducing this table

```bash
cmake -S benchmarks -B build-bench -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DCTRLPP_FETCH_BENCHMARK_DEPS=ON
cmake --build build-bench -j6 --target bench_pid bench_lqr bench_dare bench_kalman \
  bench_ekf bench_ukf bench_mrac bench_l1 bench_sysid bench_trajectory bench_dsp bench_so3
mkdir -p run && cd run && taskset -c 2 ../build-bench/internal/bench_pid
```

Each benchmark writes a fixed filename into the working directory, so give every
run a directory of its own.  Substitute any other target name on the last line.

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

Taken in the same window as the internal table, under the conditions the
[Environment](#environment) section states in full.  Identical problem data is
handed to both sides of every comparison.  Medians of three serial passes.

Two columns carry the disclosure that a bare ratio cannot:

- **Compiler.**  Every row below was built with GCC 16.1.1 except the nonlinear
  predictive control row, which needs GCC 13.4.1 and says so.  Absolute times
  are not comparable across that boundary.
- **Agreement.**  How far apart the two implementations' answers are.  A speed
  ratio between two arms that did not reach the same answer is not a
  comparison, so this figure is what makes each row legitimate.  Beside it,
  where one is definable, each arm carries a figure for how far its own answer
  is from satisfying the equation it claims to solve.  **Two implementations
  that agree perfectly can be wrong together**, which is why both figures
  appear rather than only the first.

`instr x` is competitor instructions over ctrlpp instructions, so above one
means ctrlpp does less work.  **Instructions are the archival-grade metric**;
wall time is directional because processor boost is on.

### Trajectory generation: ctrlpp vs ruckig

Building a time-optimal profile and reading a state off a built profile are
different operations, and they are timed apart.

| Benchmark | ctrlpp | ruckig | wall x | **instr x** | Compiler | Agreement |
|-----------|-------:|-------:|-------:|------------:|----------|----------:|
| Jerk-limited profile **synthesis** | 24.7 ns (314 ins) | 297.3 ns (3309 ins) | 12.0 | **10.5** | GCC 16.1.1 | 2.2e-16 |
| Built-profile **evaluation** | 8.2 ns (128 ins) | 5.7 ns (101 ins) | 0.70 | **0.79** | GCC 16.1.1 | 8.9e-16 |

Cross-pass spread at or below 4.35% on all four arms.

**ctrlpp is about 10.5x fewer instructions on synthesis; ruckig is about 1.27x
fewer on evaluation.**  The synthesis row races
`ctrlpp::online_planner_3rd::update` against `ruckig::Ruckig::calculate`; the
evaluation row races `online_planner_3rd::sample` against
`ruckig::Trajectory::at_time`.

Both sides produce the **same** time-optimal profile: their total durations
agree to 2.2e-16, one unit in the last place.  Agreement on synthesis is that
duration difference; on evaluation it is the largest relative deviation of the
two arms' sampled position, velocity and acceleration over the common horizon.

Each arm's own answer is checked too.  On synthesis, the absolute position error
of its own profile against the commanded target: ctrlpp 6.7e-16, ruckig 0.  On
evaluation, the point-symmetry residual of its own profile about its own
midpoint: ctrlpp 8.9e-16, ruckig 4.4e-16.  All four are a few units in the last
place of the commanded displacement, and a difference of three or four units in
the last place does not distinguish the two libraries.

Both profiles are certified admissible against the shared velocity,
acceleration and jerk limits.  The worst relative margin, **net of the scan
resolution it was observed at**, is -6.7e-09 for ctrlpp and -4.4e-09 for
ruckig.  Both negative, so neither profile violates a limit by more than the
sampling could explain.

```bash
cmake -S benchmarks -B build-bench -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DCTRLPP_FETCH_BENCHMARK_DEPS=ON -DCTRLPP_BENCH_RUCKIG=ON
cmake --build build-bench -j6 --target bench_trajectory_vs_ruckig
mkdir -p run && cd run && taskset -c 2 ../build-bench/comparison/ruckig/bench_trajectory_vs_ruckig
```

### MPC: ctrlpp vs libmpc++

| Benchmark | ctrlpp | libmpc++ | wall x | **instr x** | Compiler | Agreement |
|-----------|-------:|---------:|-------:|------------:|----------|----------:|
| Linear MPC solve (NX=4, NU=2, N=10) | 57.2 us (863K ins) | 7553 us (127M ins) | 132 | **147** | GCC 16.1.1 | 2.3e-12 |

Cross-pass spread at or below 1.73%.

**ctrlpp is about 147x fewer instructions.**  Both use OSQP as the QP backend.
The difference is in the QP formulation layer: ctrlpp builds the sparse QP
directly from the state-space model, while libmpc++ has a heavier abstraction
with logging and intermediate allocations.

**The two arms are run at matched settings, and the ratio depends on that.**
libmpc++'s parameters are set to ctrlpp's tolerances, to an iteration budget of
4000 rather than its own default of 100, to warm start on rather than its own
default of off, and to the same polish setting; and the terminal weight is
pinned to the stage weight so both arms minimize the same functional, which
libmpc++'s per-step diagonal output weight cannot otherwise express.  With the
weights matched, both libraries report the same objective to nine figures.

**Left at its own defaults libmpc++ is far faster and returns an infeasible
answer.**  Its returned plan then violates the input box by 8.6e-02 against a
bound of 1.0 while reporting itself feasible, and its reported cost is lower
than ctrlpp's precisely because it is optimizing over a set it has not reached.
Each arm's own plan is therefore checked against the plant dynamics and the box
bounds: at matched settings ctrlpp violates by 1.2e-16 and libmpc++ by 1.2e-12.
That column is the only thing on this page that could have caught the defaulted
comparison, since the timings looked reasonable and the competitor's own status
flag said the answer was feasible.

```bash
cmake -S benchmarks -B build-bench -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DCTRLPP_FETCH_BENCHMARK_DEPS=ON -DCTRLPP_BENCH_LIBMPC=ON
cmake --build build-bench -j6 --target bench_mpc_vs_libmpc
mkdir -p run && cd run && taskset -c 2 ../build-bench/comparison/libmpc/bench_mpc_vs_libmpc
```

### QP solve: ctrlpp vs osqp-eigen

| Benchmark | ctrlpp | osqp-eigen | wall x | **instr x** | Compiler | Agreement |
|-----------|-------:|-----------:|-------:|------------:|----------|----------:|
| QP solve (MPC-sized, 2 of 5 rows active) | 8.00 us (90.9K ins) | 7.84 us (89.7K ins) | 0.98 | **0.99** | GCC 16.1.1 | 2.8e-17 |

Cross-pass spread at or below 3.23%.

**Parity: ctrlpp runs 1.3% more instructions.**  Both are thin Eigen wrappers
around the same OSQP 1.0.0 solver, given the same operands and the same
settings, so bit-identical answers are the correct outcome rather than a
defect.  The near-parity confirms that ctrlpp's QP layer adds negligible
overhead, which is the only claim this row supports.

Each arm's own primal-dual answer is checked against the optimality conditions:
relative residual 6.9e-17 on both.  That is the load-bearing figure, because the
agreement column cannot detect a shared error -- two identical wrappers around a
broken solver would agree perfectly.

**The constraints bind.**  At the solution 2 of the 5 rows are active, one at
its upper bound with multiplier +0.5509 and one at its lower bound with
-0.2414, the active rows are full rank, and the condition number of the active
optimality system is 8.67.  This matters because a corpus whose rows are all
slack races quadratic-programming solvers on an effectively unconstrained dense
program, which is not what such a solver is for.  The active set is determined
by exact enumeration rather than by thresholding a slack, so no tolerance enters
the classification, and both quadratic-programming comparisons refuse to run on
a corpus that has gone slack.

```bash
cmake -S benchmarks -B build-bench -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DCTRLPP_FETCH_BENCHMARK_DEPS=ON -DCTRLPP_BENCH_OSQP_EIGEN=ON
cmake --build build-bench -j6 --target bench_qp_vs_osqp_eigen
mkdir -p run && cd run && taskset -c 2 ../build-bench/comparison/osqp_eigen/bench_qp_vs_osqp_eigen
```

**This block was executed verbatim after the measurement window closed**, into a
fresh working directory on the same idle machine.  It returned 7.99 us for
ctrlpp against the 8.00 us published above (0.18% apart, against a stated spread
of 1.61%) and 7.89 us for osqp-eigen against 7.84 us (0.72% apart, against 3.23%).
Both instruction counts and both accuracy figures reproduced exactly.

### Nonlinear predictive control: ctrlpp vs ETH Control Toolbox

| Benchmark | ctrlpp | ct_optcon | wall x | **instr x** | Compiler | Agreement |
|-----------|-------:|----------:|-------:|------------:|----------|----------:|
| Closed-loop nonlinear MPC (damped oscillator, 20 steps) | 28.6 ms (486M ins) | 0.376 ms (4.66M ins) | 0.013 | **0.0096** | **GCC 13.4.1** | 9.1e-08 |

Cross-pass spread at or below 1.98%.

**ct_optcon is about 104x fewer instructions on this problem.**  The two sides
run different algorithm families: ctrlpp solves a sequential quadratic program
over the whole horizon, while ct_optcon runs a Gauss-Newton multiple-shooting
sweep.  The timing is a legitimate comparison only because they reach the same
solution, and they do: the first applied inputs agree to 9.1e-08 and the costs
each arm's own closed loop realizes agree to 1.7e-4 relative, with ctrlpp's
marginally lower (1.04935 against 1.04953) under one objective the benchmark
harness owns rather than either library's internal cost.

**This is the one row on this page built with a different compiler.**  The
competitor's umbrella header assigns to a `const` member at
`ct/optcon/dms/dms_core/TimeGrid.h:53`, which GCC 13 accepts and later releases
reject, so this target gets its own build tree and its own compiler.  Do not
compare its absolute times against any other row here.

Two further disclosures the ratio depends on: the ctrlpp arm is the
**runtime-horizon** controller, not the compile-time-horizon default, and the
competitor's iteration budget is 10.  Controller construction is inside the
timed region on both arms, so a repeated call repeats cold work rather than
warm-starting off the previous one.

```bash
cmake -S benchmarks -B build-bench-gcc13 -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DCTRLPP_FETCH_BENCHMARK_DEPS=ON -DCMAKE_CXX_COMPILER=g++-13 -DCTRLPP_BENCH_CT=ON
cmake --build build-bench-gcc13 -j6 --target bench_nmpc_vs_ct
mkdir -p run && cd run && taskset -c 2 ../build-bench-gcc13/comparison/ct/bench_nmpc_vs_ct
```

### Continuous Riccati and LQR gain: ctrlpp vs ETH Control Toolbox

Both sides solve the **same** continuous algebraic Riccati equation on identical
data, swept over problem size.  Compiler GCC 16.1.1 on every rung.

| NX | ctrlpp `care` | ct `CARE` | **instr x** | ctrlpp own residual | ct own residual |
|---:|--------------:|----------:|------------:|--------------------:|----------------:|
| 2 | 1.42 us (17.0K) | 1.55 us (21.8K) | **1.28** | 1.4e-16 | 1.2e-15 |
| 4 | 6.08 us (82.0K) | 6.36 us (102K) | **1.24** | 1.9e-16 | 1.4e-15 |
| 6 | 14.2 us (221K) | 18.0 us (306K) | **1.39** | 3.3e-16 | 2.3e-15 |
| 8 | 27.2 us (466K) | 33.6 us (637K) | **1.37** | 5.2e-16 | 3.0e-15 |
| 12 | 67.3 us (1.22M) | 80.6 us (1.60M) | **1.31** | 4.9e-16 | 5.3e-15 |
| 16 | 126 us (2.42M) | 153 us (3.22M) | **1.33** | 5.6e-16 | 4.8e-15 |
| 20 | 211 us (4.17M) | 267 us (5.60M) | **1.34** | 4.2e-16 | 5.0e-15 |
| 24 | 321 us (6.48M) | 407 us (8.58M) | **1.32** | 5.0e-16 | 4.9e-15 |
| 30 | 586 us (11.8M) | 762 us (15.8M) | **1.34** | 5.1e-16 | 9.1e-15 |

Cross-pass spread at or below 2.71% on every row.  Agreement between the two
arms' solutions runs 1.9e-15 at NX=2 to 3.0e-14 at NX=30.

**ctrlpp is 1.24x to 1.39x fewer instructions across a fifteen-fold size
sweep, and separately its own residual is flat where the competitor's rises.**
Those are independent facts: ctrlpp's relative residual stays between 1.4e-16
and 5.6e-16 over the whole sweep, while ct's climbs from 1.2e-15 to 9.1e-15.
The agreement column says nothing about this, and cannot: it necessarily tracks
the less exact of the two arms.

The gain form of the same comparison, `ctrlpp::lqr_gain_continuous` against
`ct::optcon::LQR`, tracks the solve almost exactly: **1.25x to 1.40x fewer
instructions** over the same nine rungs, with the two arms' gains agreeing to
between 6.4e-15 and 1.4e-13.  Both arms' gains produce the identical closed-loop
spectral abscissa to twelve or more significant figures on every rung, and every
value is negative, so both gains stabilize the plant.

**On the gain rows, accuracy is one shared upper bound and not a per-arm
discriminator.**  Both arms measure between 1.4e-16 and 9.5e-16, and the reason
is algebraic rather than a matter of judgement: the residual is reconstructed
from the gain, and that reconstruction is **quadratic** in the gain error.  At a
gain error near 1e-15 the true suboptimality term is of order 1e-30, far below
what double precision represents, so what is measured is the reconstruction's
own numerical floor, common to both arms and passing identically through both.
The honest statement is a single bound: **no worse than about 1e-15, identically
for both libraries.**  The solve rows above are the opposite case, where the
per-arm split is a real and reportable difference.

The sweep stops at NX=30 because that is where the library as shipped stops.
The stack-allocation ceiling is a compile-time property of the default build,
and **no compile definition raising it appears in any benchmark source or build
file**, so every rung above is reproducible against ctrlpp as a user receives
it.  A figure that needs a non-default compile definition is a figure a reader
cannot reproduce.

```bash
cmake -S benchmarks -B build-bench -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DCTRLPP_FETCH_BENCHMARK_DEPS=ON -DCTRLPP_BENCH_CT=ON
cmake --build build-bench -j6 --target bench_care_vs_ct bench_lqr_continuous_vs_ct
mkdir -p run && cd run && taskset -c 2 ../build-bench/comparison/ct/bench_care_vs_ct
```

### Discrete Riccati: ctrlpp vs ETH Control Toolbox

Both sides solve the **same** discrete algebraic Riccati equation.  Compiler
GCC 16.1.1 on every rung.

| NX | ctrlpp `dare` | ct `DARE` | **instr x** | Agreement | ctrlpp own residual | ct own residual |
|---:|--------------:|----------:|------------:|----------:|--------------------:|----------------:|
| 2 | 2.63 us (33.5K) | 8.22 us (46.9K) | **1.40** | 5.0e-04 | 2.6e-15 | 3.5e-06 |
| 4 | 11.1 us (159K) | 28.6 us (299K) | **1.89** | 5.1e-04 | 3.9e-15 | 3.6e-06 |
| 6 | 31.9 us (496K) | 52.5 us (644K) | **1.30** | 1.2e-03 | 2.4e-15 | 3.3e-06 |
| 8 | 62.5 us (1.04M) | 145 us (2.13M) | **2.05** | 2.7e-03 | 3.4e-14 | 3.2e-06 |
| 12 | 208 us (3.77M) | 604 us (9.65M) | **2.56** | 2.7e-03 | 3.9e-14 | 3.2e-06 |

**ctrlpp is 1.30x to 2.56x fewer instructions, and part of that advantage is
bought with accuracy that is not ctrlpp's to spend.  Read the two residual
columns before the ratio.**

The agreement column here is **not a rounding figure**.  At 5e-04 to 2.7e-03 it
is eleven orders of magnitude above what the continuous pair reports, and the
own-residual columns say where it comes from.  ct's `DARE` is a **value
iteration** on the time-varying recursion, seeded at `P = Q` and stopped when
the largest coefficient change falls below a **caller-settable default of
1e-6**, with a 1000-iteration cap
(`/usr/include/ct/optcon/lqr/riccati/DARE-impl.hpp:41-58`).  Its relative
residual sits at 3.2e-06 and is flat in problem size, which is the signature of
a stopping tolerance rather than of conditioning.  `ctrlpp::dare` performs a
direct symplectic solve, and its residual is 2.4e-15 to 3.9e-14.

So the ratio is partly a function of a default the reader would otherwise not
see.  At a tighter tolerance ct would do more work, and nothing measured here
says how much more.

The discrete sweep stops at NX=12 for the same shipped-ceiling reason the
continuous sweep stops at NX=30.

One measurement artifact, named rather than smoothed: `ct::optcon::DARE` at NX=4
showed a 30.8% cross-pass spread in wall time (28.6, 37.1 and 28.3 us) on an
instruction count identical to eight significant figures in all three passes.
That is processor boost variance on identical work, and it is the clearest
single reason on this page to read the instruction columns rather than the
times.  Every other row in the window is at or below 4.75%.

```bash
cmake -S benchmarks -B build-bench -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release \
  -DCTRLPP_FETCH_BENCHMARK_DEPS=ON -DCTRLPP_BENCH_CT=ON
cmake --build build-bench -j6 --target bench_dare_vs_ct
mkdir -p run && cd run && taskset -c 2 ../build-bench/comparison/ct/bench_dare_vs_ct
```

### What separates ctrlpp and ct_optcon on the Riccati rows

The two libraries do measurably different amounts of work around the solve, and
these are the differences, stated without attributing any particular share of
the measured ratio to any one of them.

- **Arithmetic**: both factor a `2n x 2n` matrix in real arithmetic.  On the
  continuous rows both apply `Eigen::RealSchur` to the Hamiltonian.  On the
  discrete row `ctrlpp::dare` applies `Eigen::RealSchur` to the symplectic
  matrix (`lib/ctrlpp/include/ctrlpp/control/dare.h`) while ct iterates the
  recursion rather than factoring at all.
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
- **Problem formulation**: the continuous rows compare both sides on the
  continuous equation and the discrete row compares both sides on the discrete
  equation, so each row races like against like.  The symplectic matrix the
  discrete solver builds requires an inverse of `A^T` during setup that the
  Hamiltonian does not; both sides invert `R`.
- **Termination**: ctrlpp's solvers are direct and terminate when the
  factorization is complete.  ct's discrete solver is iterative and terminates
  on a caller-settable tolerance, which is why its discrete residual is flat in
  problem size while ctrlpp's tracks conditioning.
- **Acceptance**: `ctrlpp::dare` verifies the matrix it is about to return
  before returning it -- a definiteness factorization, a forward-error
  estimate against a half-significand criterion, and a closed-loop eigenvalue
  check -- and repeats that verification at the caller's weight scale whenever
  the scale is not one.  ct returns the extracted solution unverified.  That
  verification is a measurable share of the discrete entry point's cost;
  [Discrete Riccati: the acceptance check and the equilibrated
  path](#discrete-riccati-the-acceptance-check-and-the-equilibrated-path)
  measures it on its own.

### What these ratios do not claim

Each row above is one problem family, one dimension ladder, one machine, one
compiler and one dependency set.  None of them is a general statement about
either library, and the direction of the comparison is not uniform: ctrlpp is
ahead on the Riccati family, the predictive-control row and trajectory
synthesis; behind on trajectory evaluation and by two orders of magnitude on
nonlinear predictive control; and at parity on the quadratic program.

## Environment

The [Internal benchmark results](#internal-benchmark-results) and the
[Competitive comparisons](#competitive-comparisons) were taken in a single
measurement window on 2026-08-10, under these conditions.  The
[Discrete Riccati deep dive](#discrete-riccati-the-acceptance-check-and-the-equilibrated-path)
was taken separately and states its own machine settings and toolchain in place.

| | |
|---|---|
| Processor | AMD Ryzen 7 5800X3D, 8 cores / 16 threads |
| Governor | `performance`, read from the machine by the measurement itself rather than assumed |
| Boost | on (`cpufreq/boost = 1`), likewise observed |
| Simultaneous multithreading | active; the pinned core's sibling was idle throughout |
| Pinning | `taskset -c 2`, one core, every run |
| Exclusivity | full.  Load average 0.27 before the first run and 0.82 after the last; processor pressure at or below 0.01 throughout.  No compile or competing process ran during the window. |
| Passes | 3 serial passes per target, medians reported |
| Compiler | GCC 16.1.1 (`-O2`, Release) for every row except the nonlinear predictive control row, which needs GCC 13.4.1 and is marked as such |
| OS | EndeavourOS (Arch Linux), kernel 6.18.42-1-lts |
| Eigen | 3.4.1 |
| OSQP | 1.0.0 |
| osqp-eigen | 0.11.0 |
| ETH Control Toolbox | control-toolbox-optcon 3.0.2 |
| ruckig | 0.9.2 |
| libmpc++ | 1.0.0 |
| nanobench | 4.3.11 |

Warm-up and iteration budget are per target rather than global, because the
targets span six orders of magnitude in cost.  The microbenchmarks use 100
warm-up iterations and a 10000-iteration epoch floor; the Riccati and competitor
comparisons use 50 and 100; the discrete Riccati sweep uses 10 warm-up
iterations and 51 fixed epochs; the nonlinear predictive control row, whose
single operation is a twenty-step closed loop, uses 2 and 5.

**Instructions are the metric to quote.**  They are frequency invariant, so
boost being on does not move them, and they reproduced across the three passes
to 0.00% on nearly every row.  Wall times carry boost variance and are given for
scale, not for precision.

To reproduce with a clock held steadier than the numbers above were taken with,
pin the governor and disable boost before running:

```bash
sudo cpupower frequency-set -g performance
echo 0 | sudo tee /sys/devices/system/cpu/cpufreq/boost
```

Or use `pyperf system tune` as nanobench suggests.  Note that disabling boost
changes the wall times relative to what is published here, which was measured
with boost on; the instruction counts are unaffected.
