# Benchmarks

Performance measurements for ctrlpp components and comparison comparisons
against other C++ control libraries.

## Running benchmarks

The benchmark suite is a standalone CMake project under `benchmarks/`.

### Prerequisites

- C++23 compiler
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
| `n4sid` | `identify()` | 578 us | 1.7K | 6.9M |

System dimensions: NX=2 for estimators and sysid, NX=4/NU=2 for LQR/DARE,
5-tap for FIR, 5-knot natural cubic spline.

## Competitive comparisons

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

ct_optcon's LQR solver is faster because it uses a different algebraic
Riccati equation solver.  The key differences:

- **Arithmetic**: ct uses real Schur decomposition (`Eigen::RealSchur`) on
  the Hamiltonian matrix with real `double` arithmetic.  ctrlpp uses complex
  Schur decomposition (`Eigen::ComplexSchur`) on the symplectic matrix with
  `std::complex<double>` — double the memory footprint and more expensive
  multiply-accumulate operations.
- **Schur reordering**: ct calls LAPACK `dtrsen` directly, a
  decades-optimized Fortran routine for reordering real Schur forms.  ctrlpp
  implements its own Givens rotation swaps in C++ for complex Schur
  reordering.
- **Problem formulation**: ct solves the continuous-time ARE (CARE) while
  ctrlpp solves the discrete-time ARE (DARE).  The symplectic matrix for
  DARE requires an additional matrix inverse during setup that the
  Hamiltonian for CARE does not.

This is a known gap tracked as PERF-F01.  The planned fix is to switch the
DARE solver to real Schur decomposition with LAPACK `dtrsen` for reordering.

## Environment

Results above were collected on:

- CPU: AMD (Zen-class, frequency scaling enabled &mdash; results may vary)
- OS: EndeavourOS (Arch Linux), kernel 6.18
- Compiler: GCC 15.2.1, `-O2` (Release)
- Eigen: 3.4+
- OSQP: 1.0.0
- nanobench: 4.3.11

For stable results, disable CPU frequency scaling before benchmarking:

```bash
sudo cpupower frequency-set -g performance
```

Or use `pyperf system tune` as nanobench suggests.
