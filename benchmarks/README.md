# ctrlpp Benchmarks

## Overview

Standalone benchmark suite for ctrlpp. Measures internal hot-path performance
and competitive comparisons against external C++ control libraries.

Uses [nanobench](https://nanobench.ankerl.com/) (v4.3.11) for microbenchmarking
with automatic epoch tuning, statistical analysis, and Linux perf counter
integration.

## Prerequisites

### Required

- CMake >= 3.28
- C++23 compiler (GCC 15+, Clang 18+, MSVC 2022+)
- Eigen3

### Optional (Competitive Benchmarks)

Each competitor is gated behind its own CMake option, all OFF by default.

| Library | CMake Option | Install (Arch/EndeavourOS) | Install (other) |
|---------|-------------|---------------------------|-----------------|
| libmpc++ 1.0.0 | `CTRLPP_BENCH_LIBMPC` | auto (FetchContent) | auto (FetchContent) |
| osqp-eigen 0.11.0 | `CTRLPP_BENCH_OSQP_EIGEN` | `yay -S osqp-eigen` | build from source |
| HPIPM 0.1.3 | `CTRLPP_BENCH_HPIPM` | `yay -S hpipm blasfeo` | build from source |
| ct_optcon 3.0.2 | `CTRLPP_BENCH_CT` | `yay -S control-toolbox` | build from source (catkin heritage) |
| ruckig 0.17.3 | `CTRLPP_BENCH_RUCKIG` | `yay -S ruckig` / auto (FetchContent) | auto (FetchContent) |
| Drake | `CTRLPP_BENCH_DRAKE` | manual install (no AUR, Bazel-only) | pre-built tar.gz (Ubuntu/macOS) |

## Quick Start

```bash
./bench.sh
```

## Building Manually

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

## Running Individual Benchmarks

```bash
./build/internal/bench_pid
./build/internal/bench_lqr
```

## Enabling Competitive Benchmarks

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
    -DCTRLPP_BENCH_RUCKIG=ON \
    -DCTRLPP_BENCH_OSQP_EIGEN=ON
cmake --build build -j$(nproc)
```

Or via bench.sh:

```bash
./bench.sh --competitive ruckig,osqp_eigen
```

To run only internal benchmarks (skip competitive):

```bash
./bench.sh --internal-only
```

## Performance Counters (Linux)

nanobench reports instructions, branches, and cache misses via `perf_event`.

Requires:

```bash
sudo sysctl kernel.perf_event_paranoid=-1
```

Or run benchmarks with `sudo`.

If `perf_event` is unavailable, benchmarks gracefully degrade to time-only.

## Output

- **Console:** nanobench default table (ns/op, op/s, error%)
- **CSV:** comma-separated files in `results/` directory after `bench.sh` run

## Hardware Reference

Benchmarks developed on AMD Ryzen 7 5800X3D, GCC 15.2.1, Eigen 5.0.1.
Results will vary by hardware -- always compare on the same machine.

## Notes on Drake

Drake uses Bazel as its build system and has no CMake FetchContent support.
It is not available via AUR. The `CTRLPP_BENCH_DRAKE` gate requires a manual
Drake installation with CMake `find_package` support. On Arch Linux, this is
impractical without significant effort.
