# Determinism and RT-Safety Matrix

Per-module real-time safety and determinism guarantees for the ctrlpp hot
paths.  Every YES cell cites the test or build artifact that proves it; no
cell in this table is self-reported.  The evidence was re-run in full before
this table was written.

## Reading the columns

| Column | Meaning |
|--------|---------|
| allocation-free? | The steady-state hot path (compute/predict/update/evaluate/sample) performs zero heap allocation. Construction and setup may allocate. |
| bounded-iterations? | No unbounded loops; any internal iteration carries a fixed cap. |
| wall-clock-free? | The compute path never reads a clock. |
| exceptions-off-clean? | Compiles under `-fno-exceptions -fno-rtti -DCTRLPP_NO_EXCEPTIONS`. |
| deterministic(seeded)? | The output is a pure function of the inputs; any randomness is injectable and seed-deterministic. |
| evidence | The artifact that proves each YES. |

## Evidence artifacts

- **`*_nomalloc_test` targets** (`tests/unit/`): belt-and-suspenders
  allocation guards.  Each test arms a throwing `eigen_assert` (so
  `EIGEN_RUNTIME_NO_MALLOC` cannot be elided) plus a process-global
  `operator new`/`operator delete` counter (so `aligned_malloc` cannot slip
  past), warms the object up, then asserts zero allocations across a
  steady-state loop.  All six targets are registered in
  `tests/unit/CMakeLists.txt` and run serially.  Verified green:
  `ctest --test-dir build-tests -R nomalloc` passes 6 of 6.
- **`scripts/cross_compile_check.sh` leg 1**: host build of the core surface
  with `-fno-exceptions -fno-rtti -DCTRLPP_NO_EXCEPTIONS` at Scalar `double`
  and `float`.  This is the authoritative exceptions-off evidence; verified
  PASS.  Leg 2 (host `-std=c++20`, exceptions on) also passes and covers the
  `ctrlpp::detail::expected` fallback.
- **`tests/compile/embedded_core_float.cpp`**: the exception-free
  instantiation witness.  It compiles the `control`, `estimation`, `dsp`,
  `trajectory`, and `state_space` umbrellas and instantiates at least one
  representative `float` type per module (every estimator individually), is
  built by every leg of the cross-compile script, and passes as the
  `embedded_core_float` ctest.
- **Bare-metal cross-compile (leg 3)**: the same script carries an
  `arm-none-eabi-g++` Cortex-M7 leg (`-mcpu=cortex-m7 -mfpu=fpv5-d16
  -mfloat-abi=hard -fno-exceptions -fno-rtti -DCTRLPP_NO_EXCEPTIONS`) that
  compiles the `float` witness translation unit against the toolchain's own
  hosted libstdc++ subset (no `-nostdinc++`, no `-ffreestanding`).  It passes,
  proving the embedded header subset compiles bare-metal.  The leg requires
  the toolchain's C library headers (the `arm-none-eabi-newlib` package); when
  they are absent the script reports the missing package and exits without
  running leg 3.
- **Clock audit**: a sweep of `lib/ctrlpp/include/` for `<chrono>` and clock
  reads matches only the two optimization solver backends
  (`mpc/nlopt_solver.h`, `mpc/argmin_solver.h`), which carry the opt-in
  `max_time` budget.  Every other core header is clock-free.

## The matrix

| module | allocation-free? | bounded-iterations? | wall-clock-free? | exceptions-off-clean? | deterministic(seeded)? | evidence |
|--------|------------------|---------------------|------------------|-----------------------|------------------------|----------|
| `pid` | YES | YES (no loop) | YES | YES | YES | `pid_nomalloc_test` (position form, velocity form, composed anti-windup and derivative filter); leg 1 + `embedded_core_float` |
| `lqr` / `lqr_time_varying` steady-state control law | YES | YES (no loop in the gain application) | YES | YES | YES | lqr TEST_CASEs in `pid_nomalloc_test.cpp` ("lqr steady-state control law", "lqr_time_varying steady-state control law"); leg 1 witness instantiates `lqr_gain` |
| `kalman_filter` | YES | YES (closed form) | YES | YES | YES | `estimation_nomalloc_test` (kalman_filter case); leg 1 + `embedded_core_float` |
| `ekf` | YES | YES (closed form) | YES | YES | YES | `estimation_nomalloc_test` (ekf case); leg 1 + `embedded_core_float` |
| `ukf` | YES | YES (fixed sigma-point set) | YES | YES | YES | `estimation_nomalloc_test` (ukf case); leg 1 + `embedded_core_float` |
| `mekf` | YES | YES (closed form) | YES | YES | YES | `estimation_nomalloc_test` (mekf case); leg 1 + `embedded_core_float` |
| `manifold_ukf` | YES | YES (fixed sigma-point set) | YES | YES | YES | `estimation_nomalloc_test` (manifold_ukf case); leg 1 + `embedded_core_float` |
| `particle_filter` | YES (guard covers resample and roughening) | YES (fixed particle count) | YES | YES | YES (injected seeded RNG) | `estimation_nomalloc_test` (particle_filter case forces resampling every update; twin filters seeded `std::mt19937_64{42}` must agree bitwise over 64 steps) |
| `complementary_filter` | YES | YES (closed form) | YES | YES | YES | `estimation_nomalloc_test` (complementary_filter case); leg 1 + `embedded_core_float` |
| dsp: `biquad` / `cascaded_biquad` / `vector_biquad` / `fir` | YES | YES (fixed sections and taps) | YES | YES | YES | `dsp_nomalloc_test` (one case per filter); leg 1 + `embedded_core_float` |
| Riccati steady-state solve: `dare` / `care` | YES | YES (Eigen Schur iteration bound; sign-function Newton capped at `max_iters = 40` in `detail/care_sign_function.h`) | YES | YES | YES | `dare_care_nomalloc_test` (NX = 2, 4, 8 across the Schur, sign-function, and balanced-Schur variants); leg 1 witness calls `dare` and `care` |
| trajectory evaluation (polynomial paths, velocity profiles, `cubic_spline`, `smoothing_spline`, `bspline_trajectory`) | YES | YES (closed form; B-spline recursion bounded by compile-time degree) | YES | YES | YES | `trajectory_nomalloc_test` (evaluate cases for cubic/quintic/septic, trapezoidal, double-S, modified sin/trap, cubic_spline, smoothing_spline, bspline_trajectory); leg 1 + `embedded_core_float` |
| online planners: `online_planner_2nd` / `online_planner_3rd` | YES | YES (closed-form segment logic) | YES | YES | YES | `trajectory_nomalloc_test` (update and sample cases for both planners); leg 1 + `embedded_core_float` |
| `recursive_arx` / `rls` | YES | YES (rank-one update, no loop) | YES | not covered (the leg 1 witness does not include the sysid headers) | YES | `sysid_nomalloc_test` (rls and recursive_arx update cases) |
| `mpc` / `nmpc` / `mhe` / `nmhe` | NO (soft real-time: the solve allocates and iterates) | YES when capped (`max_eval`, OSQP `max_iter`) | YES with `max_time = 0` (the default); the `max_time` budget is non-RT | not covered (opt-in OSQP/NLopt/argmin backends sit outside the embedded core witness) | solver-dependent | labeled soft real-time; caps and defaults in `mpc/nlopt_solver.h`, `mpc/argmin_policies.h`, `mpc/osqp_solver.h` |
| static-memory linear MPC | planned | planned | planned | planned | planned | not implemented; the future hard real-time path (see below) |

## Wall-clock budgets are not RT-safe

The optimization solver settings expose two kinds of budget.  The iteration
caps (`nlopt_settings::max_eval`, `argmin_settings::max_eval`, and the OSQP
`max_iter` limit) bound the solve by a deterministic count of steps.  The
wall-clock budget (`nlopt_settings::max_time`, `argmin_settings::max_time`)
bounds the solve by elapsed time, which depends on machine load, frequency
scaling, and cache state, so two identical calls can stop at different
iterates.  The wall-clock budget is therefore **not real-time safe** and must
not be selected on a real-time path.  It defaults to `0` (disabled); the
RT-safe selection is to leave it disabled so only the deterministic iteration
caps bound the solve.

## MPC, NMPC, and MHE are soft real-time

The horizon optimizers allocate and iterate inside the solve: the backends
build their problem structures on the heap and run an iterative QP or NLP
loop to convergence or to the iteration cap.  This is inherent to the
formulation, not an implementation defect, so these modules are labeled
**soft real-time**: latency is bounded in iterations when capped, but the
per-call time and allocation behavior do not meet the hard real-time
obligations the rest of the table certifies.  Callers on a deadline should
budget for the capped worst case and treat the solution status
(`solve_status::max_iterations`) as a first-class outcome.

## Static-memory linear MPC: the future hard real-time path

A fixed-size, iteration-capped QP with statically allocated problem storage
would bring linear MPC under the same hard real-time obligations as the rest
of the matrix (zero steady-state allocation, deterministic iteration bound,
no clock).  That path is planned and explicitly **not implemented**; its row
above is a placeholder so the gap stays visible rather than being rounded up.

## Kernels are passive

ctrlpp modules expose pure `compute`/`predict`/`update`/`evaluate`/`sample`
functions and own no thread, scheduler, timer, run loop, or clock.  This
matrix certifies **RT-safety** (the necessary condition).  **RT-scheduling**
(periods, deadlines, priorities, and WCET budgeting) belongs to the caller:
a bare-metal superloop, an RTOS task, or the host application's executor.
