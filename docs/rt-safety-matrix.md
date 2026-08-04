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
  allocation guards.  Each test arms a throw-free `eigen_assert` that sets a
  pollable sentinel (`ctrlpp_test::detail::eigen_alloc_violation`, so
  `EIGEN_RUNTIME_NO_MALLOC` cannot be elided and the trap compiles under
  `-fno-exceptions`) plus a process-global `operator new`/`operator delete`
  counter (so `aligned_malloc` cannot slip past), warms the object up, then
  asserts zero allocations and no sentinel violation across a steady-state
  loop.  All six targets are registered in `tests/unit/CMakeLists.txt`, run
  serially, and compile in the default `-fno-exceptions` tree.  Verified green:
  `ctest --test-dir build/dev -R nomalloc` passes 6 of 6.
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
| `particle_filter` | YES (guard covers resample, roughening and the covariance read) | YES (fixed particle count) | YES | YES | YES (injected seeded RNG) | `estimation_nomalloc_test` (particle_filter case forces resampling every update and reads `covariance()` inside the armed window; twin filters seeded `std::mt19937_64{42}` must agree bitwise over 64 steps in both the estimate and the reported uncertainty) |
| `complementary_filter` | YES | YES (closed form) | YES | YES | YES | `estimation_nomalloc_test` (complementary_filter case); leg 1 + `embedded_core_float` |
| dsp: `biquad` / `cascaded_biquad` / `vector_biquad` / `fir` | YES | YES (fixed sections and taps) | YES | YES | YES | `dsp_nomalloc_test` (one case per filter); leg 1 + `embedded_core_float` |
| Riccati steady-state solve: `dare` / `care` | YES **on the heap, and the heap is not the binding cost here** (see the stack note below this table, which a hard-real-time caller must read before sizing a task stack: at `NX = 8` the whole chain reaches 42,760 bytes) | YES (Eigen Schur iteration bound; sign-function Newton capped at `max_iters = 40` in `detail/care_sign_function.h`, followed by a fixed-size rescale-and-factorize extraction. On every accepted solve, on **every** continuous method tag and on the discrete solver, fixed-size checks of the returned solution run before it is reported: the discrete path solves a fixed-size Stein system of dimension `NX(NX+1)/2` on the symmetric subspace to estimate the answer's own forward error, and every continuous path evaluates the counted Riccati residual and one non-accumulating real Schur factorization of the closed-loop spectrum, from the single definition in `detail/care_postconditions.h`. All three continuous tags reach that definition -- the sign-function path, the real-Schur path and the balanced-Schur path, the last verifying against the caller's Hamiltonian rather than the balanced one -- so no tag carries a bound the others do not. The balanced tag additionally runs a DGEBAL-style balance ahead of all of this. It allocates nothing, and it terminates because every accepted rescale strictly reduces that row's norm sum by the factor 19/20 using power-of-two scalings clamped to the scalar type's range. Its sweep count is data-dependent, and it is now **bounded at compile time** the way `max_iters` bounds the Newton loop: the cap is `2 * NX * (max_exponent - min_exponent)`, derived from the scalar type and the dimension rather than tuned, and reaching it stops the sweep and returns the partially balanced matrix. Stopping early is safe because balancing is a similarity preconditioner and every applied step updates `H` and `D` together, so `H_returned == D^-1 * H_original * D` holds exactly at any cut point -- a capped run is less well balanced, never wrong. The cap is a guarantee rather than a budget: over 550,000 draws with entries spread across `2^-280` to `2^+280`, the worst sweep counts observed were 2 / 12 / 43 / 45 at NX = 1 / 2 / 4 / 8 against ceilings of 4,090 / 8,180 / 16,360 / 32,720, so the worst measured run sits about three orders of magnitude below its ceiling. Budget from the measurement and treat the cap as the backstop) | YES | YES | YES | `dare_care_nomalloc_test` (NX = 2, 4, 8 across the Schur, sign-function, and balanced-Schur variants); `care_convergence_anchor_test` (fixed-seed near-axis and simple-input scale sweeps under all three method tags, plus magnitude-band anchors); `riccati_magnitude_test` (both ends of the arithmetic range); leg 1 witness calls `dare` and `care`. The allocation cell's heap claim is proved by `dare_care_nomalloc_test`; its stack claim is proved by the `-fstack-usage` frames and runtime watermarks in the stack note below, which are host `double` measurements, now corroborated on silicon by the ESP32 stack probe in `examples/embedded/esp32/main/app_main.cpp` -- the on-target `float` peaks track the host `double` prediction to within the scalar width, and that probe also establishes that at `float` the ACCEPTANCE CHECK, not the stack, is what bounds the usable state dimension |
| velocity profile construction: `trapezoidal_trajectory::create` / `double_s_trajectory::create` | YES (the profile is returned by value inside a `ctrlpp::expected`; nothing on the path owns storage) | YES (closed form on the trapezoidal path and on every double-S path but one; the cruise-free double-S rise is a bracket halved to exhaustion, bounded by one more than the significand width, so 25 evaluations for `float` and 54 for `double`) | YES | YES | YES (no random source; identical inputs exhaust the bracket at the identical step) | `trajectory_nomalloc_test` (construction is outside the armed window, as for every other type there); leg 1 + `embedded_core_float`, which instantiate both `create` factories |
| trajectory evaluation (polynomial paths, velocity profiles, `cubic_spline`, `smoothing_spline`, `bspline_trajectory`) | YES | YES (closed form; B-spline recursion bounded by compile-time degree) | YES | YES | YES | `trajectory_nomalloc_test` (evaluate cases for cubic/quintic/septic, trapezoidal, double-S, modified sin/trap, cubic_spline, smoothing_spline, bspline_trajectory); leg 1 + `embedded_core_float` |
| online planners: `online_planner_2nd` / `online_planner_3rd` | YES | YES (closed-form segment logic) | YES | YES | YES | `trajectory_nomalloc_test` (update and sample cases for both planners); leg 1 + `embedded_core_float` |
| trajectory time scaling: `rescale_to` / `can_rescale_to` / `synchronize` (trapezoidal, double-S) | YES (two-pass over a `std::span`, no owning copy, no allocation; the trapezoidal solve's boundary-duration helper returns a two-scalar aggregate by value and owns no storage) | YES (closed form on the trapezoidal and the rest-to-rest double-S paths -- the trapezoidal plateau and valley shapes are each a single quadratic in the cruise velocity's distance from that shape's own boundary, straight-line with no iteration; bracket exhaustion on the nonzero-boundary-velocity double-S path, bounded by one more than the significand width, so 25 evaluations for `float` and 54 for `double`, plus one step per binary exponent the bracket spans) | YES | YES | YES (no random source; identical inputs exhaust the bracket at the identical step) | `trajectory_nomalloc_test`; `trajectory_rescale_anchor_test` and the rescaling cases in `trajectory_hardening_test`; leg 1 + `embedded_core_float` |
| `recursive_arx` / `rls` | YES (the update's result is a `ctrlpp::expected<void, rls_update_error>` holding one enumerator and no owning member; the two norms feeding the resolution floor are unevaluated Eigen expressions over existing storage) | YES (rank-one update, no loop; the refusal guard adds a fixed count of reads per cycle, every dimension a compile-time template parameter -- `NP*NP + NP` for the carried-state scan, `NP` for the regressor scan, and two `NP`-term norms for the denominator's scale) | YES | not covered (the leg 1 witness does not include the sysid headers) | YES (the guard is a pure predicate on the operands and the carried members; no random source, no clock) | `sysid_nomalloc_test` (rls and recursive_arx update cases, re-run green with the guard on every cycle of both 256-cycle armed windows); the typed refusals asserted in `sysid_hardening_test` |
| `mpc` / `nmpc_dynamic` / `mhe` / `nmhe` (runtime-horizon, dynamic solver) | NO (soft real-time: the solve allocates and iterates) | YES when capped (`max_eval`, OSQP `max_iter`) | YES with `max_time = 0` (the default); the `max_time` budget is non-RT | not covered (opt-in OSQP/NLopt/argmin backends sit outside the embedded core witness) | solver-dependent | labeled soft real-time; caps and defaults in `mpc/nlopt_solver.h`, `mpc/argmin_policies.h`, `mpc/osqp_solver.h`, `mpc/argmin_qp_solver.h` |
| `nmpc` (the DEFAULT; = `nmpc_static`, compile-time horizon, argmin `nw_sqp`, bounded decision `NV` + constraint `MaxM`) | YES — strict-zero: 0.00 allocs/step in steady state | YES when capped (`max_eval`) | YES with `max_time = 0` (the default); the `max_time` budget is non-RT | not covered by this witness (argmin's `-fno-exceptions` instantiation is clean upstream; the ctrlpp-side no-exceptions dogfood is a separate witness) | YES — the constraint bound feeds only the QP result-multiplier storage, never the compute workspace, so argmin's `nw_sqp` bit-identity golden is unchanged | `nmpc_static_nomalloc_test` (double_integrator NX=2 NU=1 NH=5 → NV=17, MaxM=12; throwing `eigen_assert` + `EIGEN_RUNTIME_NO_MALLOC` + global `operator new` counter; `static_assert(strict_allocation_free)`) |
| static-memory linear MPC | planned | planned | planned | planned | planned | not implemented; the future hard real-time path (see below) |
| the estimator `update` rejection guard (`kalman_filter`, `ekf`, `ukf`, `mekf`, `manifold_ukf`, `luenberger_observer`, `complementary_filter`) | YES (an `allFinite()` scan is an unevaluated Eigen expression over existing storage; the result is a `ctrlpp::expected<void, E>` holding one enumerator and no owning member) | YES (a fixed count of reads per step: state + covariance + measurement, every dimension a compile-time template parameter -- `NX + NX*NX + NY` for the covariance filters, `NX + NY` for the observer, `4 + NB + NE*NE + NY` for the MEKF, at most `7 + 9 + 1` for the complementary filter) | YES | YES | YES (a pure predicate on the operands; no random source and no state read beyond the members it scans) | `estimation_nomalloc_test` re-run green with the guard on every step of all six converted cases (128-step armed window each, 0 allocations, sentinel clean); the four-part rejection asserted in `{kalman,ekf,ukf,mekf,manifold_ukf,luenberger,complementary_filter}_hardening_test`; `scripts/cross_compile_check.sh` all three legs PASS |
| the controller step rejection guard (`pid::compute` both overloads, `mrac_controller::evaluate`, `l1_controller::evaluate`) | YES (an `allFinite()` scan is an unevaluated Eigen expression over existing storage; the result is a `ctrlpp::expected<vector_t, E>` whose payload is the same by-value vector the surface already returned, so no owning member is added) | YES (a fixed count of reads per cycle, every dimension a compile-time template parameter: `pid` scans only the members its policy composition makes live, at most `11*NY` plus one scalar test on the step; `mrac` scans `NX + NU` vector and `NU*(NX + NU)` matrix entries; `l1` scans `2*NX + 2*NU`. The `l1` cycle additionally scans `NU` more for the pre-projection finiteness test that feeds `health()`) | YES | YES | YES (a pure predicate on the operands and the carried members; no random source, no clock) | `pid_nomalloc_test` re-run green with the guard on every cycle of all three composed cases (256-cycle armed window each, 0 allocations, sentinel clean); the four-part rejection asserted in `{pid,mrac,l1}_hardening_test`, each proven to fail with its guard deleted; both standing trees green at baseline (90/90, 129/129) |

## Stack cost of the Riccati solve

The `allocation-free?` column is a **heap** statement. For most rows in this
table that is the whole resource story, because their hot paths hold a handful
of fixed-size matrices. It is not the whole story for the Riccati solve: it
runs a fixed-size decomposition whose dimension is the caller's `NX`, its frames
grow as `NX^3`, and on the four-to-sixteen-kibibyte task stacks this milestone's
own targets run, **the stack is what will break a caller, not the heap.**

Two numbers per state dimension. The first is the acceptance check's own frame,
as the compiler's `-fstack-usage` report gives it. The second is the deepest
disturbed word measured at runtime over the whole `dare` chain, which is the
number a task stack must actually cover. Both are given for the shipped header,
alongside the same chain with no acceptance check on it at all, so the cost of
the check is separable from the cost of the solve.

| `NX` | acceptance-check frame | whole-chain peak, shipped | whole-chain peak, no check at all |
|---:|---:|---:|---:|
| 2 | 864 | 5,352 | 5,352 |
| 4 | 2,400 | 11,688 | 9,416 |
| 6 | 6,432 | 23,144 | 17,112 |
| 8 | 15,056 | **42,760** | 28,168 |

Supported maximum state dimension, from the whole-chain peak, strict: no margin
for the caller's own frames and none for RTOS overhead.

| task stack | supported `NX`, shipped | supported `NX`, no check at all |
|---|---|---|
| 4 KiB | **none** | **none** |
| 8 KiB | `NX <= 2` | `NX <= 2` |
| 16 KiB | `NX <= 4` | `NX <= 4` |
| 32 KiB | `NX <= 6` | `NX <= 8` |
| 48 KiB | `NX <= 8` | `NX <= 8` |
| 64 KiB | `NX <= 8` | `NX <= 8` |

**Read the first three rows before the last three.** Across the whole
four-to-sixteen-kibibyte band the supported maximum is the same under every
construction including none at all, so it is a property of the solve rather than
of the check: `dare` needs 17,112 bytes at six states and 28,168 at eight with
no acceptance check whatsoever. Turning the check off does not buy a small-stack
target a larger problem. It changes the answer in exactly one band, at 32 KiB.

**Provenance.** These are HOST measurements: `g++ (GNU) 16.1.1 20260728`,
`-std=c++20 -O2`, no `-march` (driver default `-mtune=generic -march=x86-64`),
x86-64 Linux, `double`, runtime watermarks taken on a pthread with a 64 MiB
stack against a zero-byte harness floor. A target's own frames differ with its
ABI, register file and calling convention.

### On silicon, and what it does to the table above

The host figures above were a PREDICTION. They have now been tested on a board,
and they hold once the one variable that separates the two measurements is
accounted for: **the host table is `double` and the board runs `float`.**

Measured on an ESP32-WROOM (Xtensa LX6, 240 MHz, single-precision FPU), ESP-IDF
v6.0.2, `xtensa-esp32-elf-g++` 15.2.0 at `-Os` with exceptions and RTTI off,
Eigen 3.4.0, one FreeRTOS task per dimension so each figure is that dimension's
own peak rather than a running minimum over a sweep:

| `NX` | on-silicon peak, `float` | host peak, `double` | host / board | accepted? |
|---:|---:|---:|---:|:--|
| 2 | 3,100 | 5,352 | 1.73x | solved |
| 3 | 4,160 | -- | -- | refused |
| 4 | 5,932 | 11,688 | 1.97x | solved |
| 5 | 8,336 | -- | -- | refused |
| 6 | 11,500 | 23,144 | 2.01x | refused |
| 8 | 20,560 | 42,760 | 2.08x | refused |

**The ratio converges on exactly the scalar width.** It is 1.73x at `NX = 2` and
climbs to 2.08x at `NX = 8`, which is what a halved scalar predicts for a chain
whose cost is dominated by `NX`-dimensioned arrays and diluted at small `NX` by
fixed overhead that does not scale with the scalar. The host table is therefore
usable as written for `double` and halves for `float`; it is not refuted and it
is not to be replaced by these figures.

**Supported maximum on this board: `NX = 4`, leaving 2,256 bytes of an 8,192-byte
task stack.** Established twice and by different means. A pass giving every
dimension 49,152 bytes recorded `NX = 5` needing 8,336 bytes, which is 144 more
than the task stack has; a second pass at the control task's own 8,192 bytes then
walked the same dimensions and the board reported

    Debug exception reason: Stack canary watchpoint triggered (probe_nx5)

so the overflow is OBSERVED AND NAMED at the instruction that caused it, not
inferred from a reset. Both stack-overflow detection modes are enabled in
`examples/embedded/esp32/sdkconfig.defaults`; the end-of-stack watchpoint is what
fired, and it is the more precise of the two because the canary check only runs
at a context switch. Enabling the watchpoint costs up to 60 bytes of every task's
usable stack, so these figures and any taken without it are not interchangeable.

**THE STACK IS NOT THE BINDING LIMIT FOR `float` ON THIS BOARD. THE ACCEPTANCE
CHECK IS.** Only `NX = 2` and `NX = 4` were accepted. Every other probed
dimension returned `arithmetic_limit` -- "solution is not reliable at this
precision" -- including `NX = 6`, whose 11,500 bytes would fit a 16 KiB task
comfortably. At `float` the half-significand margin is `sqrt(eps) = 3.45e-4`, and
on this corpus the estimated forward error exceeds it from `NX = 5` up and at
`NX = 3`. A caller sizing a 16 KiB task for a six-state `float` plant would find
the stack sufficient and the answer refused. That is the gate working as
designed, and it is the number to plan against.

Witness: `examples/embedded/esp32/main/app_main.cpp`, whose control loop still
designs its gain, runs its 201 steps, streams them over UART2 and reports
`golden diff PASS` in the same boot that carries the probe.

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

The runtime-horizon optimizers allocate and iterate inside the solve: the
backends build their problem structures on the heap and run an iterative QP or
NLP loop to convergence or to the iteration cap.  This is inherent to the
formulation, not an implementation defect, so `nmpc_dynamic` and the runtime-
horizon `mpc` / `mhe` / `nmhe` are labeled **soft real-time**: latency is bounded
in iterations when capped, but the per-call time and allocation behavior do not
meet the hard real-time obligations the rest of the table certifies.  Callers on
a deadline should budget for the capped worst case and treat the solution status
(`solve_status::max_iterations`) as a first-class outcome.

The **default** `nmpc` is the exception (it aliases the compile-time-horizon
`nmpc_static`): pinning the horizon at compile time fixes the decision dimension
`NV`, and binding a compile-time constraint cap `MaxM` moves argmin's `nw_sqp`
result-multiplier storage inline, so the steady-state solve allocates nothing
(measured 0.00 allocs/step) once warmed up.  It remains iteration-bounded by
`max_eval` and clock-free with `max_time = 0`, so on those axes it meets the hard
real-time obligations the runtime-horizon `nmpc_dynamic` path cannot.

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

## Exceptions posture: consumers flag-agnostic, own tests -fno-exceptions-enforced

The library is **consumer-flag-agnostic**: no `-fno-exceptions` / `-fno-rtti`
is forced onto any installed, interface, or exported target, and `config.h`'s
`__cpp_exceptions` auto-detection adapts to whatever the consumer compiles with.
The controller and estimator construction paths converted for runtime
validation use fallible factories returning `ctrlpp::expected`; those factories
are available in every build mode.  The optional OSQP, NLopt and argmin backend
adapters likewise expose a single fallible `setup(problem)` returning
`ctrlpp::expected<void, E>`.  These claims are scoped to the converted
construction and setup APIs, not to every public type or accessor in the
library.

As a self-imposed compatibility guarantee, ctrlpp's **own** tests and benches
dogfood the throw-free discipline: the default build tree
(`CTRLPP_TESTS_WITH_EXCEPTIONS=OFF`, the `dev` preset) compiles ctrlpp's test
targets **and** Catch2 under `-fno-exceptions -fno-rtti`, with Catch2 built
`CATCH_CONFIG_DISABLE_EXCEPTIONS`.  The blocking `no-exceptions` CI job builds
that tree and runs the full ctest **including the allocation-free (no-malloc)
suite**, so the allocation guards and the compiler-enforced prohibition on
throw expressions are exercised together in the same shipping configuration.

Two caveats follow from the disabled-exceptions Catch2:

- `REQUIRE_THROWS*` / `CHECK_THROWS*` are unavailable and a failing assertion
  reports-and-aborts instead of throwing.  The `[!shouldfail]` meta-test and any
  test that reaches a throwing third-party backend therefore live in the
  separate **exceptions carve-out tree** (`CTRLPP_TESTS_WITH_EXCEPTIONS=ON`, the
  `exceptions` preset), never deleted, and are built+run by the `exceptions` CI
  job.  Converted construction and solver setup stay fallible in both trees,
  but `ctrlpp::expected::value()` remains a checked convenience whose bad-access
  path throws in this exceptions carve-out.
- The OSQP and NLopt solver backends throw internally, so every OSQP/NLopt-linked
  test and comparison bench requires the exceptions build; argmin's static NMPC
  path is throw-free and runs in the default `-fno-exceptions` tree.

## The expected result type is owned, not std::expected

`ctrlpp::expected` is a single, always-on C++20 implementation the library owns on
every toolchain; it never aliases `std::expected`, even where the standard library
ships it.  This keeps one controlled embedded floor: storage is a raw discriminated
union rather than `std::variant` (no `<variant>`, no `bad_variant_access` /
valueless-by-exception machinery), `operator*` and `error()` are unchecked, and the
only throw site — `value()` on an error — is gated behind `__cpp_exceptions` with a
`std::abort()` fallback, so the header compiles clean under `-fno-exceptions -fno-rtti`
and is built and exercised by the default `dev` tree, not only the exceptions tree.
This does not make `value()` exception-free: erroneous access throws when
exceptions are enabled and aborts when they are disabled.  Callers must branch on
the `expected` before using `operator*` or `operator->`; code that needs a checked
accessor must keep `value()` out of a real-time failure path.

It is a faithful-API result type, **not** bit-for-bit `std::expected`: the owned copy
and move special members make it non-trivially-copyable even when `T` and `E` are both
trivial, so triviality is not propagated to the ABI.  In practice `T` is usually a
non-trivial Eigen type, so the cost rarely bites; the conditional-triviality path is
deliberately deferred until a caller needs a trivially-relocatable expected.  Boundary
interop with `std::expected` / `std::unexpected` is provided through explicit converting
constructors and conversion operators, guarded by `__cpp_lib_expected`, so results still
cross any `std` boundary without forcing the owned type onto that boundary.
