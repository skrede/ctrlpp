# Determinism and RT-Safety Matrix

Per-module real-time safety and determinism guarantees for the ctrlpp hot
paths.  Every YES cell cites the test or build artifact that proves it; no
cell in this table is self-reported.  The evidence was re-run in full before
this table was written.

## Reading the columns

| Column | Meaning |
|--------|---------|
| allocation-free? | The steady-state hot path (compute/predict/update/evaluate/sample) performs zero heap allocation. Construction and setup may allocate. **THIS COLUMN IS A HEAP STATEMENT AND SAYS NOTHING ABOUT THE STACK.** For a row whose hot path holds a handful of fixed-size matrices that is the whole resource story. For a row that runs a fixed-size decomposition whose dimension the CALLER chooses, it is not: those frames grow with that dimension, and on a small task stack the stack is what breaks first. Every such row has a stack section below the table -- see "Stack cost of the Riccati solve", "Stack cost of the estimator rows" and "Stack cost of the predictive controller row" -- and a `YES` here must not be read as a resource claim on its own. |
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
  counter, warms the object up, then asserts an exact allocation count across a
  steady-state loop.  All eleven targets are registered in
  `tests/unit/CMakeLists.txt`, run serially because the counter is
  process-global, and compile in the default `-fno-exceptions` tree; a twelfth,
  `nmpc_static_nomalloc_test`, needs the optional nonlinear-programming backend
  and is registered with it.  Verified green: `ctest --test-dir build/dev -R
  nomalloc` passes 11 of 11.
  - **Neither mechanism substitutes for the other, and that is measured on this
    library rather than argued.**  The allocations the caller-dimensioned rows
    make above the boundary described below go through Eigen's `aligned_malloc`,
    which calls `std::malloc` directly, so **the counter reads zero at every one
    of them** and a guard carrying only the counter would report every one
    allocation-free.  The counter covers the converse case, every non-Eigen heap
    allocation Eigen's own bookkeeping cannot see.
  - **Both mechanisms are demonstrated to fire.**  `estimation_nomalloc_test`
    and `dare_care_nomalloc_test` each carry a negative control with one case
    per mechanism: a direct call to the replaced allocation function, which
    unlike a new-expression cannot be elided, and an Eigen allocation inside an
    armed window.  Without it a grid can report zero at every point while
    neither mechanism is capable of firing.
  - **The whole family is SKIPPED under MemorySanitizer.**  That runtime ships
    its own strong `operator new`/`delete`, which collide at link with the
    counting replacements, so the targets are excluded there rather than having
    the harness drop the replacements and pass without testing anything.  A
    claim that this project's evidence includes a MemorySanitizer leg must say
    that the no-allocation grid is not in it.
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
  `max_time` budget.  `argmin_solver` instantiates a clock read only when it is
  built with `argmin_step_and_time_budget`; under its default step budget it
  reads none.  Every other core header is clock-free.

## The matrix

| module | allocation-free? | bounded-iterations? | wall-clock-free? | exceptions-off-clean? | deterministic(seeded)? | evidence |
|--------|------------------|---------------------|------------------|-----------------------|------------------------|----------|
| `pid` | YES | YES (no loop) | YES | YES | YES | `pid_nomalloc_test` (position form, velocity form, composed anti-windup and derivative filter); leg 1 + `embedded_core_float` |
| `lqr` / `lqr_time_varying` steady-state control law (`lqr_time_varying::compute` refuses a step at or past its horizon on the `ctrlpp::expected` failure channel) | YES | YES (no loop in the gain application) | YES | YES | YES | lqr TEST_CASEs in `pid_nomalloc_test.cpp` ("lqr steady-state control law", "lqr_time_varying steady-state control law"); leg 1 witness instantiates `lqr_gain` |
| `kalman_filter` | YES below a measurement dimension of 48, **NO at and above it** at two states or more, where the gain solve's Householder sequence goes blocked and the linear-algebra library heap-allocates its block reflector factor once per update (see "The heap claim's boundary" below). A single-state filter stays allocation-free at every measurement dimension. | YES (closed form) | YES | YES | YES | `estimation_kalman_nomalloc_test` (26 configurations spanning 1 to 128 on both axes: the 13 its measurement table prints, the 6 values only its supported-maximum table names, the state axis at its instantiation limit, both sides of the boundary at 47 and 48, and the single-state exemption); `estimation_nomalloc_test` (the fixed instantiation at two states and one output); leg 1 + `embedded_core_float`. On silicon at `double` and that same instantiation, the NUCLEO-H753ZI evidence image reports 0 allocations over 201 steps behind a passing allocation canary, with the covariance entry `P(1,1)`, the velocity variance, after its last step within the derived bound of the host golden; the state estimate is not compared (`scripts/board_capture_fixtures/evidence.txt`, `arm-none-eabi-gcc` 16.2.0; see "On silicon: the NUCLEO-H753ZI leg" below). At `float` and the same dimensions, the ESP32-POE-EA evidence image reports 0 allocations over 201 steps behind its passing allocation and heap-view canaries, with a clean heap view in the window, and its covariance entry `P(1,1)` after the last step within the single-precision bound derived for it (`scripts/board_capture_fixtures/esp32-evidence.txt`, `xtensa-esp-elf-gcc` 15.2.0 under ESP-IDF v6.0.2; see "On silicon, and what it does to the table above") |
| `ekf` | YES below a measurement dimension of 48, **NO at and above it** at two states or more, for the same reason and with the same single-state exemption as the row above. | YES (closed form) | YES | YES | YES | `estimation_ekf_nomalloc_test` (25 configurations spanning 1 to 128 on both axes: the 13 its measurement table prints, the 6 values only its supported-maximum table names, the state axis at its instantiation limit, both sides of the boundary at 47 and 48, and the single-state exemption); `estimation_nomalloc_test` (the fixed instantiation at two states and one output); leg 1 + `embedded_core_float` |
| `ukf` | YES at every measured configuration, including 128 states by 128 outputs, **because its gain decomposition defaults to LDLT and builds no Householder sequence at all**. Selecting the QR option puts it on the boundary the rows beside it have, at a measurement dimension of 48. | YES (fixed sigma-point set) | YES | YES | YES | `estimation_ukf_nomalloc_test` (25 configurations spanning 2 to 128 on both axes, including both instantiation limits and the corner, plus 2 more with the QR option selected, walking both sides of the boundary that option brings back); `estimation_nomalloc_test` (the fixed instantiation at two states and one output); leg 1 + `embedded_core_float` |
| `mekf` | YES below a measurement dimension of 48, **NO at and above it** (see "The heap claim's boundary" below). The boundary is on the measurement dimension and not on the error state: 125 bias states against 4 outputs is allocation-free. This row has no single-state exemption, since its right-hand side is the error state and is never one column. | YES (closed form) | YES | YES | YES | `estimation_mekf_nomalloc_test` (20 configurations: the 11 of the 13 its measurement table prints that instantiate, the 4 values only its supported-maximum table names, both axes at their instantiation limits with the corner, and both sides of the boundary at 47 and 48); `estimation_nomalloc_test` (the fixed instantiation at three bias states and three outputs); leg 1 + `embedded_core_float` |
| `manifold_ukf` | YES below a measurement dimension of 48, **NO at and above it** (see "The heap claim's boundary" below). No single-state exemption: the right-hand side is the rotation's three columns. | YES (fixed sigma-point set) | YES | YES | YES | `estimation_manifold_ukf_nomalloc_test` (12 configurations spanning 2 to 128 on its one caller axis: the 5 its measurement table prints, the 4 values only its supported-maximum table names, the axis at its instantiation limit, and both sides of the boundary at 47 and 48); `estimation_nomalloc_test` (the fixed instantiation at three outputs); leg 1 + `embedded_core_float` |
| `particle_filter` | YES (guard covers resample, roughening and the covariance read). `update` returns a typed `ctrlpp::expected<void, pf_update_error>`, and the new weights are staged in a member array rather than a local one, so the update's stack frame does not grow with the particle count | YES (fixed particle count) | YES | YES | YES (injected seeded RNG; a refused update draws nothing from it) | `estimation_nomalloc_test` (particle_filter case forces resampling every update, checks every update's result and reads `covariance()` inside the armed window; twin filters seeded `std::mt19937_64{42}` must agree bitwise over 64 steps in both the estimate and the reported uncertainty); `embedded_core_float` steps it once at `float` |
| `complementary_filter` | YES | YES (closed form) | YES | YES | YES | `estimation_nomalloc_test` (complementary_filter case); leg 1 + `embedded_core_float` |
| dsp: `biquad` / `cascaded_biquad` / `vector_biquad` / `fir` | YES | YES (fixed sections and taps) | YES | YES | YES | `dsp_nomalloc_test` (one case per filter); leg 1 + `embedded_core_float`. On silicon at `double`, a `biquad` low-pass on the NUCLEO-H753ZI evidence image reports 0 allocations over 201 samples behind a passing allocation canary, and its last step's value equals the host golden of the section formation that image was built from (`scripts/board_capture_fixtures/evidence.txt`, `arm-none-eabi-gcc` 16.2.0; see "On silicon: the NUCLEO-H753ZI leg" below). Neither board has been captured with the current formation, whose host golden differs from that one in the last two digits. At `float`, the same `biquad` on the ESP32-POE-EA evidence image reports 0 allocations over 201 samples behind its passing allocation and heap-view canaries, with a clean heap view in the window, and its last step's value within the single-precision bound derived for it (`scripts/board_capture_fixtures/esp32-evidence.txt`, `xtensa-esp-elf-gcc` 15.2.0 under ESP-IDF v6.0.2; see "On silicon, and what it does to the table above") |
| Riccati steady-state solve: `dare` / `care` | YES **on the heap, and the heap is not the binding cost here** (see the stack note below this table, which a hard-real-time caller must read before sizing a task stack: at `NX = 8` the whole chain reaches 41,304 bytes) | YES (Eigen Schur iteration bound; sign-function Newton capped at `max_iters = 40` in `detail/care_sign_function.h`, followed by a fixed-size rescale-and-factorize extraction. On every accepted solve, on **every** continuous method tag and on the discrete solver, fixed-size checks of the returned solution run before it is reported: the discrete path solves a fixed-size Stein system of dimension `NX(NX+1)/2` on the symmetric subspace to estimate the answer's own forward error, and every continuous path evaluates the counted Riccati residual and one non-accumulating real Schur factorization of the closed-loop spectrum, from the single definition in `detail/care_postconditions.h`. All three continuous tags reach that definition -- the sign-function path, the real-Schur path and the balanced-Schur path, the last verifying against the caller's Hamiltonian rather than the balanced one -- so no tag carries a bound the others do not. The balanced tag additionally runs a DGEBAL-style balance ahead of all of this. It allocates nothing, and it terminates because every accepted rescale strictly reduces that row's norm sum by the factor 19/20 using power-of-two scalings clamped to the scalar type's range. Its sweep count is data-dependent, and it is now **bounded at compile time** the way `max_iters` bounds the Newton loop: the cap is `2 * NX * (max_exponent - min_exponent)`, derived from the scalar type and the dimension rather than tuned, and reaching it stops the sweep and returns the partially balanced matrix. Stopping early is safe because balancing is a similarity preconditioner and every applied step updates `H` and `D` together, so `H_returned == D^-1 * H_original * D` holds exactly at any cut point -- a capped run is less well balanced, never wrong. The cap is a guarantee rather than a budget: over 550,000 draws with entries spread across `2^-280` to `2^+280`, the worst sweep counts observed were 2 / 12 / 43 / 45 at NX = 1 / 2 / 4 / 8 against ceilings of 4,090 / 8,180 / 16,360 / 32,720, so the worst measured run sits about three orders of magnitude below its ceiling. Budget from the measurement and treat the cap as the backstop) | YES | YES | YES | `dare_care_nomalloc_test` (NX = 2, 4, 8 across the Schur, sign-function, and balanced-Schur variants); `care_convergence_anchor_test` (fixed-seed near-axis and simple-input scale sweeps under all three method tags, plus magnitude-band anchors); `riccati_magnitude_test` (both ends of the arithmetic range); leg 1 witness calls `dare` and `care`. The allocation cell's heap claim is proved by `dare_care_nomalloc_test`; its stack claim is proved by the `-fstack-usage` frames and runtime watermarks in the stack note below, which are host `double` measurements, now corroborated on silicon by the ESP32 stack probe in `examples/embedded/esp32/main/stack_probe.cpp` -- the on-target `float` peaks track the host `double` prediction to within the scalar width, and that probe also establishes that at `float` the ACCEPTANCE CHECK, not the stack, is what bounds the usable state dimension. At `double` the NUCLEO-H753ZI image designs a two-state gain through `dare` on every boot: the gain agrees with the host's within its derived tolerance, and the closed loop it drives is bit-identical across 21 warm resets (`scripts/board_capture_fixtures/evidence.txt` and `timing.txt`, `arm-none-eabi-gcc` 16.2.0; see "On silicon: the NUCLEO-H753ZI leg" below). That board measures no stack figure of its own for the solve. At `float` the ESP32-POE-EA image designs the same two-state gain through `dare` on every boot: the gain agrees with the host's within its derived tolerance, the closed loop's last state lies within its derived bound, and the closed loop is bit-identical across 25 restarts of the timing image, 20 by the chip's enable pin and 5 by a replug (`scripts/board_capture_fixtures/esp32-evidence.txt` and `esp32-timing.txt`, `xtensa-esp-elf-gcc` 15.2.0 under ESP-IDF v6.0.2; see "On silicon, and what it does to the table above") |
| velocity profile construction: `trapezoidal_trajectory::create` / `double_s_trajectory::create` | YES (the profile is returned by value inside a `ctrlpp::expected`; nothing on the path owns storage) | YES (closed form on the trapezoidal path and on every double-S path but one; the cruise-free double-S rise is a bracket halved to exhaustion, bounded by one more than the significand width, so 25 evaluations for `float` and 54 for `double`) | YES | YES | YES (no random source; identical inputs exhaust the bracket at the identical step) | `trajectory_nomalloc_test` (construction is outside the armed window, as for every other type there); leg 1 + `embedded_core_float`, which instantiate both `create` factories |
| trajectory evaluation (polynomial paths, velocity profiles, `cubic_spline`, `smoothing_spline`, `bspline_trajectory`) | YES | YES (closed form; B-spline recursion bounded by compile-time degree) | YES | YES | YES | `trajectory_nomalloc_test` (evaluate cases for cubic/quintic/septic, trapezoidal, double-S, modified sin/trap, cubic_spline, smoothing_spline, bspline_trajectory); leg 1 + `embedded_core_float`. On silicon at `double`, double-S evaluation on the NUCLEO-H753ZI evidence image reports 0 allocations over 201 evaluations behind a passing allocation canary, and its last step's value equals the host golden (`scripts/board_capture_fixtures/evidence.txt`, `arm-none-eabi-gcc` 16.2.0; see "On silicon: the NUCLEO-H753ZI leg" below). At `float`, double-S evaluation on the ESP32-POE-EA evidence image reports 0 allocations over 201 evaluations behind its passing allocation and heap-view canaries, with a clean heap view in the window, and its last step's value within the single-precision bound derived for it (`scripts/board_capture_fixtures/esp32-evidence.txt`, `xtensa-esp-elf-gcc` 15.2.0 under ESP-IDF v6.0.2; see "On silicon, and what it does to the table above") |
| online planners: `online_planner_2nd` / `online_planner_3rd` | YES | YES (closed-form segment logic) | YES | YES | YES | `trajectory_nomalloc_test` (update and sample cases for both planners); leg 1 + `embedded_core_float` |
| trajectory time scaling: `rescale_to` / `can_rescale_to` / `synchronize` (trapezoidal, double-S) | YES (two-pass over a `std::span`, no owning copy, no allocation; the trapezoidal solve's boundary-duration helper returns a two-scalar aggregate by value and owns no storage) | YES (closed form on the trapezoidal and the rest-to-rest double-S paths -- the trapezoidal plateau and valley shapes are each a single quadratic in the cruise velocity's distance from that shape's own boundary, straight-line with no iteration; bracket exhaustion on the nonzero-boundary-velocity double-S path, bounded by one more than the significand width, so 25 evaluations for `float` and 54 for `double`, plus one step per binary exponent the bracket spans) | YES | YES | YES (no random source; identical inputs exhaust the bracket at the identical step) | `trajectory_nomalloc_test`; `trajectory_rescale_anchor_test` and the rescaling cases in `trajectory_hardening_test`; leg 1 + `embedded_core_float` |
| `recursive_arx` / `rls` | YES (the update's result is a `ctrlpp::expected<void, rls_update_error>` holding one enumerator and no owning member; the two norms feeding the resolution floor are unevaluated Eigen expressions over existing storage) | YES (rank-one update, no loop; the refusal guard adds a fixed count of reads per cycle, every dimension a compile-time template parameter -- `NP*NP + NP` for the carried-state scan, `NP` for the regressor scan, and two `NP`-term norms for the denominator's scale) | YES | not covered (the leg 1 witness does not include the sysid headers) | YES (the guard is a pure predicate on the operands and the carried members; no random source, no clock) | `sysid_nomalloc_test` (rls and recursive_arx update cases, re-run green with the guard on every cycle of both 256-cycle armed windows); the typed refusals asserted in `sysid_hardening_test` |
| `mpc` / `nmpc_dynamic` / `mhe` / `nmhe` (runtime-horizon, dynamic solver) | NO (soft real-time: the solve allocates and iterates) | YES when capped (`max_eval`, OSQP `max_iter`) | YES with NLopt's `max_time = 0` (the default) and argmin's default step budget, which has no `max_time`; a wall-clock budget is non-RT | not covered (opt-in OSQP/NLopt/argmin backends sit outside the embedded core witness) | solver-dependent | labeled soft real-time; caps and defaults in `mpc/nlopt_solver.h`, `mpc/argmin_policies.h`, `mpc/osqp_solver.h`, `mpc/argmin_qp_solver.h` |
| `nmpc` (the DEFAULT; = `nmpc_static`, compile-time horizon, argmin `nw_sqp`, bounded decision `NV` + constraint `MaxM`) | YES below a decision dimension `NV` of 48, strict-zero at 0.00 allocs/step in steady state on a converging trajectory (`nw_sqp`) and at an exact regulation equilibrium, where every solve takes the SQP null step (`nw_sqp` and `filter_nw_sqp`, static path), and **NO at and above it**, where the solve crosses the same blocked-Householder boundary the estimator rows do (see "The heap claim's boundary" below). **The heap is not the binding cost here either** (see "Stack cost of the predictive controller row" below, which a hard-real-time caller must read before sizing a task stack: the configuration this cell's own evidence pins needs 21,480 bytes of stack in steady state and 33,432 on its first solve, and constructing it needs 48,040) | YES when capped (`max_eval`) | YES: the default step budget instantiates no clock read; the opt-in `argmin_step_and_time_budget` is non-RT | not covered by this witness (argmin's `-fno-exceptions` instantiation is clean upstream; the ctrlpp-side no-exceptions dogfood is a separate witness) | YES: the constraint bound feeds only the QP result-multiplier storage, never the compute workspace, so argmin's `nw_sqp` bit-identity golden is unchanged | `nmpc_static_nomalloc_test`: the shipped pin (double_integrator NX=2 NU=1 NH=5 → NV=17, MaxM=12; sentinel `eigen_assert` + `EIGEN_RUNTIME_NO_MALLOC` + global `operator new` counter; `static_assert(strict_allocation_free)`), plus a two-configuration bracket over the damped chain that walks the boundary, at `NV` 47 and 48. **That is 2 of the 25 configurations this row publishes a figure for**, and the reason is measured: a unit holding all 25 peaked at 21.9 GB of compiler resident set and was killed by the kernel, and one holding 8 built and ran green at all 8 but cost 17.4 GB to compile and 232 s to run against roughly a second for the pin alone. Its equilibrium case holds the pin's controller at the origin for both `nw_sqp` and `filter_nw_sqp` on the static path, asserts the sentinel and the counter, and fails against argmin revisions before the null-step fix through the sentinel alone. The NUCLEO-H753ZI image, run on the board, reports 0 allocations over 201 armed steps behind a passing allocation canary. Every one of its solves stops on the stationarity test, and its realized closed-loop cost lies within a derived bound of the exact finite-horizon law's (`scripts/board_capture_fixtures/evidence.txt`, image `e03bfd30f81f283f155217320a8aa0f7ee8ea379`, `arm-none-eabi-gcc` 16.2.0). Its stack high-water, cycle counts and compile ceiling on that board are in "On silicon: the NUCLEO-H753ZI leg" below |
| static-memory linear MPC | planned | planned | planned | planned | planned | not implemented; the future hard real-time path (see below) |
| the estimator `update` rejection guard (`kalman_filter`, `ekf`, `ukf`, `mekf`, `manifold_ukf`, `luenberger_observer`, `complementary_filter`) | YES (an `allFinite()` scan is an unevaluated Eigen expression over existing storage; the result is a `ctrlpp::expected<void, E>` holding one enumerator and no owning member) | YES (a fixed count of reads per step: state + covariance + measurement, every dimension a compile-time template parameter -- `NX + NX*NX + NY` for the covariance filters, `NX + NY` for the observer, `4 + NB + NE*NE + NY` for the MEKF, at most `7 + 9 + 1` for the complementary filter) | YES | YES | YES (a pure predicate on the operands; no random source and no state read beyond the members it scans) | `estimation_nomalloc_test` re-run green with the guard on every step of all six converted cases (128-step armed window each, 0 allocations, sentinel clean); the four-part rejection asserted in `{kalman,ekf,ukf,mekf,manifold_ukf,luenberger,complementary_filter}_hardening_test`; `scripts/cross_compile_check.sh` all three legs PASS |
| the controller step rejection guard (`pid::compute` both overloads, `mrac_controller::evaluate`, `l1_controller::evaluate`) | YES (an `allFinite()` scan is an unevaluated Eigen expression over existing storage; the result is a `ctrlpp::expected<vector_t, E>` whose payload is the same by-value vector the surface already returned, so no owning member is added) | YES (a fixed count of reads per cycle, every dimension a compile-time template parameter: `pid` scans only the members its policy composition makes live, at most `11*NY` plus one scalar test on the step; `mrac` scans `NX + NU` vector and `NU*(NX + NU)` matrix entries; `l1` scans `2*NX + 2*NU`. The `l1` cycle additionally scans `NU` more for the pre-projection finiteness test that feeds `health()`) | YES | YES | YES (a pure predicate on the operands and the carried members; no random source, no clock) | `pid_nomalloc_test` re-run green with the guard on every cycle of all three composed cases (256-cycle armed window each, 0 allocations, sentinel clean); the four-part rejection asserted in `{pid,mrac,l1}_hardening_test`, each proven to fail with its guard deleted; both standing trees green at baseline (90/90, 129/129) |

## The heap claim's boundary on the caller-dimensioned rows

The `allocation-free?` cells above read `YES` for the five estimator rows and for
the default `nmpc`. **On five of those six that holds below a boundary and not
above it, and the boundary is at 48.** It was found by arming the no-allocation
guard at every configuration this document prints a figure for, which is what the
grids in `estimation_<row>_nomalloc_test` and `nmpc_static_nomalloc_test` do; the
single-configuration tests that preceded them sat at two states and one output,
far below it, and could not have found it.

**Where the 48 comes from.** The gain solve on `kalman_filter`, `ekf`, `mekf` and
`manifold_ukf` is a column-pivoting Householder QR of the `NY x NY` innovation
covariance applied to a right-hand side with one column per state -- per error
state on `mekf`, three on `manifold_ukf`. Eigen applies a Householder sequence
BY BLOCK once the sequence is at least `BlockSize = 48` long and the destination
has more than one column (`Householder/HouseholderSequence.h`), and the blocked
application declares its block reflector factor as

    Matrix<Scalar, TFactorSize, TFactorSize, RowMajor> T(nbVecs, nbVecs)

with `TFactorSize` taken from the destination's compile-time column count
(`Householder/BlockHouseholder.h`). The destination there is a dynamically sized
block, so `TFactorSize` is `Dynamic` and **`T` is heap allocated: 18,432 bytes,
once per update, through Eigen's own aligned allocator.** The 48 is therefore
read out of the linear-algebra library rather than fitted, and it is the same in
3.4.0 and 3.4.1, which is the one axis this document otherwise has to qualify
every figure against. **The published numbers below are the walk's, not that
arithmetic's:** every row was measured on both sides.

| row | axis the boundary lies on | last allocation-free | first allocating |
|---|---|---:|---:|
| `kalman_filter` | measurement dimension, at two states or more | 47 | 48 |
| `ekf` | the same | 47 | 48 |
| `ukf` | none at the default decomposition; see below | -- | -- |
| `mekf` | measurement dimension, at any bias dimension | 47 | 48 |
| `manifold_ukf` | measurement dimension | 47 | 48 |
| `nmpc_static` | decision dimension `NV` | 47 | 48 |

**Two exemptions, measured rather than inferred from the condition.**

- A single-state `kalman_filter` or `ekf` gives that solve a one-column
  right-hand side, so the blocked path is never entered. Both are measured
  allocation-free at one state against 48 outputs, and `kalman_filter` at one
  state against 128 outputs as well.
- **`ukf` builds no Householder sequence at all**, because its gain
  decomposition defaults to LDLT. Its whole grid is allocation-free, including
  the corner at 128 states and 128 outputs. Selecting the QR option -- a public
  configuration field -- puts it back on the boundary exactly: measured clean at
  47 outputs and allocating at 48.

**THE GLOBAL ALLOCATION COUNTER CANNOT SEE ANY OF IT.** Every one of those
allocations goes through Eigen's `aligned_malloc`, which calls `std::malloc`
directly and never reaches the replaced `operator new`. The counter read **zero
at every allocating configuration**, and each of those points asserts that zero
alongside the fired sentinel, because that is the demonstration: a guard carrying
only the counter -- the guard most projects write -- would report every one of
these configurations allocation-free. The two mechanisms are not redundancy.

**What this means for a caller sizing a system.** Every entry in every
supported-maximum table in this document is at most 36, so **a configuration that
fits a 64 KiB task stack is below the boundary on every row**. Reaching it needs
a task stack much larger than any tabulated here. The boundary is not on the path
of a caller sizing from those tables; it is on the path of a caller with a large
stack who reads `YES` in the matrix and stops there.

**What is armed, and what is not.** 112 configurations across the six rows: 26 on
`kalman_filter`, 25 on `ekf`, 27 on `ukf` (two of them with the QR option), 12 on
`manifold_ukf`, 20 on `mekf` and 2 on `nmpc_static`. Each row's test asserts its
own count rather than describing it, so a point cannot be dropped without the
test failing. Two gaps, stated rather than left to be inferred:

- The five estimator rows arm **every configuration their tables print a figure
  for**, except the value one in the 4 KiB rows of the `ukf` and `manifold_ukf`
  supported-maximum tables: one state and one output, one state at four outputs
  and four states at one output on `ukf`, and one output on `manifold_ukf`. The
  interior fill behind their supported-maximum tables is 4,323 configurations and
  is **not** individually armed: nothing above is a claim about a dimension lying
  between two armed ones.
- The predictive controller row arms **2 of its 25**, both sides of its
  boundary, plus the shipped pin at `NV` 17. Each of its points instantiates the
  whole nonlinear-programming substrate and then solves at that dimension: a unit
  holding all 25 peaked at 21.9 GB of compiler resident set and was killed by the
  kernel, and one holding 8 -- the endpoint of every published line plus the
  bracket -- built and ran green at all 8 but cost 17.4 GB to compile and 232 s
  to run. Neither is a cost this suite can carry on an ordinary machine, so what
  is kept is the pair that walks the boundary. **The published endpoints of this
  row are not armed.** Armed against `g++ (GNU) 16.1.1`,
`-std=c++20 -fno-exceptions -fno-rtti`, `double`, Eigen 3.4.0 -- the release the
test tree fetches, which is not the 3.4.1 the stack figures above were taken
against.

## Stack cost of the Riccati solve

The `allocation-free?` column is a **heap** statement. For most rows in this
table that is the whole resource story, because their hot paths hold a handful
of fixed-size matrices. It is not the whole story for the Riccati solve: it
runs a fixed-size decomposition whose dimension is the caller's `NX`, its frames
grow fast in that dimension, and on the four-to-sixteen-kibibyte task stacks this
milestone's own targets run, **the stack is what will break a caller, not the
heap.** How fast is measured below rather than asserted: the acceptance check's
dominant object is the `M x M` forward-error operator with `M = NX(NX+1)/2`, so
its ASYMPTOTIC growth is the fourth power of `NX`, but over the measured range
the frame grows 17.4-fold for a four-fold increase in `NX`, which sits between
the square and the fourth power because the operator has not yet swamped the
part of the frame that does not scale with it.

Four numbers per state dimension, at input dimension 1, all four taken by
`tools/stack_watermark.sh` in one run. The first is the acceptance check's own
frame as the compiler's `-fstack-usage` report gives it, and the second is the
deepest frame anywhere in the chain with the function that owns it, which is not
the same function at every dimension. The third is the deepest disturbed word
measured at runtime over the whole `dare` chain, which is the number a task
stack must actually cover. The fourth is the same measurement on the solve
ALONE -- the symplectic operand factorization, the symplectic build, and the
real-Schur reorder and Riccati extraction, with the acceptance check and the gain
formation taken off the end -- so the cost of the check is separable from the
cost of the solve.

| `NX`, input dimension held at 1 | acceptance-check frame | deepest frame | function owning it | whole-chain peak, shipped | whole-chain peak, solve alone |
|---:|---:|---:|---|---:|---:|
| 2 | 864 | 2,304 | `swap_real_schur_2x2_general` | 4,552 | 3,368 |
| 4 | 2,400 | 2,400 | `estimate_riccati_forward_error` | 11,176 | 8,176 |
| 6 | 6,432 | 6,432 | `estimate_riccati_forward_error` | 21,912 | 13,560 |
| 8 | 15,056 | 15,056 | `estimate_riccati_forward_error` | **41,304** | 25,416 |

**At two states the acceptance check is NOT the deepest frame in the chain**, and
a reader who took the first column for the chain's maximum would be low by a
factor of 2.7 there. The two coincide from four states up.

Supported maximum state dimension, from the whole-chain peak, strict: no margin
for the caller's own frames and none for RTOS overhead.

| task stack | supported `NX`, shipped, input dimension held at 1 | supported `NX`, solve alone, input dimension held at 1 |
|---|---|---|
| 4 KiB | **none** | `NX <= 2` |
| 8 KiB | `NX <= 2` | `NX <= 4` |
| 16 KiB | `NX <= 4` | `NX <= 6` |
| 32 KiB | `NX <= 6` | `NX <= 8` |
| 48 KiB | `NX <= 8` | `NX <= 8` |
| 64 KiB | `NX <= 8` | `NX <= 8` |

**The acceptance check costs exactly one rung of this ladder on every task stack
from 4 KiB through 32 KiB, and nothing at all above that.** It is not a rounding
effect at the small end: at 4 KiB the shipped chain does not fit at any state
dimension while the solve alone fits two states, and the shipped chain reaches
41,304 bytes at eight states against the solve's own 25,416. A caller who cannot
afford that rung is choosing between a smaller plant and an unverified answer,
and this table is what that choice costs. Above 32 KiB the check is free in
these terms, because both constructions run out of ladder at eight states rather
than out of stack.

**Provenance.** These are HOST measurements: `g++ (GNU) 16.1.1 20260728`,
`-std=c++20 -O2 -fno-exceptions -fno-rtti -pthread`, no `-march` (driver default
`-mtune=generic -march=x86-64`), x86-64 Linux, `double`, Eigen 3.4.1, corpus the
discrete damped chain at **input dimension 1**, runtime watermarks taken on a
pthread with a 64 MiB stack against a zero-byte harness floor and a 1,024-byte
harness gap. A target's own frames differ with its ABI, register file and
calling convention.

**The input dimension and the linear-algebra release are part of that provenance
and not details of it, and both were varied rather than assumed fixed.** Same
instrument, same flags, same corpus:

| what was varied | `NX = 2` | `NX = 4` | `NX = 6` | `NX = 8` |
|---|---:|---:|---:|---:|
| whole-chain peak, three inputs instead of one | +80 | +256 | +448 | +80 |
| whole-chain peak, Eigen 3.4.0 instead of 3.4.1 | -16 | -96 | -176 | -176 |
| acceptance-check frame, Eigen 3.4.0 instead of 3.4.1 | **-48** | -16 | 0 | 0 |

The whole-chain shifts are 0.2% to 2.3%, small enough not to move any entry in
either ladder above and large enough that a figure quoted without its input
dimension cannot be reproduced. The frame shift is 5.6% at two states, and it
matters more than its size suggests: **the 864 above is an Eigen 3.4.1 figure and
the corresponding 3.4.0 figure is 816**, so a reader who measures against the
release the test tree fetches will not reproduce the first cell of the first
column. Re-run `tools/stack_watermark.sh` under your own pairing rather than
reading across.

### On silicon, and what it does to the table above

The host figures above were a PREDICTION. They have now been tested on a board,
and they hold once the one variable that separates the two measurements is
accounted for: **the host table is `double` and the board runs `float`.**

Measured on an Olimex ESP32-POE-EA (ESP32-D0WD-V3, chip revision v3.1, Xtensa LX6
at 160 MHz, single-precision FPU) under the released ESP-IDF v6.0.2,
`xtensa-esp-elf-gcc` 15.2.0 (crosstool-NG esp-15.2.0_20251204) at `-Os` with
exceptions and RTTI off and the framework's default C assertions on, Eigen 3.4.0,
the end-of-stack watchpoint on, one FreeRTOS task per dimension so each figure is
that dimension's own peak rather than a running minimum over a sweep. These
figures were taken while the probe ran in both images the leg builds, which then
differed only in the allocation sensor, so both columns are given; the ratio uses
the image without the sensor. The probe now runs in the allocation-evidence
image only, and the timing image is built without it and without assertions, so
the column without the sensor belongs to that earlier timing build and is not
reproduced by the current one:

| `NX` | on-silicon peak, `float`, no sensor | with the allocation sensor | host peak, `double` | host / board | accepted? |
|---:|---:|---:|---:|---:|:--|
| 2 | 2,804 | 2,700 | 4,552 | 1.62x | solved |
| 3 | 4,120 | 4,032 | -- | -- | refused |
| 4 | 5,948 | 5,924 | 11,176 | 1.88x | solved |
| 5 | 8,344 | 8,368 | -- | -- | refused |
| 6 | 11,500 | 11,492 | 21,912 | 1.91x | refused |
| 8 | 20,552 | 20,560 | 41,304 | 2.01x | refused |

These figures replace earlier ones whose framework version could not be pinned.
They moved by as much as 296 bytes (`NX = 2`, from 3,100) and by at most 40 bytes
elsewhere, and the accepted and refused dimensions did not change. A figure here
belongs to this toolchain, this framework version and this image, not to the
source. The current allocation-evidence image, whose ELF has the sha256 digest
`c7bc0e9faef5707f2ec9d327766ecc397c84312c91872c51383a3e628aa25f98`, still runs
the probe, and it reproduces the column with the sensor to within 8 bytes at
every dimension: 2,692, 4,024, 5,932, 8,360, 11,500 and 20,552 bytes. Its pass at
the control task's own 8,192 bytes accepts `NX = 4` with 2,264 bytes free and
stops at `NX = 5` on the same watchpoint line below, so the supported maximum and
every accepted and refused dimension held.

**The ratio converges on exactly the scalar width.** It is 1.62x at `NX = 2` and
climbs to 2.01x at `NX = 8`, which is what a halved scalar predicts for a chain
whose cost is dominated by `NX`-dimensioned arrays and diluted at small `NX` by
fixed overhead that does not scale with the scalar. The host table is therefore
usable as written for `double` and halves for `float`; it is not refuted and it
is not to be replaced by these figures. **The two columns do not share an input
dimension**: the host column is a single-input pose throughout, while the board
probe drives one input at `NX = 2, 3, 5` and two at `NX = 4, 6, 8`. The ratio is
therefore a statement about the scalar across that pair of pose families and not
a conversion factor to apply cell by cell.

**Supported maximum ON THE STACK for this board: `NX = 4`, leaving 2,248 bytes of
an 8,192-byte task stack (2,256 with the allocation sensor).** That is the stack answer only -- the accuracy gate
refuses `NX = 3` and everything above `NX = 4` at `float`, and clears `NX = 4`
itself by just 3.4 percent, so the usable maximum is lower than this line alone
implies. See the accuracy paragraph below. Established twice and by different
means. A pass giving every
dimension 49,152 bytes recorded `NX = 5` needing 8,344 bytes, which is 152 more
than the task stack has (8,368 and 176 with the sensor); a second pass at the control task's own 8,192 bytes then
walked the same dimensions and the board reported

    Debug exception reason: Stack canary watchpoint triggered (probe_nx5)

so the overflow is OBSERVED AND NAMED at the instruction that caused it, not
inferred from a reset. Both stack-overflow detection modes are enabled in
`examples/embedded/esp32/sdkconfig.defaults`; the end-of-stack watchpoint is what
fired, and it is the more precise of the two because the canary check only runs
at a context switch. Enabling the watchpoint costs up to 60 bytes of every task's
usable stack, so these figures and any taken without it are not interchangeable.

**THE STACK IS NOT THE BINDING LIMIT FOR `float`. THE ACCEPTANCE CHECK IS, AND
IT IS NOT A BOARD PROPERTY.** Only `NX = 2` and `NX = 4` were accepted; every
other probed dimension was refused, including `NX = 6`, whose 11,500 bytes would
fit a 16 KiB task comfortably. A caller sizing a 16 KiB task for a six-state
`float` plant would find the stack sufficient and the answer refused.

This reproduces exactly on the host -- same dimensions accepted, same refused,
with `double` answering every one of them -- so it is a property of the SCALAR
TYPE and the library, not of Xtensa. It is pinned by
`tests/unit/dare_float_precision_test.cpp`, which asserts the bands rather than
this paragraph.

**The refusing site is the forward-error gate specifically**, not the
definiteness floor and not the gain: measured on the solver's own answer, the
definiteness pivots and the gain are healthy at every dimension, and the
estimated relative forward error as a fraction of the `float` half-significand
margin `sqrt(eps) = 3.4527e-04` is

| `NX` | 2 | 3 | 4 | 5 | 6 | 8 |
|---|---:|---:|---:|---:|---:|---:|
| estimate / margin | 0.169 | 5.172 | **0.966** | 4.625 | 2.590 | 4.105 |

**THERE IS NO SUPPORTED MAXIMUM STATE DIMENSION HERE, BECAUSE THE STATE DIMENSION
IS NOT WHAT DECIDES IT.** Widening the sweep to vary the input count
independently shows the pattern is not about `NX` at all. Measured true relative
error of the `float` answer against a `double` reference, with
`group = NX / NU`, the number of states each input must reach through:

| pose | `group` | true error | / margin | verdict |
|---|---:|---:|---:|---|
| `NX=2, NU=1` | 2 | 5.8e-05 | 0.169 | accepted |
| `NX=8, NU=4` | 2 | 2.5e-04 | 0.725 | accepted |
| `NX=4, NU=2` | 2 | 3.3e-04 | 0.966 | accepted |
| `NX=12, NU=6` | 2 | 3.9e-04 | 1.122 | refused |
| `NX=6, NU=3` | 2 | 4.4e-04 | 1.274 | refused |
| `NX=6, NU=2` | 3 | 8.9e-04 | 2.590 | refused |
| `NX=8, NU=2` | 4 | 1.4e-03 | 4.105 | refused |
| `NX=5, NU=1` | 5 | 1.6e-03 | 4.625 | refused |

**`NX = 8` with four inputs is ACCEPTED while `NX = 6` with three is refused.**
Two things drive the error and neither is the state dimension on its own. It
rises from `NX = 2` to `NX = 4` and then PLATEAUS near `3e-04` whether `NX` is 4,
8 or 12 -- accumulated rounding through the Schur-and-extraction chain, which
saturates. On top of that, chain length per input costs real accuracy: `group = 2`
sits near `3e-04`, `group = 3` near `9e-04`, `group = 4` to `5` near `1.5e-03`,
because at `dt = 0.01` the far end of a long chain is weakly controlled.

**The margin is `sqrt(eps) = 3.4527e-04`, and the `group = 2` family sits at 0.7x
to 1.3x of it.** The error and the criterion coincide, which is why the accept
and refuse outcomes look erratic across that family: those poses straddle the
line and which side they land on is rounding noise. `sqrt(eps)` is the
half-significand criterion -- float has 24 bits and the gate asks that 12 survive
-- and this family genuinely spends about half of float's significand from
`NX = 4` up.

**THE GATE IS NOT OVER-REJECTING.** The estimated error matches the true error
against a `double` reference to three significant figures at every pose above,
`est / true = 1.0` throughout. The refused answers really do carry more than
half-significand error; refusing them is the gate working, not misfiring.

**EVERY FIGURE IN THE TABLE ABOVE IS TOOLCHAIN-SPECIFIC, BY MORE THAN A FACTOR OF
TWO.** Measured on g++ 16.1.1, x86-64, `-O2`. A test that asserted the
`NX = 8, NU = 2` pose exceeds the margin by more than a factor of two -- 4.105
here -- FAILED on Apple clang, where the same expression fell below 2.0. So the
`float` column is a reading of one compiler, not a property of the library, and
no individual accept-or-refuse verdict in it is portable. What IS portable, and
what `tests/unit/dare_float_precision_test.cpp` asserts instead, is the
comparison between scalars on the same pose: `float` spends about 55,000 times
the fraction of its own margin that `double` spends of its own. Treat the table
as an illustration of the mechanism and measure your own toolchain.

What a caller should take from this is a MEASUREMENT INSTRUCTION rather than a
dimension: at `float`, on a plant whose inputs must act through chains of more
than two states, expect the solver to refuse, and expect poses near `group = 2`
to sit on the margin either way. `double` answers every pose in the table.

Witness: `examples/embedded/esp32/main/stack_probe.cpp` for the probe, with each
dimension's solve in `examples/embedded/esp32/main/stack_probe_dimension_impl.h`, started from
`examples/embedded/esp32/main/app_main.cpp`, whose control loop still
designs its gain, runs its 201 steps, streams them over UART2 and reports
`golden diff PASS` in the same boot that carries the probe.

**The families, the allocation sensor and the cycle counter on this board.** The
leg builds two images from the same sources, the stack probe's aside, which only
the evidence image compiles, and every figure below names the one it came from:

- the allocation-evidence image, `ctrlpp_esp32`, whose ELF has the sha256 digest
  `c7bc0e9faef5707f2ec9d327766ecc397c84312c91872c51383a3e628aa25f98`. It keeps
  the framework's default assertions and task watchdog on, as the demonstration
  ships, and runs the stack probe after its report. Its boot report is committed
  as `scripts/board_capture_fixtures/esp32-evidence.txt`.
- the timing image, `ctrlpp_esp32_timing`, digest
  `a7d830451598f39f56a1b579900c59613d420dc5f89e872a18d0e2fcbdb7fbd5`, built
  without the allocation sensor and without the probe, with assertions off, and
  with the task watchdog released before it reports, so neither a check nor a
  watchdog interrupt lands inside a timed region. Its boot report is committed as
  `scripts/board_capture_fixtures/esp32-timing.txt`.

Both were built with ESP-IDF v6.0.2 and `xtensa-esp-elf-gcc` 15.2.0
(crosstool-NG esp-15.2.0_20251204) at `-Os`, and both report running at 160 MHz,
on core 0 at priority 5 under a 100 Hz FreeRTOS tick, with code in flash behind
the cache and interrupts left enabled. Both were captured with
`scripts/board_capture.sh --board=esp32`, which applies the NUCLEO-H753ZI leg's
refusals to this stream, and its `--self-test` replays those refusals over both
reports on every change. **Every figure in this part belongs to that toolchain,
that framework version, `float` and that posture. None of it converts to or
from the NUCLEO-H753ZI's figures below, which are `double`, a different
compiler, a bare superloop and caches off.**

**Every family passes.** The leg runs four of the NUCLEO-H753ZI's five families
at `float`, against the same host goldens and against bounds derived for single
precision:

| family | quantity judged | departure | derived bound | departure / bound |
|---|---|---:|---:|---:|
| control | final state norm | 1.183e-09 | 3.370e-06 | 3.5e-04 |
| estimation | the covariance entry `P(1,1)`, the velocity variance, after the last step | 2.764e-07 | 5.976e-03 | 4.6e-05 |
| signal processing | the last step's value | 2.569e-08 | 8.694e-03 | 3.0e-06 |
| trajectory | the last step's value | 1.589e-07 | 6.182e-03 | 2.6e-05 |

The signal-processing row and the `biquad` timing row below, like the
NUCLEO-H753ZI's, were captured from images of the earlier `biquad` section
formation, and the row is judged against that formation's golden and bound. The current formation's single-precision bound is
`6.718e-03`, and no image of it has been captured on this board.

The designed gain departs from the host's by `[8.943e-04, 2.766e-04]` against a
tolerance of `2.292e-02`, which at 3.9 percent of its bound is the tightest
margin on the board. Both images print the same verdicts. The predictive family
is declared absent on every report: the backend's solver computes in `double`
only, and its first solve's measured stack on the NUCLEO-H753ZI, 29,724 bytes,
exceeds this leg's 8,192-byte control task. These figures replace those of the
first capture on this leg, taken under an unreleased framework, on a different
ESP32 module and with an older library: the gain departure was then
`[2.157e-03, 7.549e-04]`, and the final-norm departure was 2.249e-08, judged
against a hand-set tolerance of 2.0e-03 rather than a derived bound. Both fell,
by 59 and 63 percent for the gain and by 95 percent for the final norm. Since
the library and the module changed with the framework, no part of the move can
be attributed to the framework alone.

**Allocation, on the evidence image.** Two views are printed side by side,
because neither sees everything. The first is a counter on the replaced C++
`operator new`, plain and aligned, beside the linear-algebra library's
no-allocation sentinel. Its canary allocates once through each `new` and once
through an Eigen heap vector, and reads `observed=2 expected=2 eigen_trips=1`.
This leg does not wrap the C allocation family, so a C allocation is visible only
to the second view, the free heap. That view's own start call allocates a 24-byte
snapshot inside every window, so it can never read zero: the image measures that
self-count on an empty window, and the view's canary must see the 32 bytes it
requests beyond it (`retained=100 requested=32 self=24`). Behind both canaries
every family reports 0 allocations and 0 sentinel trips over its 201 steps, and
every armed window's heap view, the closed loop's 201 windows included, reads
clean against the 24-byte self-count. The heap view counts the whole heap, so a
clean window also says no other task allocated during it; if one had, the view
could not say whose the allocation was.

**Stack, on both images.** The control task's high-water, read once after the
report's last step under the same watchpoint setting as the probe: 3,160 bytes
used of 8,192 on the evidence image, and 3,576 on the timing image, which runs
its timed regions on the same task before the report.

**Timing, on the timing image.** Cycles from the core's own cycle counter, which
runs from reset with no enable step. The image proves it live by measurement
before printing any figure (`cycle_counter PASS state=live`), and the timed task
is pinned to one core because each core has its own counter. Every figure is raw:
the sampling overhead is published as its own region and never subtracted. The
cold figure is the first execution after reset (`n = 1`); the warm figures cover
201 repeats, and the last column counts the warm samples during which the
scheduler tick advanced:

| region | block | cold | warm min | warm median | warm max | warm max per iteration | at 160 MHz | tick crossed |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| sampling overhead | 1 | 24 | 24 | 24 | 24 | 24 | 0.15 us | 0 of 201 |
| control: gain design (`lqr_gain`) | 1 | 951,561 | 506,512 | 509,052 | 614,744 | 614,744 | 3.84 ms | 65 of 201 |
| control: closed-loop step | 201 | 48,021 | 41,030 | 41,030 | 42,301 | 210.5 | 1.32 us | 5 of 201 |
| estimation: `kalman_filter` step | 201 | 1,033,128 | 953,131 | 954,382 | 964,545 | 4,798.7 | 29.99 us | 120 of 201 |
| dsp: `biquad` sample | 201 | 10,017 | 7,476 | 7,476 | 9,375 | 46.6 | 0.29 us | 1 of 201 |
| trajectory: double-S evaluation | 201 | 41,602 | 26,400 | 26,400 | 28,299 | 140.8 | 0.88 us | 4 of 201 |

Interrupts and the scheduler are left as shipped and nothing is masked, so the
warm band carries what a task on this core sees: the tick interrupt, flash-cache
misses and the other core's flash traffic, which the image does not separate.
Outside the overhead region, the cold figure runs from 1.08 times the warm
median (estimation) to 1.87 times (gain design). The evidence image, the demonstration as shipped, keeps the
assertions these figures run without. An earlier timing build that kept them,
also kept the probe and left the task watchdog armed measured warm medians up to
4.8 times higher, and the watchdog's interrupt landed inside an estimation
sample at 2.3 times that region's median. It differed in all three ways at once,
and its capture does not separate them, so its figures are retired rather than
attributed.

**Determinism on this board.** The control family's closed loop is printed as
201 rows of binary32 bit patterns. Both series ran on the timing image whose ELF
has the sha256 digest
`a7d830451598f39f56a1b579900c59613d420dc5f89e872a18d0e2fcbdb7fbd5` (its
allocation-evidence twin is
`c7bc0e9faef5707f2ec9d327766ecc397c84312c91872c51383a3e628aa25f98`), built with
assertions off and with the task watchdog released before it reports, and
compared with `scripts/board_capture.sh compare` at a bound of 0 cycles stated
before collection:

- 20 runs, each restarted by a pulse on the chip's enable pin through the console
  bridge, with no operator action;
- 5 runs, each after the operator unplugged and replugged the board's own USB
  cable. Two of them were accepted on a second replug, because the first
  attempt's report had lost its head and was refused.

Every run in both series decoded `reset=power_on reset_reason=1 reset_rom=1`.
The block was byte-identical in all 25 runs (sha256
`512113e3b1844698fca4ac90bc779044081be75e19d073b016d50ef7916ae8df`), the same
block the leg's earlier images printed with the probe and assertions on. Every
timing figure in every region had zero spread within each series, and the two
series measured the same figures.

**WHAT THE LABEL CANNOT SAY.** This chip reports a power cycle and an
enable-pin pulse with the same reset reason, so on this board, unlike the
NUCLEO-H753ZI below, the decoded label does not separate the two series. The
enable-pin series is warm-restart evidence only, whatever its label reads; its
summary says so, and the capture refuses to pool it with power-on runs. The
power-on series rests on the operator's replug, and the only check the capture
makes of it is that the console bridge's device node vanished and returned,
which shows that the USB connection dropped, not that the supply decayed. Zero
differences leaves a per-run difference rate up to about 15 percent plausible at
95 percent confidence for the 20 warm runs, and up to about 60 percent for the 5
power-on runs.

### On silicon: the NUCLEO-H753ZI leg

A second board runs five of this table's families at `double`, judges each
against a host golden within a derived bound, and reports what it measured
about itself while doing so. Measured on a NUCLEO-H753ZI (STM32H753ZI, Cortex-M7
with a double-precision FPU), built with `arm-none-eabi-gcc (Arch Repository)
16.2.0` at `-Os` with exceptions and RTTI off, newlib-nano, Eigen 3.4.0 with
`EIGEN_STACK_ALLOCATION_LIMIT=8192`, and the backend at the pinned
`30a2cd4ea2fdab5b548080644908db640542eacd`. The image reads its posture off the
hardware and prints it: the core at 64 MHz from the internal oscillator with no
PLL, both L1 caches off, 7 flash wait states, code in flash, data and stack in
DTCM, the heap in AXI SRAM, no interrupt enabled, one bare superloop.

**EVERY FIGURE IN THIS SUBSECTION BELONGS TO THAT TOOLCHAIN, THAT FIXED-SIZE
ALLOCATION LIMIT, THAT SCALAR AND THAT POSTURE.** None of them carries over to a
different compiler or compiler version, a different
`EIGEN_STACK_ALLOCATION_LIMIT`, `float`, or a different clock and cache setup.
With the instruction cache off and 7 wait states, code placement alone moves the
cycle counts, so a timing figure belongs to the image it was read from and not
to the source. The ESP32 figures above were taken on `xtensa-esp-elf-gcc` 15.2.0
under ESP-IDF v6.0.2 at `float`; the two legs are not interchangeable.

| family | what the image runs | row it bears on |
|---|---|---|
| control | `lqr_gain`, through `dare`, on a two-state, one-input plant, then 201 closed-loop steps of the designed gain | Riccati steady-state solve |
| estimation | `kalman_filter` at two states, one input and one output, 201 predict-and-update steps | `kalman_filter` |
| signal processing | `biquad::low_pass`, 201 samples | dsp |
| trajectory | `double_s_trajectory::create`, then 201 evaluations | trajectory evaluation |
| predictive | the default `nmpc` at the shipped pin (`NX = 2`, `NU = 1`, `NH = 5`, so `NV = 17` and `MaxM = 12`), 20 warm-up solves, then 201 solve-and-plant steps | `nmpc` |

Two images are linked from the same sources, and every figure below names the
one it came from:

- the allocation-evidence image, `ctrlpp_nucleo_h753zi`, build id
  `e03bfd30f81f283f155217320a8aa0f7ee8ea379`, which carries an allocation
  sensor. Its boot report is committed as
  `scripts/board_capture_fixtures/evidence.txt`.
- the timing image, `ctrlpp_nucleo_h753zi_timing`, build id
  `5e0486618196df29e7e73fc391fd1103515f6d6a`, linked without the sensor so the
  timed code carries none of its checks. Its boot report is committed as
  `scripts/board_capture_fixtures/timing.txt`.

Both reports were captured with `scripts/board_capture.sh`. It refuses a report
whose fingerprint is not the built image's, whose canary did not pass, or whose
stack, timing or verdict lines are incomplete, and its `--self-test` replays
every one of those refusals over the two committed reports on every change.

**Every family passes.** The departure from the host golden, against the bound
derived for it:

| family | quantity judged | departure | derived bound |
|---|---|---:|---:|
| control | final state norm | 2.062e-17 | 6.070e-15 |
| estimation | the covariance entry `P(1,1)`, the velocity variance, after the last step | 6.661e-16 | 2.226e-11 |
| signal processing | the last step's value | 0 | 3.239e-11 |
| trajectory | the last step's value | 0 | 2.303e-11 |
| predictive | closed-loop cost over the run | 1.425e-04 | 1.134 |

**The signal-processing row and the `biquad` timing row below belong to the
section formation both images were built from**, which rounds 22 times per step:
its golden is `-0.86508002466156553` and its bound `3.239e-11`. The current
formation builds the low-pass section at unit scale and rounds 17 times per step;
its host golden is `-0.86508002466156531` and its `double` bound `2.503e-11`. No
image of the current formation has been captured on this board.

The designed gain departs from the host's by `[-1.506e-12, -4.006e-13]` against
a tolerance of `8.538e-11`. The predictive bound holds only under its own
premise, which the report checks: all 221 of the run's solves stopped on the
stationarity test.

**Allocation, on the evidence image.** A canary first proves the sensor can see
an allocation: it makes one deliberate allocation through each path the sensor
watches (a C allocation, a C++ `new` and an Eigen heap vector) and reads
`observed=3 expected=3`. Behind that canary every family reports 0 allocations
and 0 sentinel trips over its 201 steps. The canary does not extend the sensor's
sight to the C library's reentrant allocator entry points, which never pass
through the wrapped symbols. Setup heap: 1,088 bytes before the predictive
family is built, then 76,720 bytes after its construction and its 20 warm-up
solves. That is exactly the heap the image reserves: the reserve is set to the
measured requirement, and a request past it is refused by name on the console.

**Stack, on the evidence image.** The watermark paints the stack below each
region's origin, down to the end of bss, 73,012 bytes below the stack top, and
walks it after the region. Every line reports a zero harness floor and an
unsaturated window:

| region | high-water below the stack top, bytes |
|---|---:|
| startup: reset, static constructors, entry to `main` | 468 |
| predictive construction | 2,752 |
| predictive first solve | 30,616 |
| predictive steady state, the 201 window solves | 26,588 |
| whole program | 30,616 |

The whole program's high-water is the first solve, and the image's link-time
stack reserve, `_Min_Stack_Size`, is set to exactly that figure. The capture
refuses any report whose high-water exceeds the reserve it reads out of the
image. **The construction figure is not the host table's construction figure.**
The board builds the controller in static storage, so 2,752 bytes is the
constructor chain alone; the host table builds the controller as a local and so
counts the object and its inline solver. The first-solve and steady-state
figures measure the corresponding chains on a different toolchain and
instruction set; they do not replace the host's 33,432 bytes on the first solve
and 21,480 in steady state, which "Stack cost of the predictive controller row"
below gives with the backend revision they were measured against.

**Compile ceiling.** At `EIGEN_STACK_ALLOCATION_LIMIT=8192` and `double`, with
the leg's own flags, the largest decision dimension of the predictive row that
compiles is `NV = 32` (`NX = 2, NH = 10` and `NX = 10, NH = 2`, one input), and
`NV = 33` does not (`NX = 1, NH = 16` and `NX = 16, NH = 1`). The first object
to trip is an `NV x NV` fixed-size matrix in the backend's SQP state buffers,
built on the first solve, and the library's assertion reads `8712 <= 8192`.
This is the mechanism behind the host ceiling of `NV = 128` in "The largest
configuration that COMPILES" below, at a limit sixteen times smaller. No stack
figure on this board covers a configuration other than the shipped pin.

**Timing, on the timing image.** Cycles from the core's DWT cycle counter, which
the image proves live before printing any figure. Every figure is raw: the
sampling overhead is published as its own region and never subtracted. The cold
figure is the first execution after reset (`n = 1`); the warm figures cover 201
repeats. A region of block 201 times 201 consecutive steps, and its
per-iteration figure divides by that block:

| region | block | cold | warm min | warm median | warm max | warm max per iteration | at 64 MHz |
|---|---:|---:|---:|---:|---:|---:|---:|
| sampling overhead | 1 | 114 | 93 | 93 | 114 | 114 | 1.78 us |
| control: gain design (`lqr_gain`) | 1 | 319,269 | 319,137 | 319,245 | 319,245 | 319,245 | 4.99 ms |
| control: closed-loop step | 201 | 80,450 | 80,344 | 80,344 | 80,385 | 399.9 | 6.25 us |
| estimation: `kalman_filter` step | 201 | 2,382,878 | 2,382,878 | 2,382,918 | 2,382,918 | 11,855.3 | 185.2 us |
| dsp: `biquad` sample | 201 | 11,946 | 11,881 | 11,881 | 11,946 | 59.4 | 0.93 us |
| trajectory: double-S evaluation | 201 | 31,452 | 31,417 | 31,417 | 31,452 | 156.5 | 2.44 us |
| predictive: one solve and plant step | 1 | 141,964,804 | 4,862,868 | 4,882,672 | 9,676,784 | 9,676,784 | 151.2 ms |

In the cheap regions cold and warm differ by at most 132 cycles, and with both
caches off none of that is a cache effect. The predictive band is solver work:
the cold solve ran 30 iterations against 1 to 2 for a window solve. At this
posture the worst window solve, 151.2 ms, is longer than the demonstration's
100 ms sample period, the median of 76.3 ms is shorter, and the cold first solve
takes 2.22 s. A faster posture is not measured here.

**Determinism across warm resets.** The control family's closed loop is printed
as 201 rows of binary64 bit patterns, so a one-bit difference anywhere in the
trajectory is a difference in the report. Over 21 runs of the timing image (a
`scripts/board_capture.sh compare` series of 20 runs restarted with
`st-flash reset`, plus the mass-storage flash that loaded the image), every run
decoded its own reset as a pin reset, the block was byte-identical in all 21,
and it is identical to the evidence image's (sha256
`6ec241518d693cfe6c411d54573478f2a8508e53b63cc08841ee601876e92702`). Every
timing figure in every region had zero spread across the same 21 runs, so the
cold figures above are not one boot's luck. Zero differences in 21 runs leaves a
per-run difference rate up to about 14 percent plausible at 95 percent
confidence; the count is part of the result.

**Determinism across power cycles.** Both images wait at boot, before printing
anything, for at most 30 s, for a byte the capture sends once its reader is
attached, so a report read after a replug keeps its head, and a report whose
wait expired is refused. The power-on series ran on the timing image above. The
operator unplugged and replugged the board's USB cable once per run. All 5 runs
decoded their own reset as a power-on reset,
`reset=power_on reset_flags=cpu+d1+d2+bor+pin+por reset_rsr=0x00FE0000`, no run
decoded another class, every run was accepted on its first attempt, and the
image waited 75 to 134 ms for the handshake byte. The block was byte-identical in
all 5 runs, and it is the same block as the warm series above. Every timing
figure in every region had zero spread across the 5 runs, against a bound of 0
cycles stated before collection, and every figure equals the one in the table
above. Zero differences in 5 runs leaves a per-run difference rate up to about
60 percent plausible at 95 percent confidence, so the power-on half of this
claim is much weaker than the warm half.

## Stack cost of the estimator rows

The `allocation-free?` column is a heap statement for every row in this document,
not only the Riccati one. Six further rows run a fixed-size decomposition whose
dimension the CALLER chooses, so the same question the Riccati row was forced to
answer applies to them: `kalman_filter` and `ekf` factorize the `NY x NY`
innovation covariance, `ukf` and `manifold_ukf` do that and additionally hold a
sigma-point set live across the propagation, `mekf` does the same at the
error-state dimension, and static `nmpc` is bounded by its caller-controlled
decision and constraint dimensions.

**Two dimensions, not one, and they do not have the same names on every row.**
`kalman_filter`, `ekf` and `ukf` take a state dimension and a measurement
dimension. `mekf` takes a BIAS dimension `NB` and a measurement dimension, and
its error state is the DERIVED `NE = 3 + NB`, so the grid runs over what a caller
picks and reports what that induces. `manifold_ukf` takes a measurement dimension
and nothing else: its state is a rotation, so its state dimension is three by
construction and is not a caller's to choose. Four rows below therefore carry a
two-dimensional measurement and the fifth carries a one-dimensional one, and the
fifth says so rather than presenting five readings of one configuration as a
trend.

### What each column is, and what the tables can and cannot resolve

Every row below carries a pair of tables. The first is the measurement:

- **deepest single frame** -- the largest frame the compiler's `-fstack-usage`
  report attributes to any one function in the chain, with the function named.
  This is a LOWER BOUND on what a task stack must hold, because the report
  attributes nothing to callees.
- **whole-chain peak** -- the deepest disturbed word measured at runtime over
  `predict` followed by `update`, painted as one block. This is the number a task
  stack must actually cover. The construction and one warm-up step happen outside
  the painted window, matching the definition the `allocation-free?` column uses.
- The peak is given on three lines through the grid: both dimensions equal, the
  measurement dimension swept with the other held, and the other swept with the
  measurement dimension held. **Every column heading names what it holds.**

The second table is the supported maximum by task stack, and it is computed from
the whole-chain peak of an INTERIOR FILL rather than from those three lines. The
fill measures **every value of one axis against every value of the other** at a
**stride of four**, from the smallest instantiable value to the largest, so an
entry below names the largest measured value that fits. **The stride is four and
nothing is measured between two adjacent values of it**, so an entry is not a
claim about the three integers above it. It is never computed from the frame
column: that column excludes every library frame beneath the named function and
would overstate what fits.

Two limits are published per row, and **they answer different questions**. The
first is the largest dimension that COMPILES AT ALL: past it the instantiation
does not exist and no task stack changes that. The second is the largest whose
whole-chain peak FITS a given task stack: past it the instantiation exists and
the task overflows. A caller who does not fit a dimension needs to know which of
the two it is, because a bigger stack answers one of them and nothing answers the
other.

The harness floor read zero on every one of the sixty-eight measurements behind
the three-line tables and on every one of the 4,323 behind the fill, printed
beside each figure by the instrument, so every peak is attributable to the call
rather than to the harness. The three-line figures reproduced byte-identically
across two independent runs of the whole grid, and **the fill reproduces all 44
of the published cells it covers, byte for byte**, having been taken at a
different painted window and a different build parallelism.

**One caveat at the very bottom of the grid.** The measurement has a resolution
floor: a chain shallower than the harness gap of 1,024 bytes disturbs nothing the
walk can see and reports zero, which is the value the harness floor reports.
Two points in the fill are there: `kalman_filter` and `ekf` at one state and
one output. Re-measured at gaps of 512, 256, 128 and 64 bytes, with a zero floor
at every one, they resolve to **720 and 880 bytes** and are byte-identical across
all four gaps. `ukf` at one state and one output needs 1,040 bytes, which the
fill resolves at the standard gap. No other point in the grid is inside the
gap.

### `kalman_filter`

The frame column belongs to the equal-dimension build.

| dimension | deepest single frame | function owning it | peak, `NX` = `NY`, nothing held | peak, `NY` swept, `NX` held at 4 | peak, `NX` swept, `NY` held at 4 |
|---:|---:|---|---:|---:|---:|
| 2 | 1,504 | `Eigen::internal::make_block_householder_triangular_factor` | 2,056 | 2,376 | 2,424 |
| 4 | 1,504 | `Eigen::internal::make_block_householder_triangular_factor` | 3,416 | 3,416 | 3,416 |
| 8 | 6,624 | `kalman_filter::update_covariance` | 11,928 | 6,264 | 6,792 |
| 16 | 25,088 | `kalman_filter::update_covariance` | 33,128 | 13,144 | 20,184 |
| 32 | 98,848 | `kalman_filter::update_covariance` | **158,712** | 41,688 | 71,224 |

Strict: no margin for the caller's own frames and none for RTOS overhead.

From the interior fill, at a stride of four. **The largest dimension that
compiles at all is 128 on both axes**, and the entries below are far inside it.

| task stack | largest fitting value, `NX` = `NY`, nothing held | largest fitting `NY`, `NX` held at 4 | largest fitting `NX`, `NY` held at 4 |
|---|---|---|---|
| 4 KiB | 4 | 4 | 4 |
| 8 KiB | 4 | 8 | 8 |
| 16 KiB | 8 | 16 | 12 |
| 32 KiB | 12 | 24 | 16 |
| 48 KiB | 16 | 32 | 24 |
| 64 KiB | 16 | 36 | 28 |

### `ekf`

The frame column belongs to the equal-dimension build.

| dimension | deepest single frame | function owning it | peak, `NX` = `NY`, nothing held | peak, `NY` swept, `NX` held at 4 | peak, `NX` swept, `NY` held at 4 |
|---:|---:|---|---:|---:|---:|
| 2 | 1,504 | `Eigen::internal::make_block_householder_triangular_factor` | 1,832 | 2,664 | 2,408 |
| 4 | 1,568 | `ekf::update` | 3,368 | 3,368 | 3,368 |
| 8 | 10,704 | `ekf::update` | 15,688 | 7,000 | 9,048 |
| 16 | 49,968 | `ekf::update` | 62,952 | 23,240 | 26,440 |
| 32 | 201,392 | `ekf::update` | **245,752** | 64,088 | 85,672 |

Strict: no margin for the caller's own frames and none for RTOS overhead.

From the interior fill, at a stride of four. **The largest dimension that
compiles at all is 128 on both axes**, and the entries below are far inside it.

| task stack | largest fitting value, `NX` = `NY`, nothing held | largest fitting `NY`, `NX` held at 4 | largest fitting `NX`, `NY` held at 4 |
|---|---|---|---|
| 4 KiB | 4 | 4 | 4 |
| 8 KiB | 4 | 8 | 4 |
| 16 KiB | 8 | 12 | 8 |
| 32 KiB | 12 | 20 | 16 |
| 48 KiB | 12 | 24 | 20 |
| 64 KiB | 16 | 32 | 24 |

### `ukf`

The frame column belongs to the equal-dimension build. The sigma-point set is
`2 * NX + 1` points held live across the propagation, so this row pays for the
state dimension twice: once in the covariance it factors and once in the set.

| dimension | deepest single frame | function owning it | peak, `NX` = `NY`, nothing held | peak, `NY` swept, `NX` held at 4 | peak, `NX` swept, `NY` held at 4 |
|---:|---:|---|---:|---:|---:|
| 2 | 1,504 | `Eigen::internal::make_block_householder_triangular_factor` | 2,376 | 3,192 | 3,144 |
| 4 | 2,896 | `ukf::update` | 4,280 | 4,280 | 4,280 |
| 8 | 9,872 | `ukf::update` | 13,376 | 7,832 | 10,416 |
| 16 | 38,768 | `ukf::update` | 50,384 | 20,520 | 33,104 |
| 32 | 150,544 | `ukf::update` | **227,320** | 60,280 | 150,888 |

Strict: no margin for the caller's own frames and none for RTOS overhead.

From the interior fill, at a stride of four. **The largest dimension that
compiles at all is 128 on both axes** -- the sigma-point set is a `std::array` of
vectors rather than one fixed-size object, so the allocation limit does not bind
on it and this row's limit is the same as the two above.

| task stack | largest fitting value, `NX` = `NY`, nothing held | largest fitting `NY`, `NX` held at 4 | largest fitting `NX`, `NY` held at 4 |
|---|---|---|---|
| 4 KiB | 1 | 1 | 1 |
| 8 KiB | 4 | 8 | 4 |
| 16 KiB | 8 | 12 | 8 |
| 32 KiB | 12 | 20 | 12 |
| 48 KiB | 12 | 28 | 16 |
| 64 KiB | 16 | 32 | 20 |

The 4 KiB entries are 1 on all three lines: four states and four outputs need
4,280 bytes against 4,096, and the next value the fill measures below four, one
state and one output, needs 1,040. The sigma-point set carries the `NX x NX`
factor it was generated from, so that after a repair the update reduces the
repaired covariance `factor * factor^T` rather than the carried one, and that
factor is part of every figure in both tables.

### `manifold_ukf`

**This row carries ONE caller-controlled dimension.** Its state is a rotation, so
the covariance it propagates is 3 x 3 whatever the caller does, and there is no
second axis to sweep or hold. Every figure below is taken with the rotation state
at its structural three and the measurement dimension swept.

| `NY` | deepest single frame | function owning it | whole-chain peak, rotation state held at 3 |
|---:|---:|---|---:|
| 2 | 1,792 | `manifold_ukf::update` | 3,312 |
| 4 | 2,400 | `manifold_ukf::update` | 4,128 |
| 8 | 3,440 | `manifold_ukf::update` | 5,776 |
| 16 | 13,584 | `manifold_ukf::update` | 18,736 |
| 32 | 31,328 | `manifold_ukf::update` | **49,760** |

Strict: no margin for the caller's own frames and none for RTOS overhead.

From the interior fill, at a stride of four. **The largest `NY` that compiles at
all is 128.**

| task stack | largest fitting `NY`, rotation state held at 3 |
|---|---|
| 4 KiB | 1 |
| 8 KiB | 8 |
| 16 KiB | 12 |
| 32 KiB | 20 |
| 48 KiB | 28 |
| 64 KiB | 36 |

### `mekf`

**The caller picks the BIAS dimension `NB`; the error state `NE = 3 + NB` is
derived from it** and is what the covariance recursion is sized by. The tables
grid over `NB` and report `NE` beside it, because a grid over `NE` would describe
configurations no caller can reach. `NB < 3` does not instantiate at all -- the
propagation subtracts the leading three bias elements from the gyro rate -- so
the ladder's first rung exists on the measurement axis and not on the bias axis.
The frame column belongs to the equal-dimension build.

| `NB` | `NE` | deepest single frame | function owning it | peak, `NB` = `NY`, nothing held | peak, `NY` swept, `NB` held at 4 | peak, `NB` swept, `NY` held at 4 |
|---:|---:|---:|---|---:|---:|---:|
| 2 | 5 | -- | does not instantiate | does not instantiate | 4,840 | does not instantiate |
| 4 | 7 | 4,064 | `mekf::update` | 6,776 | 6,776 | 6,776 |
| 8 | 11 | 21,392 | `mekf::update` | 29,272 | 12,760 | 20,376 |
| 16 | 19 | 68,816 | `mekf::update` | 84,904 | 25,112 | 56,104 |
| 32 | 35 | 262,784 | `mekf::update` | **310,840** | 74,200 | 152,648 |

Strict: no margin for the caller's own frames and none for RTOS overhead.

From the interior fill, at a stride of four. **The largest `NB` that compiles at
all is 125** -- the error state `NE = 3 + NB` is what the covariance is sized by,
so this row's bias limit sits three below the other rows' -- **and the largest
`NY` is 128.**

| task stack | largest fitting value, `NB` = `NY`, nothing held | largest fitting `NY`, `NB` held at 4 | largest fitting `NB`, `NY` held at 4 |
|---|---|---|---|
| 4 KiB | no measured value | no measured value | no measured value |
| 8 KiB | 4 | 4 | 4 |
| 16 KiB | 4 | 8 | 4 |
| 32 KiB | 8 | 16 | 8 |
| 48 KiB | 12 | 20 | 12 |
| 64 KiB | 12 | 28 | 16 |

"No measured value" is not "nothing fits": `NB = 3` is a legal configuration --
`NB = 2` is refused at compile time and `NB = 3` is accepted, both checked -- and
the fill's stride does not land on it, so it is not measured. The smallest
measured configuration on the measurement axis needs 4,840 bytes against 4 KiB's
4,096.

### The two dimensions INTERACT, and that is measured rather than assumed

The grid exists to answer whether growth in one dimension depends on the value
held in the other. It does, on all four two-dimensional rows, and the evidence is
in the tables above rather than in a further sweep.

Take the null hypothesis that the cost is additively separable -- that a state
dimension and a measurement dimension each contribute a term and nothing depends
on the pair. Under it, the equal-dimension peak at rung `d` is predictable from
the two held-dimension lines as `W(d, 4) + W(4, d) - W(4, 4)`. Measured against
prediction, at the top rung:

| row | predicted at rung 32 | measured at rung 32 | measured / predicted |
|---|---:|---:|---:|
| `kalman_filter` | 109,496 | 158,712 | 1.45x |
| `ekf` | 146,392 | 245,752 | 1.68x |
| `ukf` | 206,888 | 227,320 | 1.10x |
| `mekf` | 220,072 | 310,840 | 1.41x |

**Every row exceeds its separable prediction on the equal-dimension line**, so a
caller cannot size one dimension from a table taken at another value of the
second, which is exactly why every heading above names what it holds.

### The shape of that interaction, from the interior fill

The four two-dimensional rows were then filled over their whole interior at a
stride of four -- 1,089 points each for `kalman_filter`, `ekf` and `ukf` and
1,023 for `mekf` -- and the same separability ratio evaluated at every one of
them:

| row | interior points | ratio, smallest | ratio, largest | where the largest sits | above the separable prediction |
|---|---:|---:|---:|---|---:|
| `kalman_filter` | 1,088 | 0.700 | 1.543 | `NX` 52, `NY` 52 | 87.4% |
| `ekf` | 1,088 | 0.598 | 2.101 | `NX` 52, `NY` 116 | 88.3% |
| `ukf` | 1,088 | 0.533 | 1.477 | `NX` 1, `NY` 1 | 82.0% |
| `mekf` | 1,022 | 0.806 | 1.906 | `NB` 48, `NY` 128 | 90.9% |

**Four fifths to nine tenths of the interior lies above the separable
prediction**, so the interaction is not a feature of the diagonal. Below the
prediction the ratio falls mostly where one dimension is at the very bottom of
its range, where a fixed part that does not scale with either dimension
dominates. On `ekf` that is the only place it falls below. Away from the first
value of either axis it falls below at one point on `mekf` (12 bias states, 8
outputs), at 12 on `kalman_filter` (8 states with 88 outputs or more, and 12
states with 128), and at 70 on `ukf`, all with the state dimension between 8 and
24 and the measurement dimension anywhere from 8 to 128, none of them below
0.834. The interior points exclude the reference point at 4 and 4, where the
ratio is one by construction, and the smallest ratios exclude the one-state,
one-output points the caveat above places inside the harness gap. `ukf`'s largest
ratio sits at that bottom corner, where the separable prediction is a sum of two
small figures; away from the first value of either axis it is 1.337, at `NX` 24
and `NY` 72.

**AND THE INTERACTION DOES NOT KEEP GROWING WITH THE DIMENSION.** The three-line
tables reach dimension 32, and up to there the excess rises monotonically. The
fill continues to 128 and it does not: on the equal-dimension line
`kalman_filter` runs 1.45 at 32, climbs to 1.48 around 96 and falls back to 1.43
at 128, with one isolated point at 52 reaching 1.54; `ukf` runs 1.10 at 32, rises
to 1.19 at 108 and falls to 1.15 from 112 to 128. Away from `ukf`'s bottom
corner, the largest excess on every row sits at a state dimension of 24 to 52
rather than at the top corner. A reader extrapolating "the excess grows with
dimension" past 32 would be extrapolating past where it holds.

What the fill still does not resolve is anything between two adjacent values of
its stride: it says nothing about the three integers between 96 and 100.

### The largest dimension that COMPILES, per row

This answers the other of the two questions, and it is not the one the tables
above answer. Past the value below, the configuration **does not exist as an
instantiation**: the linear-algebra library refuses to place the object on the
stack and the translation unit does not build, whatever task stack the caller
has.

Each was **walked to, one integer at a time**, upward from a point the interior
fill had already compiled AND run, and each refusal was checked to be that
library's `OBJECT_ALLOCATED_ON_STACK_IS_TOO_BIG` assertion rather than an
unrelated instantiation error at the same dimension. Both sides of every boundary
are recorded: the value below it that builds, and the value at it that does not.

| row | axis | largest that compiles | first refused |
|---|---|---:|---:|
| `kalman_filter` | `NX` with `NY` at 4, `NY` with `NX` at 4, and both equal | 128 | 129 |
| `ekf` | the same three | 128 | 129 |
| `ukf` | the same three | 128 | 129 |
| `manifold_ukf` | `NY`, rotation state at its structural 3 | 128 | 129 |
| `mekf` | `NB`, i.e. error state `NE = 3 + NB` | 125 | 126 |
| `mekf` | `NY` with `NB` at 4 | 128 | 129 |

**The arithmetic these figures agree with, stated so the reader can see why 128
and not some other number.** Each row's dominant fixed-size object is a square
`double` matrix of the axis dimension -- the state covariance for `kalman_filter`
and `ekf`, the same for `ukf` beside a sigma-point set that is a `std::array` of
vectors rather than one object, the innovation covariance on the measurement
axis, and the error-state covariance for `mekf`. The library's fixed-size
allocation limit is 131,072 bytes and its check is `size * sizeof(T) <= limit`,
so `d * d * 8 <= 131,072` gives `d <= 128` exactly, and `mekf`'s selector is
three below its error state, giving `NB <= 125`. **The published numbers are the
walk's, not this arithmetic's**; the arithmetic is here because a number without
one is the kind of constant this document should not carry, and the agreement is
a check rather than a derivation.

The corner of each rectangle is measured rather than assumed: `NX = NY = 128`
builds and runs on all three two-dimensional rows and `NB = 124, NY = 128` on
`mekf`, so the region below the limits is a full rectangle and not only its axes.

**The second limit binds first, by two orders of magnitude.** No entry in any
supported-maximum table above reaches 40, and every instantiation limit is 125 or
128. On these rows a caller who does not fit a dimension has a stack problem, not
an instantiation problem -- which is the opposite of nothing, because it means a
larger task stack is a real answer.

### Growth, stated from the measurement

On the equal-dimension line the whole-chain peak grows by a factor of 3.7 to 4.8
per doubling at rung 32 -- an exponent of 1.87 to 2.26 -- so the quadratic growth
is measured to hold for the whole chain and not only for the dominant object.
`manifold_ukf`, whose only axis is the measurement dimension, grows at 1.41 there.
The exponent is not constant along the axis: every row is markedly sub-quadratic
between the first two rungs, where a fixed overhead that does not scale with the
dimension still dominates, and it settles slightly BELOW quadratic at the top.
Measured over the fill's last doubling, from 64 to 128 on the equal-dimension
line, the exponent is 1.93 (`kalman_filter`), 1.92 (`ekf`), 1.94 (`ukf`), 1.95
(`mekf` from 60 to 120) and 1.64 (`manifold_ukf`).

The deepest points the grid holds, all at the instantiation limit: `ukf`
3,425,464 bytes at 128 states and 128 outputs, `ekf` 3,221,912, `mekf` 3,195,880
at 124 bias states and 128 outputs, `kalman_filter` 2,495,304, and
`manifold_ukf` 565,936 at 128 outputs. **Roughly three megabytes of stack**, on
rows whose `allocation-free?` cell says `YES`.

Against the Riccati row measured under the same instrument, at eight states and
eight outputs the four two-dimensional rows sit at 11,928 to 29,272 bytes where
the discrete Riccati solve at eight states sits at 41,304. The gap closes with
dimension and reverses: at thirty-two those four are past 150,000 bytes and the
Riccati row was not measured there.

### Provenance

`g++ (GCC) 16.2.1 20260810`, `-std=c++20 -O2 -fno-exceptions -fno-rtti -pthread`,
no `-march` (driver default `-mtune=generic -march=x86-64`), x86-64 Linux,
`double`, Eigen 3.4.1, taken on 2026-10-20. Every figure in this section belongs
to that compiler: the three-line tables, the fill, the limits and every ratio and
exponent drawn from them. The same three-line commands run on the library as it
was before the sigma-point factor and the covariance checks reproduced every
cell published under `g++ (GNU) 16.1.1 20260728` exactly, so the compiler change
moves none of them. HOST measurements. Runtime watermarks taken on a pthread
against a **zero-byte harness floor** and a 1,024-byte harness gap. The
three-line tables were taken with a 4 MiB painted region on a 64 MiB thread
stack; **the interior fill raises those to 32 MiB on a 512 MiB stack, and that is
load-bearing rather than cautious.**

**Why the window had to be raised, in one number.** The painted region is a
CEILING on what the instrument can report: a chain that reaches its bottom
disturbs the last word the walk looks at, and the walk then returns the window's
own size rather than the chain's depth -- a figure indistinguishable from a
measurement. The deepest chain in this grid is 3,425,464 bytes, which is 82
percent of the 4 MiB the three-line tables were taken with, and the deepest in
the predictive controller section below is 4,135,224 bytes, which is **98.6
percent of it**. Every record in the fill carries whether it saturated; none did.

**Reproducing the fill and the limits.** The fill is
`tools/stack_watermark.sh --rows estimators --interior`, 4,323 points, journalled
point by point so that a run cut off part-way resumes rather than restarting;
this one was cut off and resumed. `--rows` also takes a comma-separated list of
row names, so one row can be filled alone; without `--rows` the stage fills the
predictive controller row as well. The limits are
`tools/stack_watermark.sh --rows estimators --ceilings`. The fill ran at `-j6`
and `-j3`, and at that parallelism a per-point compile time is a campaign wall
time under the driver's own concurrency rather than a build cost; the driver's
station-quiet gate therefore applies at `-j1` only and the station is censused on
both sides of the stage instead.

**Corpus.** The three vector-state rows are driven over the same forward-Euler
damped chain the discrete Riccati row is measured on -- `-0.5` on the diagonal,
`1.0` on the superdiagonal, step `0.01` -- with output `i` reading state
`i mod NX`, process noise `0.01 I`, measurement noise `I` and initial covariance
`10 I`, at **input dimension 1**. The two attitude rows propagate a constant body
rate on SO(3) and measure the gravity direction repeated to the output dimension,
with process noise scaled by `1e-6`. The two families are not comparable at equal
state dimension, because they are not measurements of one plant.

Reproduce with `tools/stack_watermark.sh --rows estimators --diagonal --frames`
and `--held-dimension`. A target's own frames differ with its ABI, register file
and calling convention -- as the ESP32 capture above shows, where the `float`
peaks run at roughly half the host `double` figures.

**THE LINEAR-ALGEBRA RELEASE IS AN AXIS HERE.** The tables above are Eigen
3.4.1, which is the system release; the project's test tree fetches 3.4.0.
Re-running the whole equal-dimension line against 3.4.0 moves **16 of the 24
measured points**, by -144 to +768 bytes, at most 3.3% of the figure, on
`manifold_ukf`, whose shift is the largest and grows monotonically with the
measurement dimension. That is small, and it is not harmless: it is the size of
the margin by which a supported-maximum entry can clear its task stack, so an
entry close to its stack can read differently under one patch release of a
header-only dependency and the next. `ukf` at four states and four outputs needs
4,280 bytes against 3.4.1 and 4,296 against 3.4.0, past 4 KiB under both. A
frame size is an output of inlining a header-only template library, and a patch
release of that library is therefore part of the provenance of every figure
here. **The held-dimension sweeps were NOT re-measured against 3.4.0**; only the
equal-dimension line was, so nothing is claimed about the other two columns under
that release.

### What is still outstanding here

- **The fill's stride is four, so nothing is measured between two adjacent values
  of it.** The interaction is described rather than located, and the supported
  maxima are exact to within the stride rather than exact. Closing that is an
  exhaustive fill, which is 65,024 points against this one's 4,323.
- **The frame column was not re-measured over the interior.** It belongs to the
  equal-dimension build at the five rungs it is printed at, and taking it costs a
  second compile per point.
- **The linear-algebra release was varied only on the equal-dimension line.** It
  moves 16 of 24 points there, and the fill was taken against 3.4.1 alone.
- **`arm64` and MSVC are absent**, and every figure is `double` at `-O2`. A frame
  size is an optimizer output, so these figures belong to that optimization level
  alone.
- **Construction is not watermarked.** Every row's factory checks each noise and
  initial covariance for symmetry and then with a symmetric eigendecomposition,
  outside the painted window. Its frames are measured but its whole-chain depth
  is not: at dimension 32 the deepest are the decomposition's 17,824 bytes and
  the check's 16,992, and at the `mekf` error state of 35 they are 21,232 and
  20,240 (g++ 16.2.1, Eigen 3.4.1, `-O2`).

## Stack cost of the predictive controller row

The default `nmpc` (`nmpc_static`) is the last row whose `allocation-free?` cell
is a heap statement about a hot path whose frames the caller's dimensions size.
It is measured here rather than in the section above for two reasons: its build
needs the optional nonlinear-programming backend, and **the two dimensions its
frames are sized by are not the two a caller picks.**

**The relation, which is why this row's grid is shaped differently from the
five above.** A caller chooses a state dimension `NX`, an input dimension `NU`
and a horizon `NH`. Those induce

    NV   = (NH + 1) * NX + NH * NU
    MaxM = NX * (NH + 1)

the decision dimension and the constraint bound, and the dominant fixed-size
object is the `NV x NV` Hessian. `NV` and `MaxM` are DERIVED, not selected: most
pairs of them correspond to no reachable configuration, so a grid over them
would describe a surface a caller cannot reach. **This row therefore grids the
three chosen dimensions and reports the induced pair beside every figure. That
departure is forced by the shape of the API and not by what the measurement
costs** -- the five estimator rows grid over both of their caller-controlled
dimensions unchanged, and nothing here is a precedent for taking a cheaper
measurement on a row whose axes a caller can actually select.

### Three numbers per configuration, because three different chains run

Unlike every other row in this document, this one carries three whole-chain
peaks, and at the larger configurations the deepest and the shallowest of them
differ by a factor of more than four:

- **construction** -- `nmpc_static`'s constructor, which builds the formulation
  and sets the backend up. It is offline and it is excluded from both figures
  below, exactly as the `allocation-free?` column excludes it. It is measured
  and published anyway, because a caller that constructs the controller on the
  task's own stack pays it, and because from `NV = 8` upward it is the largest
  of the three. The figure includes the controller object itself, built as a
  local, and with it the solver's inline workspace.
- **first solve** -- the backend's solver instance is emplaced lazily on the
  first solve and only reset on every later one, so the first solve enters a
  setup path that no steady-state solve enters.
- **steady state** -- one warm-started `solve` followed by the caller's own
  plant step, which is the same loop body `nmpc_static_nomalloc_test` arms for
  its allocation proof, after twenty solves have walked the controller into
  steady state.

From `NV = 36` upward the deepest single frame on this row belongs to the
backend's setup on the first solve, `nw_sqp_policy::init`. At `NV` 3, 10 and 21
it belongs to `argmin_solver::emplace_step_solver`, which emplaces the
step-budget driver on the first solve, the path the default `max_time = 0`
takes. Both frames run inside the first-solve chain, so, as on the estimator
rows, the column is a lower bound on the first-solve peak beside it. Each frame
is printed with its owner named so a reader can see which.

| `NX` | `NU` | `NH` | `NV` | `MaxM` | deepest single frame | function owning it | peak, construction | peak, first solve | peak, steady state |
|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| 1 | 1 | 1 | 3 | 2 | 6,736 | `argmin_solver::emplace_step_solver` | 10,848 | 14,048 | 14,048 |
| 2 | 2 | 2 | 10 | 6 | 15,296 | `argmin_solver::emplace_step_solver` | 25,088 | 17,320 | 16,776 |
| 3 | 3 | 3 | 21 | 12 | 44,064 | `argmin_solver::emplace_step_solver` | 67,680 | 46,408 | 26,072 |
| 4 | 4 | 4 | 36 | 20 | 109,904 | `nw_sqp_policy::init` | 171,120 | 117,544 | 49,976 |
| 5 | 5 | 5 | 55 | 30 | 250,544 | `nw_sqp_policy::init` | 372,864 | 259,256 | 94,624 |
| 6 | 6 | 6 | 78 | 42 | 498,496 | `nw_sqp_policy::init` | 726,704 | 507,848 | 172,896 |
| 7 | 7 | 7 | 105 | 56 | 897,152 | `nw_sqp_policy::init` | 1,291,984 | 908,312 | 296,552 |
| 8 | 8 | 8 | 136 | 72 | -- | does not instantiate | does not instantiate | does not instantiate | does not instantiate |

The frame column belongs to the equal-dimension build. The three held-dimension
lines, steady-state peak first, first-solve peak second and construction third:

| rung | `NX` swept, `NU` held at 1, `NH` held at 5 | `NU` swept, `NX` held at 2, `NH` held at 5 | `NH` swept, `NX` held at 2, `NU` held at 1 |
|---:|---|---|---|
| 1 | `NV` 11: 16,184 / 18,920 / 26,352 | `NV` 17: 21,480 / 33,432 / 48,040 | `NV` 5: 14,720 / 14,720 / 14,336 |
| 2 | `NV` 17: 21,480 / 33,432 / 48,040 | `NV` 22: 26,896 / 50,424 / 72,488 | `NV` 8: 15,928 / 15,928 / 19,856 |
| 4 | `NV` 29: 36,248 / 79,864 / 115,776 | `NV` 32: 40,440 / 95,288 / 137,384 | `NV` 14: 19,568 / 25,736 / 36,816 |
| 8 | `NV` 53: 88,968 / 241,656 / 351,056 | `NV` 52: 87,064 / 232,952 / 335,184 | `NV` 26: 32,224 / 66,536 / 95,104 |
| 16 | `NV` 101: 276,680 / 841,656 / 1,221,712 | `NV` 92: 231,688 / 700,696 / 1,003,856 | `NV` 50: 82,432 / 216,200 / 309,056 |
| 20 / 23 / 32 | `NV` 125: 409,592 / 1,279,088 / 1,855,872 | `NV` 127: 420,880 / 1,319,688 / 1,888,816 | `NV` 98: 261,328 / 792,776 / 1,123,904 |
| 21 / 24 / 42 | `NV` 131: does not instantiate | `NV` 132: does not instantiate | `NV` 128: 429,792 / 1,340,280 / 1,895,248 |
| 43 | -- | -- | `NV` 131: does not instantiate |

**The first solve costs up to 3.14 times the steady-state solve**, and the
ratio climbs to between 3.06 and 3.14 at the top of all four lines, while the
two are equal at the smallest configurations, `NV` 3, 5 and 8. **Construction
costs between 4.36 and 4.53 times the steady-state solve, and between 1.41 and
1.45 times the first solve,** at the top of every line. At the two smallest
configurations, `NV` 3 and 5, it is the cheapest of the three.

### Supported maximum by task stack

Strict: no margin for the caller's own frames and none for RTOS overhead. The
ladders are geometric, so an entry names the largest MEASURED rung that fits and
says nothing whatever about the dimensions between rungs.

The entries come from the FIRST-SOLVE peak rather than the steady-state one,
because a caller that solves at all runs the first solve once. Twelve of the
twenty-four cells move between the two readings, and every one of them moves
down: read from the steady-state figure, the 64 KiB row would be 4, 4, 4 and 8.

| task stack | largest fitting equal rung | largest fitting `NX`, `NU` at 1, `NH` at 5 | largest fitting `NU`, `NX` at 2, `NH` at 5 | largest fitting `NH`, `NX` at 2, `NU` at 1 |
|---|---|---|---|---|
| 4 KiB | none exists | none exists | none exists | none exists |
| 8 KiB | none exists | none exists | none exists | none exists |
| 16 KiB | 1 | no measured rung | no measured rung | 2 |
| 32 KiB | 2 | 1 | no measured rung | 4 |
| 48 KiB | 3 | 2 | 1 | 4 |
| 64 KiB | 3 | 2 | 2 | 4 |

**THE FIRST TWO ROWS SAY "NONE EXISTS" RATHER THAN "NO MEASURED RUNG", AND THE
DIFFERENCE IS THE POINT.** The smallest configuration this row has is
`NX = NU = NH = 1`, which induces `NV = 3`, and it needs 14,048 bytes. No
configuration of the predictive controller fits a 4 KiB or an 8 KiB task stack,
and that is a statement about every configuration rather than about the rungs
that were measured.

The shipped allocation proof's own configuration, `NX = 2`, `NU = 1`, `NH = 5`,
so `NV = 17` and `MaxM = 12`, needs 21,480 bytes in steady state and 33,432 on
its first solve. It wants a 48 KiB task and does not fit a 32 KiB one: the first
solve exceeds 32 KiB by 664 bytes, while the steady state alone would fit. This
is the row where the heap statement and the stack statement diverge most
sharply: its steady-state solve allocates exactly nothing and still needs more
than six times the stack a `kalman_filter` needs at four states and four
outputs, and its first solve nearly ten times.

**If the controller is CONSTRUCTED on the same task stack, the ladder above is
not the answer.** The task then has to hold the larger of the construction and
first-solve peaks, which is the construction peak from `NV = 8` upward. Read
that way: at 16 KiB the smallest equal rung (14,048 bytes, its first solve) and
a one-step horizon (14,720, its first solve) fit; at 32 KiB the equal line
reaches rung 2 (25,088), the state line rung 1 (26,352) and the horizon line
rung 2 (19,856), and the input line has nothing; at 48 KiB the state line
reaches rung 2 (48,040), the input line rung 1 (48,040) and the horizon line
rung 4 (36,816), while the equal line stays at rung 2; at 64 KiB none of the
four lines moves further, because the equal line's rung 3 needs 67,680.
Construct off the real-time task, or size the task for construction rather than
for the solve.

### The largest configuration that COMPILES, which this row reaches

This is the other of the two limits, and it answers a different question from the
table above: past it the configuration **does not exist as an instantiation** and
no task stack changes that.

The dominant fixed-size object is the `NV x NV` Hessian, so the limit is a
statement about the decision dimension and the linear-algebra library's
fixed-size allocation limit fixes it exactly: `NV * NV * 8 <= 131,072` gives
`NV <= 128`. **Measured, and approached from three different directions:
`NV = 128` instantiates and `NV = 131` does not.** The refusal is that library's
`OBJECT_ALLOCATED_ON_STACK_IS_TOO_BIG` assertion, and it arrives at
`NX = 21, NU = 1, NH = 5`, at `NX = 2, NU = 24, NH = 5` and at
`NX = 2, NU = 1, NH = 43` -- three configurations with nothing in common except
the decision dimension they induce. A separate probe at `NV = 129` is refused
too, which closes the bracket on the derived bound. Each of the three boundaries
was re-walked, one integer at a time from the value below it, with the refusal
classified: `NX` 20 builds and 21 does not, `NU` 23 builds and 24 does not, `NH`
42 builds and 43 does not.

So a caller can tell the two failures apart on this row: above `NV = 128` the
configuration **does not exist as an instantiation**, and well below it the
configuration exists and **does not fit a small task stack**. The second limit
binds first by a wide margin.

### The supported maximum over the reachable configurations, from the fill

The three ladders above hold two dimensions fixed while sweeping the third. The
interior fill instead enumerates **354 reachable triples** -- every combination of
the three chosen dimensions at a stride of four whose induced decision dimension
is at most 128 -- and reports the supported maximum against the FIRST-SOLVE peak,
for the same reason the ladders do.

| task stack | largest fitting `NV` | every measured configuration below this `NV` fits |
|---|---|---|
| 4 KiB | none | `NV` 3 already does not fit |
| 8 KiB | none | `NV` 3 already does not fit |
| 16 KiB | 9 | 10 |
| 32 KiB | 16 | 17 |
| 48 KiB | 21 | 22 |
| 64 KiB | 25 | 26 |

**THE TWO COLUMNS AGREE AT EVERY TASK STACK IN THE TABLE, AND THAT IS A PROPERTY
OF THESE STACKS RATHER THAN OF THE ROW.** The peak is dominated by the decision
dimension and is not a function of it, so "the largest `NV` that fits" and "the
`NV` below which everything fits" are different questions. Over the 44 decision
dimensions the fill reaches by more than one route, the deepest and shallowest
configurations at the same `NV` differ by up to **1.03x**, at `NV = 9`:
`NX = 1, NU = 1, NH = 4` needs 15,944 bytes and `NX = 4, NU = 1, NH = 1` needs
15,448. Below `NV = 97` no measured configuration needs more stack than one at
a larger decision dimension, so the two columns meet. From `NV = 97` upward ten
pairs of decision dimensions have a configuration at the smaller `NV` that needs
more stack than one at the larger: `NX = 48, NU = 1, NH = 1` induces
`NV = 97` and needs 794,656 bytes, while `NX = 1, NU = 96, NH = 1` induces
`NV = 98` and needs 794,360. The pattern is consistent: of the 335 pairs of
routes to one decision dimension in which one route has the longer horizon and
the smaller state, the longer horizon is the cheaper in 334. The exception is
the `NV = 9` pair above.

Size a task from the configuration, not from the decision dimension it induces.

### The three chosen dimensions INTERACT

The peak is dominated by the decision dimension and is not a function of it.

Asymptotically the whole-chain steady-state peak is about `27 bytes * NV * NV`:
over the seven measured points with `NV >= 92` the coefficient lies between 26.1
and 27.4 while `NX` ranges over 2 to 20 and `MaxM` over 12 to 120. At that end
the decision dimension predicts the peak to within five percent whichever
chosen dimension produced it. The first-solve peak, which the supported maxima
are read from, is about `82 bytes * NV * NV` over the same seven points, between
81.8 and 82.8.

At moderate dimensions the decision dimension predicts the steady-state peak
nearly as well. `NX = 8, NU = 1, NH = 5` induces `NV = 53` and needs 88,968
bytes; `NX = 2, NU = 8, NH = 5` induces `NV = 52` and needs 87,064. A 1.9
percent difference in the decision dimension comes with a 2.2 percent difference
in the peak, and the coefficient there lies between 31.3 and 33.0: 32.2 to 33.0
for the configurations with `NX = 2` and 31.3 to 31.7 for those with `NX >= 5`.
The route matters most at the small end: over the fill, the steady-state peaks
at one decision dimension differ by 4.9 percent at `NV = 9` and by at most 2.1
percent at every larger one.

The additive-separability test the estimator rows use refuses, and not with one
sign. Predicting the equal-dimension peak at rung `d` from the three held lines
as `W(d,1,5) + W(2,d,5) + W(2,1,d) - 2*W(2,1,5)` gives 9,424 against 14,048
measured at `d = 1`, 21,344 against 16,776 at `d = 2` and 53,296 against 49,976
at `d = 4`, so the measured peak is 1.49, 0.79 and 0.94 times the prediction. A
cost that grows quadratically in a sum of the three dimensions is not a sum of
three one-dimensional costs, and the three-axis null over-subtracts the fixed
part at the small end.

The interior fill of 354 reachable triples describes that dependence rather than
only locating it, and the numbers above under the supported maximum are its
statement: at one decision dimension the first-solve peak varies by up to 1.03x
with the route taken to it. What the fill still does not resolve is anything
between two adjacent values of its stride of four.

### Provenance

`g++ (GCC) 16.2.1 20260810`, `-std=c++20 -O2 -fno-exceptions -fno-rtti -pthread`,
no `-march` (driver default `-mtune=generic -march=x86-64`), x86-64 Linux,
`double`, Eigen 3.4.1, taken on 2026-10-20. Every figure in this section belongs
to that compiler. HOST measurements. Runtime watermarks taken on a pthread
against a **zero-byte harness floor** and a 1,024-byte harness gap. The three
ladders were taken with a 4 MiB painted region on a 64 MiB thread stack; the
interior fill raises those to 32 MiB on a 512 MiB stack. This row's deepest
construction chain in the fill is 2,508,920 bytes, at `NX = 4, NU = 120,
NH = 1`: 59.8 percent of the 4 MiB window and 7.5 percent of the 32 MiB one. No
point on this row saturated either window.

**THE BACKEND REVISION IS PART OF EVERY FIGURE IN THIS SECTION.** All of them
were measured against argmin at commit
`0a638cb87b6e959f7cb7b49f0afa34f137b93b69`, fetched from an empty tree at the
revision this project pins and read back out of the fetched source rather than
out of the pin that selected it. The frames above are an output of inlining that
backend's templates, so **figures taken against two revisions of it are not
comparable**, and a figure from this section quoted without the commit beside it
is one a reader cannot re-derive. Re-run the instrument under your own pairing
rather than reading across.

**Corpus.** The same forward-Euler damped chain the Riccati and estimator rows
are measured over -- `-0.5` on the diagonal, `1.0` on the superdiagonal, step
`0.01` -- used as the controller's prediction model, with a state weight of
`10 I`, an input weight of `0.1 I`, a zero reference, an initial state whose
first component is one, and no path or terminal constraint. Every measured
configuration solved.

Reproduce the ladders with `tools/stack_watermark.sh --rows controller
--diagonal` and `--held-dimension`, the fill with `--interior --rows controller`
and the two limits with `--ceilings --rows controller`; without `--rows`, those
two stages cover the estimator rows as well. The driver configures a tree from
empty to fetch the backend when it does not already have one.

### What is still outstanding here

- **The fill's stride is four**, so it reaches 354 of the 7,611 reachable
  configurations and says nothing between two adjacent values of a stride.
- **The frame column belongs to the equal-dimension ladder** and was not
  re-measured over the fill.
- **`arm64` and MSVC are absent**, every figure is `double` at `-O2`, and one
  backend revision was measured. A frame size is an optimizer output.
- **The linear-algebra release was not varied for this row.** It moves published
  cells on the estimator rows by up to 3.4 percent, and nothing here bounds what
  it does to these.

## Wall-clock budgets are not RT-safe

The optimization solver settings expose two kinds of budget.  The iteration
caps (`nlopt_settings::max_eval`, `argmin_settings::max_eval`, and the OSQP
`max_iter` limit) bound the solve by a deterministic count of steps.  The
wall-clock budget (`nlopt_settings::max_time`, `argmin_timed_settings::max_time`)
bounds the solve by elapsed time, which depends on machine load, frequency
scaling, and cache state, so two identical calls can stop at different
iterates.  The wall-clock budget is therefore **not real-time safe** and must
not be selected on a real-time path.  The two backends select it differently.
NLopt's is a run-time field that defaults to `0` (disabled), so the RT-safe
selection is to leave it disabled.  argmin's is a compile-time choice: the
default `argmin_step_budget` has no `max_time` field and instantiates no clock
read, and only an `argmin_solver` built with `argmin_step_and_time_budget`
carries the field, which also defaults to `0`.

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
once warmed up: measured at 0.00 allocs/step on a converging trajectory with
`nw_sqp`, and at an exact regulation equilibrium, where every solve takes the
SQP null step, with both `nw_sqp` and `filter_nw_sqp` on the static path.  It
remains iteration-bounded by `max_eval` and clock-free under its default step
budget, so on those axes it meets the hard real-time obligations the
runtime-horizon `nmpc_dynamic` path cannot.

**Measured on silicon, the default `nmpc` does not fit its demonstration's
period at the worst case.** No runtime-horizon optimizer runs on either board.
The one predictive figure on silicon is the default `nmpc` at the shipped pin
(`NV = 17`, `MaxM = 12`), on the NUCLEO-H753ZI's timing image (build id
`5e0486618196df29e7e73fc391fd1103515f6d6a`,
`scripts/board_capture_fixtures/timing.txt`, `arm-none-eabi-gcc` 16.2.0, 64 MHz
with both caches off; see "On silicon: the NUCLEO-H753ZI leg" above). There, one
solve and its plant step took a median of 76.3 ms and a worst of 151.2 ms over
the 201 window solves, and the cold first solve after reset took 2.22 s. The
demonstration's control period is 100 ms (`kPredictiveDt` in
`examples/embedded/shared/golden_reference.h`). **The worst window solve
overruns that period, and so does the cold first solve; only the median fits
within it.** A bounded iteration count and zero steady-state allocation do not
bound the solve's time inside a given period, so a caller on a deadline sizes it
against the measured worst case at its own posture. The ESP32-POE-EA has no
predictive figure: both of its images declare the predictive family absent
(`scripts/board_capture_fixtures/esp32-evidence.txt`, ELF sha256
`c7bc0e9faef5707f2ec9d327766ecc397c84312c91872c51383a3e628aa25f98`, and
`scripts/board_capture_fixtures/esp32-timing.txt`, ELF sha256
`a7d830451598f39f56a1b579900c59613d420dc5f89e872a18d0e2fcbdb7fbd5`), because
the backend's solver computes in `double` only and its first solve's stack
exceeds that leg's 8,192-byte control task.

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
