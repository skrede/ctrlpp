# Determinism and RT-Safety Matrix

Per-module real-time safety and determinism guarantees for the ctrlpp hot
paths.  Every YES cell cites the test or build artifact that proves it; no
cell in this table is self-reported.  The evidence was re-run in full before
this table was written.

## Reading the columns

| Column | Meaning |
|--------|---------|
| allocation-free? | The steady-state hot path (compute/predict/update/evaluate/sample) performs zero heap allocation. Construction and setup may allocate. **THIS COLUMN IS A HEAP STATEMENT AND SAYS NOTHING ABOUT THE STACK.** For a row whose hot path holds a handful of fixed-size matrices that is the whole resource story. For a row that runs a fixed-size decomposition whose dimension the CALLER chooses, it is not: those frames grow with that dimension, and on a small task stack the stack is what breaks first. Every such row has a stack section below the table -- see "Stack cost of the Riccati solve" and "Stack cost of the estimator rows" -- and a `YES` here must not be read as a resource claim on its own. |
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
| Riccati steady-state solve: `dare` / `care` | YES **on the heap, and the heap is not the binding cost here** (see the stack note below this table, which a hard-real-time caller must read before sizing a task stack: at `NX = 8` the whole chain reaches 41,304 bytes) | YES (Eigen Schur iteration bound; sign-function Newton capped at `max_iters = 40` in `detail/care_sign_function.h`, followed by a fixed-size rescale-and-factorize extraction. On every accepted solve, on **every** continuous method tag and on the discrete solver, fixed-size checks of the returned solution run before it is reported: the discrete path solves a fixed-size Stein system of dimension `NX(NX+1)/2` on the symmetric subspace to estimate the answer's own forward error, and every continuous path evaluates the counted Riccati residual and one non-accumulating real Schur factorization of the closed-loop spectrum, from the single definition in `detail/care_postconditions.h`. All three continuous tags reach that definition -- the sign-function path, the real-Schur path and the balanced-Schur path, the last verifying against the caller's Hamiltonian rather than the balanced one -- so no tag carries a bound the others do not. The balanced tag additionally runs a DGEBAL-style balance ahead of all of this. It allocates nothing, and it terminates because every accepted rescale strictly reduces that row's norm sum by the factor 19/20 using power-of-two scalings clamped to the scalar type's range. Its sweep count is data-dependent, and it is now **bounded at compile time** the way `max_iters` bounds the Newton loop: the cap is `2 * NX * (max_exponent - min_exponent)`, derived from the scalar type and the dimension rather than tuned, and reaching it stops the sweep and returns the partially balanced matrix. Stopping early is safe because balancing is a similarity preconditioner and every applied step updates `H` and `D` together, so `H_returned == D^-1 * H_original * D` holds exactly at any cut point -- a capped run is less well balanced, never wrong. The cap is a guarantee rather than a budget: over 550,000 draws with entries spread across `2^-280` to `2^+280`, the worst sweep counts observed were 2 / 12 / 43 / 45 at NX = 1 / 2 / 4 / 8 against ceilings of 4,090 / 8,180 / 16,360 / 32,720, so the worst measured run sits about three orders of magnitude below its ceiling. Budget from the measurement and treat the cap as the backstop) | YES | YES | YES | `dare_care_nomalloc_test` (NX = 2, 4, 8 across the Schur, sign-function, and balanced-Schur variants); `care_convergence_anchor_test` (fixed-seed near-axis and simple-input scale sweeps under all three method tags, plus magnitude-band anchors); `riccati_magnitude_test` (both ends of the arithmetic range); leg 1 witness calls `dare` and `care`. The allocation cell's heap claim is proved by `dare_care_nomalloc_test`; its stack claim is proved by the `-fstack-usage` frames and runtime watermarks in the stack note below, which are host `double` measurements, now corroborated on silicon by the ESP32 stack probe in `examples/embedded/esp32/main/app_main.cpp` -- the on-target `float` peaks track the host `double` prediction to within the scalar width, and that probe also establishes that at `float` the ACCEPTANCE CHECK, not the stack, is what bounds the usable state dimension |
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

Measured on an ESP32-WROOM (Xtensa LX6, 240 MHz, single-precision FPU), ESP-IDF
v6.0.2, `xtensa-esp32-elf-g++` 15.2.0 at `-Os` with exceptions and RTTI off,
Eigen 3.4.0, one FreeRTOS task per dimension so each figure is that dimension's
own peak rather than a running minimum over a sweep:

| `NX` | on-silicon peak, `float` | host peak, `double` | host / board | accepted? |
|---:|---:|---:|---:|:--|
| 2 | 3,100 | 4,552 | 1.47x | solved |
| 3 | 4,160 | -- | -- | refused |
| 4 | 5,932 | 11,176 | 1.88x | solved |
| 5 | 8,336 | -- | -- | refused |
| 6 | 11,500 | 21,912 | 1.91x | refused |
| 8 | 20,560 | 41,304 | 2.01x | refused |

**The ratio converges on exactly the scalar width.** It is 1.47x at `NX = 2` and
climbs to 2.01x at `NX = 8`, which is what a halved scalar predicts for a chain
whose cost is dominated by `NX`-dimensioned arrays and diluted at small `NX` by
fixed overhead that does not scale with the scalar. The host table is therefore
usable as written for `double` and halves for `float`; it is not refuted and it
is not to be replaced by these figures. **The two columns do not share an input
dimension**: the host column is a single-input pose throughout, while the board
probe drives one input at `NX = 2, 3, 5` and two at `NX = 4, 6, 8`. The ratio is
therefore a statement about the scalar across that pair of pose families and not
a conversion factor to apply cell by cell.

**Supported maximum ON THE STACK for this board: `NX = 4`, leaving 2,256 bytes of
an 8,192-byte task stack.** That is the stack answer only -- the accuracy gate
refuses `NX = 3` and everything above `NX = 4` at `float`, and clears `NX = 4`
itself by just 3.4 percent, so the usable maximum is lower than this line alone
implies. See the accuracy paragraph below. Established twice and by different
means. A pass giving every
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

Witness: `examples/embedded/esp32/main/app_main.cpp`, whose control loop still
designs its gain, runs its 201 steps, streams them over UART2 and reports
`golden diff PASS` in the same boot that carries the probe.

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

### What each column is, and what the ladder can and cannot resolve

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

The second table is the supported maximum by task stack. The dimension ladder is
GEOMETRIC, so an entry names the largest MEASURED rung that fits and says nothing
whatever about the dimensions between rungs. It is not a claim that the next
integer up does not fit.

The harness floor read zero on every one of the sixty-eight measurements behind
the tables below, printed beside each figure by the instrument, so every peak is
attributable to the call rather than to the harness. Each figure reproduced
byte-identically across two independent runs of the whole grid.

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

| task stack | largest fitting rung, `NX` = `NY`, nothing held | largest fitting `NY`, `NX` held at 4 | largest fitting `NX`, `NY` held at 4 |
|---|---|---|---|
| 4 KiB | 4 | 4 | 4 |
| 8 KiB | 4 | 8 | 8 |
| 16 KiB | 8 | 16 | 8 |
| 32 KiB | 8 | 16 | 16 |
| 48 KiB | 16 | 32 | 16 |
| 64 KiB | 16 | 32 | 16 |

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

| task stack | largest fitting rung, `NX` = `NY`, nothing held | largest fitting `NY`, `NX` held at 4 | largest fitting `NX`, `NY` held at 4 |
|---|---|---|---|
| 4 KiB | 4 | 4 | 4 |
| 8 KiB | 4 | 8 | 4 |
| 16 KiB | 8 | 8 | 8 |
| 32 KiB | 8 | 16 | 16 |
| 48 KiB | 8 | 16 | 16 |
| 64 KiB | 16 | 32 | 16 |

### `ukf`

The frame column belongs to the equal-dimension build. The sigma-point set is
`2 * NX + 1` points held live across the propagation, so this row pays for the
state dimension twice: once in the covariance it factors and once in the set.

| dimension | deepest single frame | function owning it | peak, `NX` = `NY`, nothing held | peak, `NY` swept, `NX` held at 4 | peak, `NX` swept, `NY` held at 4 |
|---:|---:|---|---:|---:|---:|
| 2 | 1,504 | `Eigen::internal::make_block_householder_triangular_factor` | 2,344 | 2,952 | 3,112 |
| 4 | 2,704 | `ukf::update` | 4,088 | 4,088 | 4,088 |
| 8 | 8,704 | `ukf::update` | 12,272 | 7,544 | 8,640 |
| 16 | 35,376 | `ukf::update` | 46,976 | 20,344 | 26,832 |
| 32 | 135,568 | `ukf::update` | **212,488** | 60,008 | 128,648 |

Strict: no margin for the caller's own frames and none for RTOS overhead.

| task stack | largest fitting rung, `NX` = `NY`, nothing held | largest fitting `NY`, `NX` held at 4 | largest fitting `NX`, `NY` held at 4 |
|---|---|---|---|
| 4 KiB | 4 | 4 | 4 |
| 8 KiB | 4 | 8 | 4 |
| 16 KiB | 8 | 8 | 8 |
| 32 KiB | 8 | 16 | 16 |
| 48 KiB | 16 | 16 | 16 |
| 64 KiB | 16 | 32 | 16 |

The 4 KiB entry in the first column clears by **8 bytes** at four states and four
outputs: 4,088 against 4,096. Read it as "does not fit" unless the whole rest of
the task is the empty frame this measurement was taken in -- and see the
linear-algebra release paragraph below, which turns that eight bytes into a
measured sign flip rather than a caution.

### `manifold_ukf`

**This row carries ONE caller-controlled dimension.** Its state is a rotation, so
the covariance it propagates is 3 x 3 whatever the caller does, and there is no
second axis to sweep or hold. Every figure below is taken with the rotation state
at its structural three and the measurement dimension swept.

| `NY` | deepest single frame | function owning it | whole-chain peak, rotation state held at 3 |
|---:|---:|---|---:|
| 2 | 1,504 | `Eigen::internal::make_block_householder_triangular_factor` | 3,008 |
| 4 | 2,048 | `manifold_ukf::update` | 3,776 |
| 8 | 3,312 | `manifold_ukf::update` | 5,648 |
| 16 | 12,288 | `manifold_ukf::update` | 17,936 |
| 32 | 31,040 | `manifold_ukf::update` | **49,472** |

Strict: no margin for the caller's own frames and none for RTOS overhead.

| task stack | largest fitting `NY`, rotation state held at 3 |
|---|---|
| 4 KiB | 4 |
| 8 KiB | 8 |
| 16 KiB | 8 |
| 32 KiB | 16 |
| 48 KiB | 16 |
| 64 KiB | 32 |

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
| 16 | 19 | 68,816 | `mekf::update` | 84,904 | 25,112 | 57,592 |
| 32 | 35 | 262,784 | `mekf::update` | **310,840** | 74,200 | 152,648 |

Strict: no margin for the caller's own frames and none for RTOS overhead.

| task stack | largest fitting rung, `NB` = `NY`, nothing held | largest fitting `NY`, `NB` held at 4 | largest fitting `NB`, `NY` held at 4 |
|---|---|---|---|
| 4 KiB | no measured rung | no measured rung | no measured rung |
| 8 KiB | 4 | 4 | 4 |
| 16 KiB | 4 | 8 | 4 |
| 32 KiB | 8 | 16 | 8 |
| 48 KiB | 8 | 16 | 8 |
| 64 KiB | 8 | 16 | 16 |

"No measured rung" is not "nothing fits": `NB = 3` is a legal configuration and
was not measured, and the smallest measured configuration on the measurement axis
needs 4,840 bytes against 4 KiB's 4,096.

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
| `ukf` | 184,568 | 212,488 | 1.15x |
| `mekf` | 220,072 | 310,840 | 1.41x |

**Every row exceeds its separable prediction, and the excess grows with the
dimension**: `ekf` runs 1.08x above prediction at rung 2 and 1.68x above it at
rung 32. So a caller cannot size one dimension from a table taken at another
value of the second, which is exactly why every heading above names what it
holds.

What this does NOT establish is the SHAPE of that interaction away from the three
lines measured. Each axis was swept at one held value, so the tables locate the
interaction and do not describe it. Describing it is what an interior fill would
do, and no interior point was measured here.

### Growth, stated from the measurement

On the equal-dimension line the whole-chain peak grows by a factor of 3.7 to 4.8
per doubling at the top rung -- an exponent of 1.87 to 2.26 -- so the quadratic
growth is measured to hold for the whole chain and not only for the dominant
object. `manifold_ukf`, whose only
axis is the measurement dimension, grows at 1.46. The exponent is not constant
down the ladder: every row is markedly sub-quadratic between the first two rungs,
where a fixed overhead that does not scale with the dimension still dominates.

Against the Riccati row measured under the same instrument, at eight states and
eight outputs the four two-dimensional rows sit at 11,928 to 29,272 bytes where
the discrete Riccati solve at eight states sits at 41,304. The gap closes with
dimension and reverses: at thirty-two those four are past 150,000 bytes and the
Riccati row was not measured there.

### Provenance

`g++ (GNU) 16.1.1 20260728`, `-std=c++20 -O2 -fno-exceptions -fno-rtti -pthread`,
no `-march` (driver default `-mtune=generic -march=x86-64`), x86-64 Linux,
`double`, Eigen 3.4.1. HOST measurements. Runtime watermarks taken on a pthread
with a 64 MiB stack against a **zero-byte harness floor** and a 1,024-byte
harness gap; the painted region is 4 MiB.

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

**THE LINEAR-ALGEBRA RELEASE IS AN AXIS HERE, AND IT FLIPS A PUBLISHED CELL.**
The tables above are Eigen 3.4.1, which is the system release; the project's test
tree fetches 3.4.0. Re-running the whole equal-dimension line against 3.4.0 moves
**15 of the 24 measured points**, by -144 to +768 bytes -- at most 3.4% of the
figure, on `manifold_ukf`, whose shift is the largest and grows monotonically
with the measurement dimension. That is small, and it is not harmless: `ukf` at
four states and four outputs needs 4,088 bytes against 3.4.1 and **4,104 against
3.4.0**, so its 4 KiB entry above reads 4 under one patch release of a
header-only dependency and 2 under the next. A frame size is an output of
inlining a header-only template library, and a patch release of that library is
therefore part of the provenance of every figure here. **The held-dimension
sweeps were NOT re-measured against 3.4.0**; only the equal-dimension line was,
so nothing is claimed about the other two columns under that release.

### What is still outstanding here

- **The interior of the grid is not measured.** Each axis was swept at one held
  value. The interaction is established; its shape is not.
- **`nmpc` carries no stack figure at all.** Its decision dimension `NV` and
  constraint bound `MaxM` are derived from the horizon, the state dimension and
  the input dimension rather than chosen directly, and its build needs an
  optional backend, so it is measured separately and is not in this section.
- **The largest dimension that COMPILES is not published for any row here.** The
  ladder stops at thirty-two because that is where it stops, not because
  thirty-three fails. A caller must be able to tell "does not fit your stack"
  from "does not exist as an instantiation", and only the first of those is
  answered above.
- **`arm64` and MSVC are absent**, and every figure is `double` at `-O2`. A frame
  size is an optimizer output, so these figures belong to that optimization level
  alone.

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
