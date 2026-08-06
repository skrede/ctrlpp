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
  `max_time` budget.  Every other core header is clock-free.

## The matrix

| module | allocation-free? | bounded-iterations? | wall-clock-free? | exceptions-off-clean? | deterministic(seeded)? | evidence |
|--------|------------------|---------------------|------------------|-----------------------|------------------------|----------|
| `pid` | YES | YES (no loop) | YES | YES | YES | `pid_nomalloc_test` (position form, velocity form, composed anti-windup and derivative filter); leg 1 + `embedded_core_float` |
| `lqr` / `lqr_time_varying` steady-state control law | YES | YES (no loop in the gain application) | YES | YES | YES | lqr TEST_CASEs in `pid_nomalloc_test.cpp` ("lqr steady-state control law", "lqr_time_varying steady-state control law"); leg 1 witness instantiates `lqr_gain` |
| `kalman_filter` | YES below a measurement dimension of 48, **NO at and above it** at two states or more, where the gain solve's Householder sequence goes blocked and the linear-algebra library heap-allocates its block reflector factor once per update (see "The heap claim's boundary" below). A single-state filter stays allocation-free at every measurement dimension. | YES (closed form) | YES | YES | YES | `estimation_kalman_nomalloc_test` (26 configurations spanning 1 to 128 on both axes: the 13 its measurement table prints, the 6 values only its supported-maximum table names, the state axis at its instantiation limit, both sides of the boundary at 47 and 48, and the single-state exemption); `estimation_nomalloc_test` (the fixed instantiation at two states and one output); leg 1 + `embedded_core_float` |
| `ekf` | YES below a measurement dimension of 48, **NO at and above it** at two states or more, for the same reason and with the same single-state exemption as the row above. | YES (closed form) | YES | YES | YES | `estimation_ekf_nomalloc_test` (25 configurations spanning 1 to 128 on both axes: the 13 its measurement table prints, the 6 values only its supported-maximum table names, the state axis at its instantiation limit, both sides of the boundary at 47 and 48, and the single-state exemption); `estimation_nomalloc_test` (the fixed instantiation at two states and one output); leg 1 + `embedded_core_float` |
| `ukf` | YES at every measured configuration, including 128 states by 128 outputs, **because its gain decomposition defaults to LDLT and builds no Householder sequence at all**. Selecting the QR option puts it on the boundary the rows beside it have, at a measurement dimension of 48. | YES (fixed sigma-point set) | YES | YES | YES | `estimation_ukf_nomalloc_test` (25 configurations spanning 2 to 128 on both axes, including both instantiation limits and the corner, plus 2 more with the QR option selected, walking both sides of the boundary that option brings back); `estimation_nomalloc_test` (the fixed instantiation at two states and one output); leg 1 + `embedded_core_float` |
| `mekf` | YES below a measurement dimension of 48, **NO at and above it** (see "The heap claim's boundary" below). The boundary is on the measurement dimension and not on the error state: 125 bias states against 4 outputs is allocation-free. This row has no single-state exemption, since its right-hand side is the error state and is never one column. | YES (closed form) | YES | YES | YES | `estimation_mekf_nomalloc_test` (20 configurations: the 11 of the 13 its measurement table prints that instantiate, the 4 values only its supported-maximum table names, both axes at their instantiation limits with the corner, and both sides of the boundary at 47 and 48); `estimation_nomalloc_test` (the fixed instantiation at three bias states and three outputs); leg 1 + `embedded_core_float` |
| `manifold_ukf` | YES below a measurement dimension of 48, **NO at and above it** (see "The heap claim's boundary" below). No single-state exemption: the right-hand side is the rotation's three columns. | YES (fixed sigma-point set) | YES | YES | YES | `estimation_manifold_ukf_nomalloc_test` (12 configurations spanning 2 to 128 on its one caller axis: the 5 its measurement table prints, the 4 values only its supported-maximum table names, the axis at its instantiation limit, and both sides of the boundary at 47 and 48); `estimation_nomalloc_test` (the fixed instantiation at three outputs); leg 1 + `embedded_core_float` |
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
| `nmpc` (the DEFAULT; = `nmpc_static`, compile-time horizon, argmin `nw_sqp`, bounded decision `NV` + constraint `MaxM`) | YES below a decision dimension `NV` of 48, strict-zero at 0.00 allocs/step in steady state, and **NO at and above it**, where the solve crosses the same blocked-Householder boundary the estimator rows do (see "The heap claim's boundary" below). **The heap is not the binding cost here either** (see "Stack cost of the predictive controller row" below, which a hard-real-time caller must read before sizing a task stack: the configuration this cell's own evidence pins needs 33,304 bytes of stack in steady state, and constructing it needs 87,768) | YES when capped (`max_eval`) | YES with `max_time = 0` (the default); the `max_time` budget is non-RT | not covered by this witness (argmin's `-fno-exceptions` instantiation is clean upstream; the ctrlpp-side no-exceptions dogfood is a separate witness) | YES — the constraint bound feeds only the QP result-multiplier storage, never the compute workspace, so argmin's `nw_sqp` bit-identity golden is unchanged | `nmpc_static_nomalloc_test`: the shipped pin (double_integrator NX=2 NU=1 NH=5 → NV=17, MaxM=12; sentinel `eigen_assert` + `EIGEN_RUNTIME_NO_MALLOC` + global `operator new` counter; `static_assert(strict_allocation_free)`), plus a two-configuration bracket over the damped chain that walks the boundary, at `NV` 47 and 48. **That is 2 of the 25 configurations this row publishes a figure for**, and the reason is measured: a unit holding all 25 peaked at 21.9 GB of compiler resident set and was killed by the kernel, and one holding 8 built and ran green at all 8 but cost 17.4 GB to compile and 232 s to run against roughly a second for the pin alone |
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
  for**. The interior fill behind their supported-maximum tables is 4,677
  configurations and is **not** individually armed: nothing above is a claim
  about a dimension lying between two armed ones.
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
the three-line tables and on every one of the 4,677 behind the fill, printed
beside each figure by the instrument, so every peak is attributable to the call
rather than to the harness. The three-line figures reproduced byte-identically
across two independent runs of the whole grid, and **the fill reproduces all 44
of the published cells it covers, byte for byte**, having been taken at a
different painted window and a different build parallelism.

**One caveat at the very bottom of the grid.** The measurement has a resolution
floor: a chain shallower than the harness gap of 1,024 bytes disturbs nothing the
walk can see and reports zero, which is the value the harness floor reports.
Three points in the fill are there -- `kalman_filter`, `ekf` and `ukf` at one
state and one output. Re-measured at gaps of 512, 256, 128 and 64 bytes, with a
zero floor at every one, they resolve to **720, 880 and 1,008 bytes** and are
byte-identical across all four gaps. No other point in the grid is inside the
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
| 2 | 1,504 | `Eigen::internal::make_block_householder_triangular_factor` | 2,344 | 2,952 | 3,112 |
| 4 | 2,704 | `ukf::update` | 4,088 | 4,088 | 4,088 |
| 8 | 8,704 | `ukf::update` | 12,272 | 7,544 | 8,640 |
| 16 | 35,376 | `ukf::update` | 46,976 | 20,344 | 26,832 |
| 32 | 135,568 | `ukf::update` | **212,488** | 60,008 | 128,648 |

Strict: no margin for the caller's own frames and none for RTOS overhead.

From the interior fill, at a stride of four. **The largest dimension that
compiles at all is 128 on both axes** -- the sigma-point set is a `std::array` of
vectors rather than one fixed-size object, so the allocation limit does not bind
on it and this row's limit is the same as the two above.

| task stack | largest fitting value, `NX` = `NY`, nothing held | largest fitting `NY`, `NX` held at 4 | largest fitting `NX`, `NY` held at 4 |
|---|---|---|---|
| 4 KiB | 4 | 4 | 4 |
| 8 KiB | 4 | 8 | 4 |
| 16 KiB | 8 | 12 | 12 |
| 32 KiB | 12 | 20 | 16 |
| 48 KiB | 16 | 28 | 20 |
| 64 KiB | 16 | 32 | 24 |

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

From the interior fill, at a stride of four. **The largest `NY` that compiles at
all is 128.**

| task stack | largest fitting `NY`, rotation state held at 3 |
|---|---|
| 4 KiB | 4 |
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
| 16 | 19 | 68,816 | `mekf::update` | 84,904 | 25,112 | 57,592 |
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
| `ukf` | 184,568 | 212,488 | 1.15x |
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
| `ukf` | 1,088 | 0.533 | 1.393 | `NX` 28, `NY` 72 | 80.8% |
| `mekf` | 1,023 | 0.806 | 1.906 | `NB` 48, `NY` 128 | 90.8% |

**Four fifths to nine tenths of the interior lies above the separable
prediction**, so the interaction is not a feature of the diagonal. Below the
prediction the ratio falls only where one dimension is at the very bottom of its
range, where a fixed part that does not scale with either dimension dominates.

**AND THE INTERACTION DOES NOT KEEP GROWING WITH THE DIMENSION.** The three-line
tables reach dimension 32, and up to there the excess rises monotonically. The
fill continues to 128 and it does not: on the equal-dimension line
`kalman_filter` runs 1.45 at 32, peaks near 1.48 around 96 and falls back to 1.43
at 128; `ukf` runs 1.15 at 32, peaks near 1.22 around 80 to 96 and falls to 1.13
at 128. The largest excess on every row sits in the MIDDLE of the grid rather
than at its corner. A reader extrapolating "the excess grows with dimension" past
32 would be extrapolating past where it holds.

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
`manifold_ukf`, whose only axis is the measurement dimension, grows at 1.46 there.
The exponent is not constant along the axis: every row is markedly sub-quadratic
between the first two rungs, where a fixed overhead that does not scale with the
dimension still dominates, and it settles slightly BELOW quadratic at the top.
Measured over the fill's last doubling, from 64 to 128 on the equal-dimension
line, the exponent is 1.93 (`kalman_filter`), 1.92 (`ekf`), 1.87 (`ukf`), 1.95
(`mekf` from 60 to 120) and 1.65 (`manifold_ukf`).

The deepest points the grid holds, all at the instantiation limit: `ekf`
3,221,736 bytes at 128 states and 128 outputs, `ukf` 3,037,928, `mekf` 3,195,880
at 124 bias states and 128 outputs, `kalman_filter` 2,495,304, and
`manifold_ukf` 562,576 at 128 outputs. **Roughly three megabytes of stack**, on
rows whose `allocation-free?` cell says `YES`.

Against the Riccati row measured under the same instrument, at eight states and
eight outputs the four two-dimensional rows sit at 11,928 to 29,272 bytes where
the discrete Riccati solve at eight states sits at 41,304. The gap closes with
dimension and reverses: at thirty-two those four are past 150,000 bytes and the
Riccati row was not measured there.

### Provenance

`g++ (GNU) 16.1.1 20260728`, `-std=c++20 -O2 -fno-exceptions -fno-rtti -pthread`,
no `-march` (driver default `-mtune=generic -march=x86-64`), x86-64 Linux,
`double`, Eigen 3.4.1. HOST measurements. Runtime watermarks taken on a pthread
against a **zero-byte harness floor** and a 1,024-byte harness gap. The
three-line tables were taken with a 4 MiB painted region on a 64 MiB thread
stack; **the interior fill raises those to 32 MiB on a 512 MiB stack, and that is
load-bearing rather than cautious.**

**Why the window had to be raised, in one number.** The painted region is a
CEILING on what the instrument can report: a chain that reaches its bottom
disturbs the last word the walk looks at, and the walk then returns the window's
own size rather than the chain's depth -- a figure indistinguishable from a
measurement. The deepest chain in this grid is 3,221,736 bytes, which is 77
percent of the 4 MiB the three-line tables were taken with, and the deepest in
the predictive controller section below is 4,135,224 bytes, which is **98.6
percent of it**. Every record in the fill carries whether it saturated; none did.

**Reproducing the fill and the limits.** The fill is
`tools/stack_watermark.sh --interior`, 4,677 points, journalled point by point so
that a run cut off part-way resumes rather than restarting -- this one was cut
off once and did. The limits are `tools/stack_watermark.sh --ceilings`. The fill
ran at `-j6`, and at that parallelism a per-point compile time is a campaign wall
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

- **The fill's stride is four, so nothing is measured between two adjacent values
  of it.** The interaction is described rather than located, and the supported
  maxima are exact to within the stride rather than exact. Closing that is an
  exhaustive fill, which is 72,635 points against this one's 4,677.
- **The frame column was not re-measured over the interior.** It belongs to the
  equal-dimension build at the five rungs it is printed at, and taking it costs a
  second compile per point.
- **The linear-algebra release was varied only on the equal-dimension line.** It
  moves 15 of 24 points there and flips one supported-maximum entry, and the fill
  was taken against 3.4.1 alone.
- **`arm64` and MSVC are absent**, and every figure is `double` at `-O2`. A frame
  size is an optimizer output, so these figures belong to that optimization level
  alone.

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
peaks and they differ by a factor of four:

- **construction** -- `nmpc_static`'s constructor, which builds the formulation
  and sets the backend up. It is offline and it is excluded from both figures
  below, exactly as the `allocation-free?` column excludes it. It is measured
  and published anyway, because a caller that constructs the controller on the
  task's own stack pays it, and because it is by a wide margin the largest of
  the three.
- **first solve** -- the backend's solver instance is emplaced lazily on the
  first solve and only reset on every later one, so the first solve enters a
  setup path that no steady-state solve enters.
- **steady state** -- one warm-started `solve` followed by the caller's own
  plant step, which is the same loop body `nmpc_static_nomalloc_test` arms for
  its allocation proof, after twenty solves have walked the controller into
  steady state.

**The deepest single frame on this row belongs to the CONSTRUCTOR**, and it
exceeds the whole-chain peak of both solves. So unlike the estimator rows, where
that column is a lower bound on the chain beside it, here it is a statement
about a call neither solve figure covers. It is printed with its owner named so
a reader can see which.

| `NX` | `NU` | `NH` | `NV` | `MaxM` | deepest single frame | function owning it | peak, construction | peak, first solve | peak, steady state |
|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| 1 | 1 | 1 | 3 | 2 | 8,608 | `argmin_solver::emplace_timed_solver` | 19,216 | 13,936 | 13,936 |
| 2 | 2 | 2 | 10 | 6 | 20,432 | `nmpc_static::nmpc_static` | 44,576 | 16,696 | 16,232 |
| 3 | 3 | 3 | 21 | 12 | 59,456 | `nmpc_static::nmpc_static` | 123,376 | 43,688 | 43,688 |
| 4 | 4 | 4 | 36 | 20 | 153,984 | `nmpc_static::nmpc_static` | 314,352 | 107,144 | 91,256 |
| 5 | 5 | 5 | 55 | 30 | 340,192 | `nmpc_static::nmpc_static` | 688,272 | 234,632 | 213,864 |
| 6 | 6 | 6 | 78 | 42 | 666,736 | `nmpc_static::nmpc_static` | 1,344,048 | 459,144 | 319,592 |
| 7 | 7 | 7 | 105 | 56 | 1,189,808 | `nmpc_static::nmpc_static` | 2,392,656 | 819,288 | 557,880 |
| 8 | 8 | 8 | 136 | 72 | -- | does not instantiate | does not instantiate | does not instantiate | does not instantiate |

The frame column belongs to the equal-dimension build. The three held-dimension
lines, steady-state peak first and first-solve peak second:

| rung | `NX` swept, `NU` held at 1, `NH` held at 5 | `NU` swept, `NX` held at 2, `NH` held at 5 | `NH` swept, `NX` held at 2, `NU` held at 1 |
|---:|---|---|---|
| 1 | `NV` 11: 20,568 / 20,568 | `NV` 17: 33,304 / 33,304 | `NV` 5: 14,688 / 14,688 |
| 2 | `NV` 17: 33,304 / 33,304 | `NV` 22: 42,320 / 46,680 | `NV` 8: 15,592 / 15,592 |
| 4 | `NV` 29: 70,648 / 73,080 | `NV` 32: 73,112 / 87,048 | `NV` 14: 25,608 / 25,608 |
| 8 | `NV` 53: 200,088 / 218,776 | `NV` 52: 173,368 / 211,288 | `NV` 26: 53,976 / 61,112 |
| 16 | `NV` 101: 518,376 / 759,320 | `NV` 92: 436,184 / 632,952 | `NV` 50: 161,816 / 196,168 |
| 20 / 23 / 32 | `NV` 125: 783,064 / 1,154,032 | `NV` 127: 806,304 / 1,190,616 | `NV` 98: 492,792 / 715,912 |
| 21 / 24 / 42 | `NV` 131: does not instantiate | `NV` 132: does not instantiate | `NV` 128: 825,008 / 1,209,176 |
| 43 | -- | -- | `NV` 131: does not instantiate |

**The first solve costs up to 1.48 times the steady-state solve**, and the ratio
climbs to about 1.47 at the top of all four lines while the two are equal at the
smallest configurations. **Construction costs about 4.3 times the steady-state
solve** at the top of every line.

### Supported maximum by task stack

Strict: no margin for the caller's own frames and none for RTOS overhead. The
ladders are geometric, so an entry names the largest MEASURED rung that fits and
says nothing whatever about the dimensions between rungs.

The entries come from the FIRST-SOLVE peak rather than the steady-state one,
because a caller that solves at all runs the first solve once. Only one cell
moves between the two readings: at 16 KiB the equal-dimension line reads 2 from
the steady-state figure and 1 from the first solve.

| task stack | largest fitting equal rung | largest fitting `NX`, `NU` at 1, `NH` at 5 | largest fitting `NU`, `NX` at 2, `NH` at 5 | largest fitting `NH`, `NX` at 2, `NU` at 1 |
|---|---|---|---|---|
| 4 KiB | none exists | none exists | none exists | none exists |
| 8 KiB | none exists | none exists | none exists | none exists |
| 16 KiB | 1 | no measured rung | no measured rung | 2 |
| 32 KiB | 2 | 1 | no measured rung | 4 |
| 48 KiB | 3 | 2 | 2 | 4 |
| 64 KiB | 3 | 2 | 2 | 8 |

**THE FIRST TWO ROWS SAY "NONE EXISTS" RATHER THAN "NO MEASURED RUNG", AND THE
DIFFERENCE IS THE POINT.** The smallest configuration this row has is
`NX = NU = NH = 1`, which induces `NV = 3`, and it needs 13,936 bytes. No
configuration of the predictive controller fits a 4 KiB or an 8 KiB task stack,
and that is a statement about every configuration rather than about the rungs
that were measured.

The shipped allocation proof's own configuration -- `NX = 2`, `NU = 1`,
`NH = 5`, so `NV = 17` and `MaxM = 12` -- needs 33,304 bytes in steady state and
the same again on its first solve. It wants a 48 KiB task and does not fit a
32 KiB one. This is the row where the heap statement and the stack statement
diverge most sharply: that solve allocates exactly nothing and still needs
nearly ten times the stack a `kalman_filter` needs at four states and four
outputs.

**If the controller is CONSTRUCTED on the same task stack the ladder above is
not the answer**, because construction is the deepest of the three chains. Read
against the construction peak: nothing fits 16 KiB; at 32 KiB only the smallest
equal rung (19,216 bytes) and a one-step horizon (24,688 bytes) fit; at 48 KiB
the equal line reaches rung 2 (44,576), the horizon line rung 2 (35,248) and the
state line rung 1 (47,872), and the input line still has nothing; at 64 KiB none
of the four moves further. Construct off the real-time task, or size the task
for construction rather than for the solve.

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
| 32 KiB | 17 | 17 |
| 48 KiB | 22 | 24 |
| 64 KiB | 26 | 25 |

**THE TWO COLUMNS DISAGREE AT 48 AND 64 KiB, AND THAT DISAGREEMENT IS THE
RESULT.** The peak is dominated by the decision dimension and is not a function
of it, so "the largest `NV` that fits" and "the `NV` below which everything fits"
are different numbers. At 64 KiB a configuration inducing `NV = 26` fits --
`NX = 1, NU = 24, NH = 1` -- while one inducing `NV = 25` does not:
`NX = 12, NU = 1, NH = 1` needs 65,848 bytes. Over the 44 decision dimensions the
fill reaches by more than one route, the deepest and shallowest configurations at
the same `NV` differ by up to **1.16x**, and the pattern is consistent: horizon
buys the decision dimension more cheaply than state does.

Size a task from the configuration, not from the decision dimension it induces.

### The three chosen dimensions INTERACT

The peak is dominated by the decision dimension and is not a function of it.

Asymptotically the whole-chain steady-state peak is about `50 bytes * NV * NV`:
over the seven measured points with `NV >= 92` the coefficient lies between 50.0
and 51.5 while `NX` ranges over 2 to 20 and `MaxM` over 12 to 120. At that end
the decision dimension predicts the peak to within three percent whichever
chosen dimension produced it.

It does not at moderate dimensions. `NX = 8, NU = 1, NH = 5` induces `NV = 53`
and needs 200,088 bytes; `NX = 2, NU = 8, NH = 5` induces `NV = 52` and needs
173,368. **A 1.9 percent difference in the decision dimension comes with a 15.4
percent difference in the peak**, and the same coefficient splits into two
groups there: 64.1 to 64.7 for the configurations with `NX = 2` and 70.7 to 71.2
for those with `NX >= 5`. So growth along one chosen dimension depends on where
the other two are held, and the dependence fades toward the ceiling rather than
growing.

The additive-separability test the estimator rows use refuses as well, and not
with one sign. Predicting the equal-dimension peak at rung `d` from the three
held lines as `W(d,1,5) + W(2,d,5) + W(2,1,d) - 2*W(2,1,5)` gives 1,952 against
13,936 measured at `d = 1`, 24,608 against 16,232 at `d = 2` and 102,760 against
91,256 at `d = 4` -- 7.14x, 0.66x and 0.89x. A cost that grows quadratically in a
sum of the three dimensions is not a sum of three one-dimensional costs, and the
three-axis null over-subtracts the fixed part at the small end.

The interior fill of 354 reachable triples describes that dependence rather than
only locating it, and the numbers above under the supported maximum are its
statement: at one decision dimension the peak varies by up to 1.16x with the
route taken to it. What the fill still does not resolve is anything between two
adjacent values of its stride of four.

### Provenance

`g++ (GNU) 16.1.1 20260728`, `-std=c++20 -O2 -fno-exceptions -fno-rtti -pthread`,
no `-march` (driver default `-mtune=generic -march=x86-64`), x86-64 Linux,
`double`, Eigen 3.4.1. HOST measurements. Runtime watermarks taken on a pthread
against a **zero-byte harness floor** and a 1,024-byte harness gap. The three
ladders were taken with a 4 MiB painted region on a 64 MiB thread stack; the
interior fill raises those to 32 MiB on a 512 MiB stack, which this row is the
reason for: its deepest construction chain is 4,135,224 bytes, 98.6 percent of
the 4 MiB window, so the fill would have been within 59,080 bytes of reporting
the window instead of the chain.

**THE BACKEND REVISION IS PART OF EVERY FIGURE IN THIS SECTION.** All of them
were measured against argmin at commit
`30a2cd4ea2fdab5b548080644908db640542eacd`, fetched from an empty tree at the
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
--diagonal` and `--held-dimension`, the fill with `--interior` and the two limits
with `--ceilings`. The driver configures a tree from empty to fetch the backend
when it does not already have one.

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
