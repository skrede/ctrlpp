# Size exceptions

Files are held to about one hundred lines and two hundred at most, counting every line. Functions are
held to five to fifteen lines and twenty-five at most. This file is the only register of the units
sanctioned to stay over those ceilings, each with its reason, and two checkers enforce it.

A row belongs here only when decomposing the unit was measured to make the code slower or worse. The
first answer to an over-limit unit is to decompose it. There is no in-code marker and no other list.

## What the checker enforces

`scripts/check_file_size.py` measures the file lengths, and the lengths of CMake functions and macros.
It fails on:

- an over-limit unit that is not registered;
- a stale row, whose unit is at or under the ceiling or no longer exists;
- growth: a registered unit past its recorded count plus ten percent;
- an inflated row, which records more lines than its unit has and must be re-recorded down;
- a malformed register.

Each failure has its own exit code and message, because each has a different fix. A malformed row is
a failure, not a skip. There are no silent exceptions.

A C++ function's length is the line count clang-tidy's `readability-function-size` check reports: the
body from its opening brace to its closing brace, so the signature and the member initializer list
are not counted. `scripts/check_function_size.py` measures it over every tracked header and fails on
the same faults for functions. A CMake function or macro is measured from its `function(` or `macro(`
line to its `endfunction` or `endmacro` line, by the same span convention.

## Scope

The checkers measure every file under `lib/` that is tracked in version control and has a C++
extension, and every tracked CMake file (`CMakeLists.txt`, `*.cmake`, `*.cmake.in`) outside `tests/`,
`examples/` and `benchmarks/`. In those three trees the sizes are guidelines and are not measured.
Untracked files are invisible to the checkers, since the file list comes from `git ls-files`.

## Reasons

An exception is admitted only when decomposing the unit is measured to cost something: stack use,
instruction count or wall time that rises, or numerical output that moves. The reason states that
measured cost itself, in words and numbers, and names what splitting the unit would harm. Readability,
or line-by-line comparability with a reference implementation, is not a reason.

A reason that only says the unit is large, or cohesive, is rejected. A reason that restates the line
count is rejected, because the count has its own column. A reason that promises a future split is
rejected too: this file registers decisions already taken, not work outstanding. A reason holds the
measurement itself and never points at another document for it.

Each row records the unit's length when it was registered. A registered unit may grow to that count
plus ten percent, rounded up, before the checker fails. When a change exceeds that, re-record the
count deliberately or decompose the unit.

## Files

| Path | Lines at registration | Why |
| --- | --- | --- |

## Functions

| Path | Function | Lines at registration | Why |
| --- | --- | --- | --- |
| lib/ctrlpp/include/ctrlpp/control/mrac.h | evaluate | 27 | Every decomposition that brings it within the function ceiling raises the worst-path stack of one call. With the parameter update in its own helper, which GCC keeps out of line, the double, two-state dead-zone controller rises from 296 to 400 bytes on Cortex-M7 at -O3 and from 672 to 928 bytes on the ESP32 (Xtensa, -Os), and the e-modification controller from 528 to 608 bytes on Cortex-M7 at -O2. With one leakage factor shared by the sigma- and e-modification updates instead, the e-modification controller rises from 304 to 320 bytes on the host. Forming the command after the update, which the remaining lines would need, raises the single-state controller from 200 to 208 bytes on Cortex-M7 at -Os. Splitting off only the dead-zone gate and the finiteness check and commit leaves every measured root unchanged. |
| lib/ctrlpp/include/ctrlpp/estimation/detail/geodesic_mean.h | geodesic_mean_impl | 28 | Every split that brings it within the function ceiling raises the stack of the manifold UKF predict, which runs it every step. With the tangent-space average in its own function template, inline or not, the double predict with three outputs rises from 2292 to 2364 bytes on Cortex-M7 at -Os and from 2240 to 2304 bytes on the ESP32 (Xtensa, -Os), and the float one from 1216 to 1264 bytes on the ESP32; sharing the hemisphere flip with the covariance as well raises the float one to 1328 bytes. Writing the flip as one vector negation instead raises the function alone from 416 to 432 bytes on the host at -O3. |
| lib/ctrlpp/include/ctrlpp/estimation/ekf.h | update | 28 | Every split that brings it within the function ceiling raises the worst-path stack of one update. With the rollback of the saved state, covariance, innovation and NIS in its own member, the float filter with four states and two outputs rises from 2672 to 2688 bytes on the ESP32 (Xtensa, -Os), and the double one with analytic Jacobians from 3376 to 3392. With the measurement step from the innovation to the NIS in its own member, the double filter with one output rises from 876 to 948 bytes on Cortex-M7 at -Os and from 1024 to 1088 bytes on the ESP32, and the double one with four states from 2376 to 2424 bytes on Cortex-M7 at -O3. With the innovation in its own member, the float filter rises from 2696 to 2744 bytes on Cortex-M7 at -Os; with the restore written as one tuple assignment, the double filter rises from 2824 to 2840 bytes on Cortex-M7 at -O2. |
| lib/ctrlpp/include/ctrlpp/estimation/mekf.h | update | 29 | Every split that brings it within the function ceiling raises the worst-path stack of one update. With the rollback of the saved attitude, bias, covariance, innovation and state in its own member, the double filter with three outputs rises from 3824 to 3872 bytes on Cortex-M7 at -O2 and from 4464 to 4496 bytes on the ESP32 (Xtensa, -Os), and the one with one output from 1792 to 1888 bytes on Cortex-M7 at -O3. With the measurement step from the innovation to the state refresh in its own member, the one-output filter rises from 2368 to 2576 bytes on the host at -O3 and from 2432 to 2552 bytes on Cortex-M7 at -O2. With the restore written as one tuple assignment, it rises from 1792 to 1880 bytes on Cortex-M7 at -O3. |
| lib/ctrlpp/include/ctrlpp/mpc/argmin_problem.h | constraint_jacobian | 42 | Every decomposition that brings it within the function ceiling raises the worst-path stack of one constraint-Jacobian evaluation, which argmin's SQP policies run every iteration. With the row gather in a free function template and the central-difference loop in a member, the dynamic bridge rises from 152 to 184 bytes on Cortex-M7 at -O2 and from 200 to 232 bytes on the host at -O3, and the bridge with seventeen variables and twelve constraint rows from 152 to 192 bytes on Cortex-M7 at -O2. With both in members, the dynamic bridge rises from 152 to 176 bytes on Cortex-M7 at -O2. With only the central-difference loop split off, in a member it rises from 200 to 248 bytes on Cortex-M7 at -Os and from 256 to 304 bytes on the ESP32 (Xtensa, -Os), and in a free function template to 336 bytes on the ESP32. |
| lib/ctrlpp/include/ctrlpp/mpc/detail/argmin_options.h | map_argmin_status | 26 | Every split that brings it within the function ceiling raises the worst-path stack of one solve, which maps argmin's status on every return. With the inexact and failed statuses in a second switch behind a default, the unconstrained L-BFGS-B solve rises from 17168 to 17184 bytes on the host at -O3, and the two-state NMHE update from 7936 to 7944 bytes on Cortex-M7 at -O3. |
| lib/ctrlpp/include/ctrlpp/mpc/detail/argmin_options.h | map_argmin_stop_criterion | 28 | Every split that brings it within the function ceiling raises the worst-path stack of one solve, which maps argmin's status on every return. With the limit and remaining stops in a second switch behind a default, the two-state NMHE update rises from 7936 to 7944 bytes on Cortex-M7 at -O3. Splitting both status mappings that way also raises the fixed-dimension predictive solve from 33408 to 33456 bytes on the host at -O3 and the SLSQP step from 8112 to 8160 bytes on the ESP32 (Xtensa, -Os). |
| lib/ctrlpp/include/ctrlpp/sysid/rls.h | update | 42 | Every decomposition that brings it within the function ceiling raises the worst-path stack of one update. With the denominator classification and the parameter and covariance step in their own helpers, which GCC keeps out of line at -Os, the ESP32 (Xtensa, -Os) frame rises by 48 to 96 bytes in every measured estimator, the two-parameter double one from 432 to 528 bytes, and the two-parameter float one rises from 256 to 336 bytes on Cortex-M7 at -Os. With only the operand checks and the denominator classification split off, the two-parameter float estimator rises from 232 to 264 bytes on the host at -O3 and the four-parameter double one from 33552 to 33560 bytes on Cortex-M7 at -O2; splitting off only the operand checks costs the same 32 bytes on the host. |
| lib/ctrlpp/include/ctrlpp/trajectory/modified_trap_trajectory.h | eval_first_half | 46 | Every split that brings it within the function ceiling raises the worst-path stack of one `evaluate`, which runs it for every sample. With each of its three phases in a private member taking the shared constants and writing position, velocity and acceleration through references, the double evaluate rises from 88 to 96 bytes on Cortex-M7 at -O3, because GCC saves one more floating-point register, and the float evaluate from 80 to 96 bytes on the host at -O3. Keeping the ramp-up inline, recomputing T/8 in each member, or reordering the data members still raises the Cortex-M7 -O3 double evaluate from 88 to 96 bytes. Members returning the three values raise the Cortex-M7 -O2 evaluates by 8 bytes and the ESP32 (Xtensa, -Os) float evaluate from 112 to 128. Calling the split once on the mirrored time raises the Cortex-M7 -O3 double evaluate to 104 bytes. |
| lib/ctrlpp/include/ctrlpp/trajectory/detail/trapezoidal_rescale.h | solve_rescale | 140 | Every split that brings it within the function ceiling raises the worst-path stack of one retiming, which `rescale_to` and `can_rescale_to` run for every axis a synchronization retimes, or moves its output. With the valley shape's decrement solve in its own member, the float retiming rises from 192 to 200 bytes on Cortex-M7 at -Os and from 272 to 288 bytes on the ESP32 (Xtensa, -Os), and the float `can_rescale_to` from 88 to 96 bytes on Cortex-M7 at -O3. With the request checks split from the shape solve, which GCC then keeps out of line, the double `can_rescale_to` rises from 112 to 224 bytes on the host at -O3 and from 160 to 224 bytes on Cortex-M7 at -O3. With the three boundary shapes behind one dispatching member, the float retiming rises from 272 to 320 bytes on the ESP32. With the plateau's or the valley's cruise-velocity form in its own member, clang at -ffp-contract=fast keeps that member out of line and fuses a * h there, where the function shares it with the peak test, which changes the retimed profile of 2489 of 40000 seeded commands. |
