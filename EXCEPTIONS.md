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
| lib/ctrlpp/include/ctrlpp/estimation/manifold_ukf.h | 344 | Decomposing the filter by moving its predict and update algebra into detail function templates, declared inline and taking the state by reference, raises the worst-path stack of every measured update. The double filter with three outputs rises from 3448 to 3664 bytes on Cortex-M7 at -O2 and from 4048 to 4352 bytes on the ESP32 (Xtensa, -Os), and its predict from 2240 to 2352 bytes on the ESP32. Keeping the measurement step inside the update and carrying its rollback state in one struct still raises the update from 3448 to 3544 bytes on Cortex-M7 at -O2 and from 4048 to 4240 bytes on the ESP32. |
| lib/ctrlpp/include/ctrlpp/estimation/particle_filter.h | 641 | Every decomposition of the filter into weight, likelihood, sampling and estimate components that was measured raises the worst-path stack of one per-step update. Merging the log-weight and linear-weight resamplings into one function gives a linear-weight update the log-weight frame: for 100 double particles it rises from 1848 to 2664 bytes on the host at -O3 and from 1440 to 2240 bytes on the ESP32 (Xtensa, -Os). Keeping one resampling per weight mode inside the weight component raises every host update by 16 to 48 bytes, the 100-particle double one from 2648 to 2664, and the multinomial one from 3056 to 3136 bytes on Cortex-M7 at -O3. Keeping the per-mode resamplings as filter members instead raises the 10-particle log-weight update from 472 to 712 bytes on the host and from 600 to 800 bytes on Cortex-M7 at -O3, and the multinomial one from 3472 to 3584 bytes on the host. |
| lib/ctrlpp/include/ctrlpp/estimation/ukf.h | 327 | Decomposing the filter by moving its predict and update algebra into detail function templates, declared inline and taking the state by reference, raises the worst-path stack of every measured update. The double filter with four states and two outputs rises from 3416 to 3880 bytes on Cortex-M7 at -O2 and from 4976 to 5136 bytes on the ESP32 (Xtensa, -Os), and the float one from 4920 to 5176 bytes on the host at -O3. Keeping the measurement step inside the update and carrying its rollback state in one struct still raises the double update from 3416 to 3792 bytes on Cortex-M7 at -O2 and the float one from 4920 to 5224 bytes on the host: GCC keeps the free gain function out of line where it inlines the member, which alone adds 312 bytes on Cortex-M7 at -O2, and the struct adds 56 bytes to the update frame. |

## Functions

| Path | Function | Lines at registration | Why |
| --- | --- | --- | --- |
| lib/ctrlpp/include/ctrlpp/control/mrac.h | evaluate | 27 | Every decomposition that brings it within the function ceiling raises the worst-path stack of one call. With the parameter update in its own helper, which GCC keeps out of line, the double, two-state dead-zone controller rises from 296 to 400 bytes on Cortex-M7 at -O3 and from 672 to 928 bytes on the ESP32 (Xtensa, -Os), and the e-modification controller from 528 to 608 bytes on Cortex-M7 at -O2. With one leakage factor shared by the sigma- and e-modification updates instead, the e-modification controller rises from 304 to 320 bytes on the host. Forming the command after the update, which the remaining lines would need, raises the single-state controller from 200 to 208 bytes on Cortex-M7 at -Os. Splitting off only the dead-zone gate and the finiteness check and commit leaves every measured root unchanged. |
| lib/ctrlpp/include/ctrlpp/estimation/detail/geodesic_mean.h | geodesic_mean_impl | 28 | Every split that brings it within the function ceiling raises the stack of the manifold UKF predict, which runs it every step. With the tangent-space average in its own function template, inline or not, the double predict with three outputs rises from 2292 to 2364 bytes on Cortex-M7 at -Os and from 2240 to 2304 bytes on the ESP32 (Xtensa, -Os), and the float one from 1216 to 1264 bytes on the ESP32; sharing the hemisphere flip with the covariance as well raises the float one to 1328 bytes. Writing the flip as one vector negation instead raises the function alone from 416 to 432 bytes on the host at -O3. |
| lib/ctrlpp/include/ctrlpp/sysid/rls.h | update | 42 | Every decomposition that brings it within the function ceiling raises the worst-path stack of one update. With the denominator classification and the parameter and covariance step in their own helpers, which GCC keeps out of line at -Os, the ESP32 (Xtensa, -Os) frame rises by 48 to 96 bytes in every measured estimator, the two-parameter double one from 432 to 528 bytes, and the two-parameter float one rises from 256 to 336 bytes on Cortex-M7 at -Os. With only the operand checks and the denominator classification split off, the two-parameter float estimator rises from 232 to 264 bytes on the host at -O3 and the four-parameter double one from 33552 to 33560 bytes on Cortex-M7 at -O2; splitting off only the operand checks costs the same 32 bytes on the host. |
