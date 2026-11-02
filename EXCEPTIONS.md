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
