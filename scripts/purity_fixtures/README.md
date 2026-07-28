# Source purity gate fixtures

This directory holds one deliberately non-conforming file per rule of
`scripts/source_purity_check.sh`, so the gate can be shown to fail. A gate with
no evidence that it can fail certifies nothing while looking like enforcement,
which is the defect class the gate itself exists to remove.

Every file here violates its rule on purpose. None is compiled, none belongs to
a build target, and no `add_subdirectory`, `add_executable` or `add_library`
anywhere in the repository names this directory. The layout mirrors the gate's
scan roots (`lib/` here stands for the repository's `lib/`) so the script's
existing scan-root argument reaches the fixtures without any special-casing,
and the directory sits under `scripts/`, which no rule scans, so a real run
never sees it. That placement is deliberate: if the fixtures were ever reachable
from a real scan root the gate would fail permanently against its own test data,
and the likely response would be to weaken the gate.

| File | Rule it violates |
|---|---|
| `lib/fixture_unwrap.h` | 1, the unwrap accessor |
| `lib/fixture_throw.h` | 2, a throw statement |
| `lib/fixture_rtti.h` | 3, run-time type identification |
| `lib/fixture_missing_class_annotation.h` | 4a, a manifest declaration missing its class-level attribute |
| `lib/fixture_site_annotation.h` | 4b, a site-level attribute on a trivial getter |
| `lib/fixture_discarded_fallible.h` | 4c, a bare-statement call to a fallible function |
| `lib/fixture_exception_gate.h` | 5, the exception-mode macro outside its two files |
| `CMakeLists.txt` | 6, the discard-warning promotion absent |
| `.clang-tidy` | 7, the attribute-inserting analysis check enabled |

Run the self-test with:

```sh
bash scripts/source_purity_check.sh --self-test
```

It asserts a nonzero detection count and, more importantly, that each rule
reports **exactly one** violation. The exact count is what makes this a real
self-test rather than a check that something went wrong: a script that crashed
before scanning would also exit nonzero, and only the count distinguishes the
two. It is also what catches the removal of a comment filter, because the
fixtures for rules 1, 2, 3 and 4b each carry the banned pattern once in a
comment as well as once in code.
