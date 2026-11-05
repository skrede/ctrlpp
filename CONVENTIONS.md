# Coding conventions

The goal is code that is lean (small, readable and fast) with no loss of capability or speed. These
conventions are written for people, and they apply to everyone who writes ctrlpp, including the agents that
write much of it.

This file is ctrlpp's authoritative style specification. `.clang-format`, at the repository's root, records
the token-level shape. `EXCEPTIONS.md`, beside this file, is the only register of sanctioned over-limit
units. The checkers under `scripts/` enforce the mechanical rules, as the Enforcement section lists. Where
the short list of essentials in `CLAUDE.md` and this file differ, this file wins.

The rules serve readability. Where a rule that only review applies makes a particular spot harder to read,
readability wins and that spot is a fine exception. The rules a checker enforces yield only through
`EXCEPTIONS.md`, for size, or through the checker's own configuration, each entry with its reason.

## No planning artifacts

Never reference an ID, key or label produced by a planning or issue-tracking tool (phase numbers, milestone
names, plan or task numbers, requirement IDs, decision or invariant labels) anywhere in the codebase:
identifiers, comments, documentation, examples and commit messages. State the thing itself in plain words.
If a reader would have to open a planning document to decode a comment, the comment is wrong. The same goes
for narrative from how the code was built ("the slice", "the gate", "this milestone"), which means nothing to
a reader of the finished code.

## Language and portability

- ctrlpp is a header-only C++20 library. Its headers live under `lib/ctrlpp/include/ctrlpp/` in namespace
  `ctrlpp`, one folder per module: `control`, `dsp`, `estimation`, `lie`, `mhe`, `model`, `mpc`, `sysid`,
  `trajectory` and `util`, with the shared numerics in `detail/`.
- Eigen is the only required dependency. OSQP, NLopt and argmin are opt-in solver backends, each off by
  default and requested by its own CMake option; only the backend's adapter headers include it.
- The library and the examples throw nothing and use no run-time type identification. A fallible
  operation returns `ctrlpp::expected`, and every type is built through a factory that can report failure.
  The README's "Exception posture" section states what the test and benchmark trees build with.
- Use modern language features where they make the code clearer: concepts, `if constexpr`, designated
  initializers, ranges. A feature newer than C++20 is guarded by its feature-test macro with a fallback, or
  not used.
- The code compiles and runs on macOS, Linux and Windows, and on the embedded targets under
  `examples/embedded/`. Platform-specific code is kept isolated in the examples, never in the library.
- Follow the conventions of modern C++ libraries: the standard library, asio, boost. Existing ctrlpp
  conventions take precedence. When unsure, or when a popular idiom contradicts ctrlpp, ask rather than
  guess.

## Size

Lines include comments and blank lines; the budget is the unit as it sits on screen.

| Unit     | Target          | Hard ceiling |
| -------- | --------------- | ------------ |
| Function | 5 to 15 lines   | 25 lines     |
| File     | about 100 lines | 200 lines    |

Each scope does one thing (a library, a namespace, a class, a function), and the smaller that one thing is,
the cleaner the code reads. A unit over its ceiling is its one thing drawn too broad: decompose it, do not
widen the budget. Group code that is read together and split where a responsibility genuinely separates,
never to chase a number.

Avoid:

- **Salami-slicing**: carving a cohesive unit into fragments to fit a limit.
- **One-function files**: a file that exists only to hold a single function.
- **Artificial purity**: a wrapper or layer that adds indirection without separating a real
  responsibility.

Pieces split out of a unit go into the module's `detail/` folder (`ctrlpp/estimation/detail/`, for a
piece of an estimation header) and the `ctrlpp::detail` namespace. The public header stays where it is and
becomes the entry point that includes its pieces. A `detail/` file is still a real file and obeys the same
limits.

### How the ceilings are measured

A file's length is every line in it. A function's length is the count clang-tidy's
`readability-function-size` check reports: the body, from its opening brace to its closing brace. The
signature, the template header and the member initializer list fall outside the count; blank lines and
comments inside the body fall within it. Template patterns are measured whether or not anything
instantiates them, and code in an inactive `#if` branch is not measured. A CMake function or macro is
measured from its `function(` or `macro(` line to its `endfunction` or `endmacro` line.

### Scope

The hard ceilings apply to every tracked C++ file under `lib/` and to every tracked CMake file
(`CMakeLists.txt`, `*.cmake`, `*.cmake.in`) outside `tests/`, `examples/` and `benchmarks/`. In those three
trees the limits are guidelines. Shell, Python and other tooling scripts are exempt from the size limits;
the comment rule still applies to them.

### Exceptions

Decomposition is the presumption, even for a unit that reads as one piece. An exception is admitted only when
decomposing is measured to cost something: wall time, instruction count or stack use that rises, or numerical
output that gets less accurate, measured against a reference. Output that moves within rounding and is as
accurate as before is not a cost. Readability, or line-by-line comparability with a reference implementation
such as LAPACK or Reflexxes, is not a reason. An admitted unit is registered in `EXCEPTIONS.md`, and its
reason names the measured effect itself rather than pointing at a document that holds it. A registered unit
may grow to its recorded count plus ten percent before the gate fails.

There is no in-code marker and no debt baseline. A change is not done until every unit it touches is within
budget or registered.

## Comments

The codebase is written for people to read, and a comment must say something the code cannot. **The
default is no comment**: clear names and signatures carry the meaning. Write a comment only when it

- **(a)** cites a source: a paper, RFC, specification, URL, bug report or standards clause (an algorithm
  header cites the textbook section or paper it implements, by its key in `REFERENCES.bib`);
- **(b)** names a non-obvious algorithm or technique so a reader can look it up; or
- **(c)** explains a non-obvious design choice in one or two sentences: a bug workaround, a deliberate
  tradeoff, or a contract, ordering, side effect or invariant that is real but not visible from the code.

Never write comments that restate a name, signature, type or obvious control flow; tables that transcribe
the code beneath them; sign-posting or narration ("this section does X", "first we ... then we ...");
bookkeeping or process markers such as size-gate justifications, ledgers or "see the doc" pointers; or
design history ("relocated from ...", "X was the example").

Doc-comment blocks (`///`, `/** */`) are comments and are not exempt. Keep only text of the three kinds,
written as `//`; prose that documents an API belongs under `docs/api/`, not in the header.

A comment longer than a short tag goes on its own line above the code, never trailing; a brief trailing tag
on a data line, such as a unit or an enumerator's value, is fine. Keep a multi-line comment tight, with no
blank `//` lines padding it out. Comments count toward the size budget. CMake comments obey the same rule.

## Includes

Project includes (`#include "ctrlpp/..."`, rooted at `lib/ctrlpp/include`) come first, third-party
libraries (Eigen, the solver backends, Catch2 in the tests) second, the standard library third. The three
sections are separated by one blank line.

- Within a section, includes are grouped by folder, one blank line between groups, and never more than one
  blank line between any two groups.
- In the project section the file's own folder comes first, then the remaining project folders. The own
  group is the file's exact folder: for a header in `ctrlpp/estimation/detail/` it is that `detail/`
  folder, and `ctrlpp/estimation/` is one of the remaining groups.
- Within a folder group, includes sort by file-name length, shortest first; names of equal length sort
  alphabetically.

Order differs from this only where one header genuinely must precede another; each such case is listed,
with its reason, in the style lint's configuration. `.clang-format` does not sort includes, so the order is
kept by hand.

## Header guards and namespaces

- Header guards, not `#pragma once`. A library header's guard is `HPP_GUARD_` followed by its path under
  `lib/ctrlpp/include`, uppercased, with each `/` and the extension's `.` written as `_`:
  `ctrlpp/control/pid.h` is guarded by `HPP_GUARD_CTRLPP_CONTROL_PID_H`. An example header's guard is
  `HPP_GUARD_CTRLPP_EXAMPLES_` followed by its path under `examples/`, derived the same way:
  `examples/embedded/shared/dsp_demo.h` is guarded by `HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_DSP_DEMO_H`.
- No `// namespace X` comment after a closing namespace brace, and no comment after the guard's `#endif`.

## Naming

- Types, functions and variables are lowercase with underscores (`step_budget_solver`). Template
  parameters are `CamelCase` (`Scalar`, `Policy`), and dimension parameters are uppercase abbreviations
  (`NX`, `NU`, `NY`).
- Private and protected data members take an `m_` prefix (`m_state`, `m_updates_since_reset`). No name ends
  in a trailing underscore. The public fields of an aggregate take no prefix.
- `static constexpr` members are plain snake_case (`max_phases`), without `m_`.
- A member holding a callback or hook (a `std::function`, a `std::move_only_function`) takes a `_cb`
  suffix, whether it is a private member (`m_feed_forward_cb`) or a public field of an API type
  (`cost_cb`). Callbacks are grouped together; the suffix marks injected behavior, not state.

## Types and attributes

- **Spell out return types.** `auto` as a function's return type hides what a reader looks for first. A
  trailing return type (`auto f() -> T`) counts as an `auto` return: it stays only where the return type
  depends on the parameters (a `decltype` of the arguments) or cannot be written in front. A genuinely
  unspellable type, such as a closure, is the other exception. `auto` for a local is fine where the
  initializer makes the type obvious.
- **No `[[nodiscard]]`**, and no other attribute the compiler does not require. The attributes in use are
  `[[fallthrough]]`, `[[maybe_unused]]` and `[[noreturn]]`.
- **No decorative `constexpr`.** A function or variable is `constexpr` only when it can actually be
  evaluated in a constant expression.
- **Fixed-width integers** (`std::int32_t`, `std::uint16_t`, ...) rather than `int`, `long` and friends,
  and never a platform-specific alias where a fixed-width type serves. Horizons, iteration counts and
  phase counts are fixed-width too. The carve-outs:
  - a runtime dimension or index handed to or taken from Eigen is an `Eigen::Index`, the type Eigen itself
    uses;
  - a compile-time dimension template parameter keeps the type the API already uses for it: the
    `std::size_t` of ctrlpp's `NX`, `NU` and `NY`, or an `int`, as in Eigen's own templates;
  - a standard-container size is a `std::size_t`;
  - a C library interface keeps the types it declares: the `int` of the libm functions (the exponent of
    `std::frexp`), and the integer, pointer and opaque types of the OSQP and NLopt APIs (`OSQPInt`,
    `OSQPFloat`, NLopt's `unsigned` and `void*`).
- A fixed-size array is a `std::array`, not a C array.

## API design

- No raw pointers in public APIs. Pass by `const&`, by value, by `std::span`, or by smart pointer where that
  is the right ownership story.
- Distinguish **required** (no default), **required with a default** (the override is optional) and
  **`std::optional`** (absence is itself meaningful). A field with a literature default is a plain value
  carrying that default, never an `optional` standing in for one.

## Construction

- A type initializes its members in the **constructor's init list, never in-declaration**, so there is one
  source of truth and two defaults can never drift apart.
- The one exception is a configuration aggregate: a `struct` with no user-declared constructor, no private
  or protected members and no invariants, named `*_config`, `*_options` or `*_settings`, or a nested struct
  named `config` (as a policy's configuration is). It keeps its defaults in-declaration and stays an
  aggregate, so a caller designated-initializes only the fields it overrides and every other field keeps
  its default. `pid_config` and `nmpc_config` are of this kind.
- **Every other struct or class whose members carry defaults gets a constructor** that sets them in its
  init list. A pure aggregate with no defaults needs no constructor.
- A result or diagnostic type, a value a function fills in and hands back, gets a default constructor
  only; the code that produces it assigns its fields.
- A value type whose call sites construct it positionally (`return {a, b, c};`) gets a member-wise
  constructor taking every field in declaration order, plus a default constructor carrying the defaults,
  so positional construction keeps compiling and no field is left indeterminate. A trailing parameter of
  the member-wise constructor may take a default where callers commonly omit it, as the state-space pair's
  `C` and `D` default to zero. `path_point` is of this kind.
- The init list is written in member-declaration order, since members initialize in declaration order
  regardless (no `-Wreorder`). A member consumed by another member's constructor is declared before its
  consumer, and that dependency overrides the member triangle below.

## Class layout

Data comes before the operations on it:

```
class T
{
    // leading private nested helper types, if any

public:
    constructors
    public functions

private:
    data members
    private functions
};
```

The data members are ordered by the member triangle, subject to construction order. Each member function
is separated from the next by one blank line; the data members are packed together.

## Lifetimes and ownership

- Lifetimes are structural and single-owner; the owner sequences teardown.
- **Raw pointers are an anti-pattern.** `unique_ptr` for sole ownership, `shared_ptr` for genuinely shared
  ownership, and references for everything else: an owner hands its long-lived utilities out by reference,
  through a context object the parts reach. A raw pointer survives only where absence is a real, modeled
  state, such as a nullable lookup.
- A non-owning relationship that must outlive its holder is a reference, stored as
  `std::reference_wrapper<T>` when it lives in a container that needs assignability.
- No `shared_from_this` or `enable_shared_from_this`, and no liveness flags guarding a posted task against
  an expired target; fix the ownership and teardown order instead.

## Formatting

`.clang-format` owns the token-level shape, including the column limit. Every file a change touches is
formatted with clang-format 18 before the change is committed.

### Shape

Where `.clang-format` leaves the choice open, these are judgment calls for the eye, not mechanical rules.

- **Downward waterfall.** When a signature or declaration wraps, each continuation line is at least as long
  as the one above it, so the wrapped arguments slope down and the eye falls to the next one.
- **Member triangle.** A class's data members are ordered by declaration-line length, shortest first, so
  they form a downward slope. A single space precedes each name; names are not aligned into a column.
  Semantic grouping wins where it helps, such as keeping the callbacks together.
- **Logical blank lines.** A blank line inside a function separates its logical parts, the way a paragraph
  break does.

## Enforcement

The checkers run from the repository's root and measure tracked files only, since their file lists come
from `git ls-files`.

- **The size gates.** `python3 scripts/check_file_size.py` measures file lengths and CMake function and
  macro lengths. `python3 scripts/check_function_size.py` measures every C++ function through clang-tidy
  18 (`--clang-tidy`, or `CLANG_TIDY`, names the binary; `--eigen-include` the Eigen headers); the
  `CTRLPP_TIDY_EXTRA_ARG` variable carries the backend include paths and defines, since the measured
  translation unit includes every tracked header, the backend adapters among them. Both check every unit
  in scope against the ceilings and against `EXCEPTIONS.md`, and each kind of failure has its own exit
  code, listed by `--help`. `bash scripts/size_gate_negative_control.sh [files|functions|all]` proves that
  each gate can fail, in a throwaway copy of the tree.
- **The style lint.** `python3 scripts/check_style.py [--rule ID ...] [PATH ...]` lints the tracked C++
  files under `lib/` and `examples/`. It checks include order; attributes; `auto` and trailing return
  types; in-declaration member initializers outside the configuration exception; member naming (`m_`, no
  trailing underscore, `_cb`); header guards and `#pragma once`; comments after a closing namespace brace
  or the guard's `#endif`; and American spelling, which it also checks in the CMake files. It flags a
  call made as a bare statement that discards a `ctrlpp::expected`, which no attribute guards.
- **Source purity.** `bash scripts/source_purity_check.sh` holds the library and the examples to the
  error policy above.
- **The doc snippet build.** `python3 scripts/check_doc_snippets.py [PATH ...]` compiles, syntax-only,
  every fenced `cpp` block that contains an `#include` in the tracked `docs/` pages and `README.md`, each as
  a translation unit of its own; a failure is reported at its page and line. A deliberately partial block
  carries `<!-- snippet: fragment -->` on the line before its fence and is skipped, and a block that needs
  exceptions carries `<!-- snippet: exceptions -->` there. The `CTRLPP_TIDY_EXTRA_ARG` variable carries the
  backend include paths and defines, as for the function-size gate, and
  `python3 scripts/check_doc_snippets_test.py` proves the build can fail.

A rule's legitimate carve-outs (a header that must precede another, a quoted title in British spelling, the
list of fallible calls) live in the lint's configuration, `scripts/style_config.py`, each with its reason.
Comments, class layout and raw-pointer use are review-only: a convention no tool enforces still binds the
change, and a person rather than a build catches it.

## English

American spelling: "stabilizing", not "stabilising"; "center", not "centre".

## Commit messages

```
{Prefix}: {summary sentence}.

- {what was done, one line per item}
- {another item if applicable}
```

The summary line is brief and descriptive, and the bullets expand on it; a single-item commit may omit
them. The prefixes:

| Prefix          | When to use |
|-----------------|-------------|
| `Feature:`      | New user-visible capability or API surface. |
| `Fix:`          | Bug fix or correctness-affecting change. |
| `Refactor:`     | Internal reorganization with no user-visible behavior change. |
| `Build:`        | Build system, CI and tooling. |
| `Docs:`         | Documentation only. |
| `Examples:`     | Changes to the examples only. |
| `Optimization:` | Performance change with no API or correctness impact. |
| `Test:`         | Tests only. |
| `WIP:`          | A commit whose code does not compile. |

The rule against planning artifacts applies to commit messages too.
