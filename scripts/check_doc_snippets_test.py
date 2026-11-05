#!/usr/bin/env python3
"""Fixture repositories proving the doc snippet build fails on a broken or throwing snippet and on a bad marker,
skips what it should, and reads tracked pages only."""

import os
import subprocess
import sys
import tempfile
from pathlib import Path

CHECKER = Path(__file__).resolve().parent / "check_doc_snippets.py"
PROBE = {"lib/ctrlpp/include/ctrlpp/probe.h": "inline int probe(int x)\n{\n    return x + 1;\n}\n"}
GOOD = '#include "ctrlpp/probe.h"\n\nint main()\n{\n    return probe(1);\n}\n'
BROKEN = '#include "ctrlpp/probe.h"\n\nint main()\n{\n    int* p = probe(1);\n    return 0;\n}\n'
THROWS = '#include <stdexcept>\n\nint main()\n{\n    throw std::runtime_error("x");\n}\n'
NO_INCLUDE = "int main()\n{\n    int* p = 1.5;\n}\n"


def page(body: str, marker: str = "", language: str = "cpp") -> str:
    return f"# x\n\n{marker + chr(10) if marker else ''}```{language}\n{body}```\n"


def git(root: Path, *args: str) -> None:
    subprocess.run(["git", "-C", str(root), *args], check=True, capture_output=True)


def build(root: Path, files: dict[str, str], untracked: dict[str, str] | None = None) -> None:
    git(root, "init", "-q")
    for name, content in files.items() | (untracked or {}).items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    for name in files:
        git(root, "add", name)


def run_checker(root: Path, *args: str) -> tuple[int, str]:
    environment = {name: value for name, value in os.environ.items() if name != "CTRLPP_TIDY_EXTRA_ARG"}
    compiler = os.environ.get("CXX") or "c++"
    done = subprocess.run([sys.executable, str(CHECKER), "--source-root", str(root), "--compiler", compiler, *args],
                          capture_output=True, text=True, env=environment)
    return done.returncode, done.stdout + done.stderr


CASES = [
    ("a snippet that includes the library and uses it compiles", PROBE | {"docs/x.md": page(GOOD)}, (), 0,
     ("checked 1 snippets in 1 pages: 0 failed",)),
    ("a type error fails at the fence and the diagnostic names the page", PROBE | {"docs/x.md": page(BROKEN)}, (), 1,
     ("error: docs/x.md:3: does not compile", "docs/x.md:8:")),
    ("a marked fragment is skipped", PROBE | {"docs/x.md": page(BROKEN, "<!-- snippet: fragment -->")}, (), 0,
     ("1 fragments skipped", "checked 0 snippets")),
    ("a broken block without an include is not a snippet",
     PROBE | {"docs/x.md": page(GOOD) + page(NO_INCLUDE)}, (), 0, ("checked 1 snippets",)),
    ("a broken block in another language is not a snippet",
     PROBE | {"docs/x.md": page(GOOD) + page(BROKEN, language="c")}, (), 0, ("checked 1 snippets",)),
    ("an unmarked throwing snippet fails under the default flags", {"docs/x.md": page(THROWS)}, (), 1,
     ("error: docs/x.md:3: does not compile",)),
    ("a throwing snippet marked for exceptions compiles",
     {"docs/x.md": page(THROWS, "<!-- snippet: exceptions -->")}, (), 0, ("1 compiled with exceptions",)),
    ("an unknown marker is an error", PROBE | {"docs/x.md": page(GOOD, "<!-- snippet: fragmnet -->")}, (), 3,
     ("error: docs/x.md:3: unknown marker",)),
    ("a valid marker before a block without an include is an error",
     PROBE | {"docs/x.md": page(GOOD) + page(NO_INCLUDE, "<!-- snippet: fragment -->")}, (), 3,
     ("marker stands before a block that is not a snippet",)),
    ("a marker that no fence follows is an error",
     PROBE | {"docs/x.md": page(GOOD) + "<!-- snippet: fragment -->\n\ntext\n"}, (), 3,
     ("does not stand before a fenced block",)),
    ("a broken snippet in the root README is found", PROBE | {"README.md": page(BROKEN)}, (), 1,
     ("error: README.md:3: does not compile",)),
    ("a page outside docs/ is not read", PROBE | {"docs/x.md": page(GOOD), "notes/y.md": page(BROKEN)}, (), 0,
     ("checked 1 snippets in 1 pages",)),
    ("a path argument checks only the pages it names",
     PROBE | {"docs/a/x.md": page(GOOD), "docs/b/y.md": page(BROKEN)}, ("docs/a",), 0,
     ("checked 1 snippets in 1 pages",)),
    ("a path that names no tracked page is a usage error", PROBE | {"docs/x.md": page(GOOD)}, ("docs/gone.md",), 2,
     ("names no tracked documentation page",)),
    ("--list prints each snippet and its kind and compiles nothing",
     PROBE | {"docs/x.md": page(BROKEN) + page(THROWS, "<!-- snippet: exceptions -->")}, ("--list",), 0,
     ("docs/x.md:3 compile\n", "docs/x.md:15 exceptions\n")),
    ("a compiler that cannot be executed", PROBE | {"docs/x.md": page(GOOD)}, ("--compiler", "/nonexistent/c++"), 7,
     ("cannot execute the compiler",)),
    ("a repository with no snippet has nothing to check", {"docs/x.md": page(NO_INCLUDE)}, (), 6,
     ("holds a snippet",)),
    ("--help lists every exit code and the backend variable", {}, ("--help",), 0,
     ("  0  every", "  1  a snippet", "  2  the arguments", "  3  a marker", "  6  the tracked", "  7  the compiler",
      "CTRLPP_TIDY_EXTRA_ARG")),
]


def check(label: str, code: int, output: str, want_code: int, want_text: str) -> bool:
    if code == want_code and want_text in output:
        return True
    print(f"FAIL {label}: exit {code} (want {want_code}), expected {want_text!r} in:\n{output}")
    return False


def run_cases(root: Path) -> int:
    failures = 0
    for index, (label, files, args, want_code, wanted) in enumerate(CASES):
        case = root / f"case-{index}"
        case.mkdir()
        build(case, files)
        code, output = run_checker(case, *args)
        failures += not all(check(label, code, output, want_code, text) for text in wanted)
    return failures


def run_special(root: Path) -> int:
    untracked = root / "untracked"
    untracked.mkdir()
    build(untracked, PROBE | {"docs/x.md": page(GOOD)}, {"docs/y.md": page(BROKEN)})
    failures = not check("an untracked page is ignored", *run_checker(untracked), 0, "checked 1 snippets in 1 pages")
    empty = root / "empty"
    empty.mkdir()
    build(empty, {})
    return failures + (not check("an empty repository has nothing to check", *run_checker(empty), 6,
                                 "holds a snippet"))


def main() -> int:
    with tempfile.TemporaryDirectory() as workspace:
        root = Path(workspace)
        failures = run_cases(root) + run_special(root)
    total = len(CASES) + 2
    if failures:
        print(f"{failures} of {total} cases failed")
        return 1
    print(f"{total} cases passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
