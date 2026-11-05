#!/usr/bin/env python3
"""Compile the C++ examples of the tracked documentation (the pages under docs/ and the root README.md),
syntax-only. A fenced cpp block that contains an #include line is a snippet and is compiled as a translation
unit of its own; a block without one is not checked. A deliberately partial block carries
<!-- snippet: fragment --> on the line before its opening fence and is skipped, and a block that needs
exceptions carries <!-- snippet: exceptions --> there and is compiled with -fexceptions."""

import argparse
import os
import posixpath
import re
import shlex
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import NamedTuple

EXTRA_ARG_VARIABLE = "CTRLPP_TIDY_EXTRA_ARG"
MARKERS = ("fragment", "exceptions")

EXIT_FAILED = 1
EXIT_USAGE = 2
EXIT_MARKER = 3
EXIT_NO_INPUT = 6
EXIT_NO_COMPILER = 7

EXIT_HELP = f"""exit codes:
  0  every snippet compiles
  {EXIT_FAILED}  a snippet does not compile
  {EXIT_USAGE}  the arguments are wrong, or a path names no tracked page
  {EXIT_MARKER}  a marker is unknown or stands before a block that is not a snippet (wins over {EXIT_FAILED})
  {EXIT_NO_INPUT}  the tracked files cannot be listed, a page cannot be read, or no page holds a snippet
  {EXIT_NO_COMPILER}  the compiler cannot be executed

{EXTRA_ARG_VARIABLE}, when set, is split shell-style into compiler arguments appended to the defaults; it
carries the backend include paths and defines.
"""

FENCE = re.compile(r"^( *)(`{3,}|~{3,})(.*)$")
INCLUDE = re.compile(r"^\s*#\s*include\b")
MARKER = re.compile(r"^\s*<!--\s*snippet\b(.*?)-->\s*$")
MARKER_TEXT = re.compile(r"^:\s*(\S+)\s*$")


Finding = tuple[str, int, str]


class Block(NamedTuple):
    page: str
    fence: int
    language: str
    body: list[str]
    marker: str | None


class Snippet(NamedTuple):
    block: Block
    kind: str


def tracked_pages(root: Path) -> list[str]:
    listing = subprocess.run(["git", "-C", str(root), "ls-files", "-z"], check=True,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout
    names = listing.decode("utf-8", "surrogateescape").split("\0")
    return sorted(name for name in names if name == "README.md" or
                  (name.startswith("docs/") and name.endswith(".md")))


def closing(lines: list[str], start: int, opening: re.Match) -> int:
    fence = opening.group(2)
    pattern = re.compile(rf"^ *{re.escape(fence[0])}{{{len(fence)},}}\s*$")
    return next((index for index in range(start + 1, len(lines)) if pattern.match(lines[index])), len(lines))


def dedented(lines: list[str], indent: int) -> list[str]:
    return [line[min(indent, len(line) - len(line.lstrip(" "))):] for line in lines]


def parse(page: str, text: str) -> tuple[list[Block], dict[int, str]]:
    lines, blocks, markers, index = text.splitlines(), [], {}, 0
    while index < len(lines):
        opening = FENCE.match(lines[index])
        if not opening:
            marker = MARKER.match(lines[index])
            if marker:
                markers[index + 1] = marker.group(1)
            index += 1
            continue
        end = closing(lines, index, opening)
        words = opening.group(3).split()
        body = dedented(lines[index + 1:end], len(opening.group(1)))
        blocks.append(Block(page, index + 1, words[0].lower() if words else "", body, markers.pop(index, None)))
        index = end + 1
    return blocks, markers


def is_candidate(block: Block) -> bool:
    return block.language == "cpp" and any(INCLUDE.match(line) for line in block.body)


def classify(page: str, blocks: list[Block], stray: dict[int, str]) -> tuple[list[Snippet], list[Finding]]:
    snippets = []
    errors = [(page, line, "the marker does not stand before a fenced block") for line in stray]
    for block in blocks:
        text = MARKER_TEXT.match(block.marker) if block.marker is not None else None
        kind = "compile" if block.marker is None else text.group(1) if text else None
        if kind not in MARKERS + ("compile",):
            errors.append((page, block.fence - 1, f"unknown marker '{block.marker.strip()}'; the markers are "
                                                  f"'snippet: fragment' and 'snippet: exceptions'"))
        elif is_candidate(block):
            snippets.append(Snippet(block, kind))
        elif kind != "compile":
            errors.append((page, block.fence - 1, f"the '{kind}' marker stands before a block that is not a "
                                                  f"snippet (a cpp block with an #include)"))
    return snippets, errors


def scan(root: Path, pages: list[str]) -> tuple[list[Snippet], list[Finding]]:
    snippets, errors = [], []
    for page in pages:
        found, wrong = classify(page, *parse(page, (root / page).read_bytes().decode("utf-8", "surrogateescape")))
        snippets += found
        errors += wrong
    return snippets, errors


def flags(args: argparse.Namespace, kind: str) -> list[str]:
    exceptions = "-fexceptions" if kind == "exceptions" else "-fno-exceptions"
    include = args.source_root.resolve() / "lib" / "ctrlpp" / "include"
    return ["-std=c++20", "-fsyntax-only", exceptions, "-fno-rtti", "-I", str(include), "-isystem",
            str(args.eigen_include), *args.extra_args]


def translation_unit(snippet: Snippet) -> str:
    page = snippet.block.page.replace("\\", "\\\\").replace('"', '\\"')
    return f'#line {snippet.block.fence + 1} "{page}"\n' + "".join(line + "\n" for line in snippet.block.body)


def compile_snippet(args: argparse.Namespace, directory: Path, number: int,
                    snippet: Snippet) -> tuple[int | None, str]:
    path = directory / f"snippet-{number}.cpp"
    path.write_text(translation_unit(snippet), encoding="utf-8", errors="surrogateescape")
    try:
        done = subprocess.run([args.compiler, *flags(args, snippet.kind), str(path)], cwd=args.source_root,
                              capture_output=True, text=True, errors="replace")
    except OSError as exc:
        return None, str(exc)
    return done.returncode, done.stdout + done.stderr


def compile_all(args: argparse.Namespace, snippets: list[Snippet]) -> list[tuple[int | None, str]]:
    with tempfile.TemporaryDirectory(prefix="ctrlpp-doc-snippets-") as directory, \
            ThreadPoolExecutor(max_workers=args.jobs) as pool:
        return list(pool.map(lambda item: compile_snippet(args, Path(directory), *item), enumerate(snippets)))


def requested(root: Path, paths: list[str], pages: list[str], parser: argparse.ArgumentParser) -> list[str]:
    selected = set()
    for given in paths:
        path = Path(given)
        try:
            name = path.resolve().relative_to(root.resolve()).as_posix() if path.is_absolute() else \
                posixpath.normpath(path.as_posix())
        except ValueError:
            name = given
        under = [page for page in pages if name == "." or page == name or page.startswith(name.rstrip("/") + "/")]
        if not under:
            parser.error(f"{given} names no tracked documentation page (docs/**/*.md or README.md)")
        selected.update(under)
    return sorted(selected) if paths else pages


def report(snippets: list[Snippet], results: list[tuple[int | None, str]]) -> int:
    failed = [(snippet, output) for snippet, (code, output) in zip(snippets, results) if code]
    for snippet, output in failed:
        print(f"error: {snippet.block.page}:{snippet.block.fence}: does not compile\n{output}", end="",
              file=sys.stderr)
    compiled = [snippet for snippet in snippets if snippet.kind != "fragment"]
    pages = len({snippet.block.page for snippet in snippets})
    print(f"checked {len(compiled)} snippets in {pages} pages: {len(failed)} failed, "
          f"{len(snippets) - len(compiled)} fragments skipped, "
          f"{sum(snippet.kind == 'exceptions' for snippet in snippets)} compiled with exceptions")
    return EXIT_FAILED if failed else 0


def run(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    try:
        pages = requested(args.source_root, args.paths, tracked_pages(args.source_root), parser)
        snippets, errors = scan(args.source_root, pages)
    except (OSError, subprocess.CalledProcessError) as exc:
        print(f"error: cannot read the tracked pages under {args.source_root}: {exc}", file=sys.stderr)
        return EXIT_NO_INPUT
    for page, line, reason in sorted(errors):
        print(f"error: {page}:{line}: {reason}", file=sys.stderr)
    if not snippets:
        print(f"error: no tracked page under {args.source_root} holds a snippet", file=sys.stderr)
        return EXIT_MARKER if errors else EXIT_NO_INPUT
    if args.list:
        print("".join(f"{snippet.block.page}:{snippet.block.fence} {snippet.kind}\n" for snippet in snippets), end="")
        return EXIT_MARKER if errors else 0
    compiled = iter(compile_all(args, [snippet for snippet in snippets if snippet.kind != "fragment"]))
    results = [(0, "") if snippet.kind == "fragment" else next(compiled) for snippet in snippets]
    unrunnable = next((output for code, output in results if code is None), None)
    if unrunnable is not None:
        print(f"error: cannot execute the compiler {args.compiler}: {unrunnable}", file=sys.stderr)
        return EXIT_NO_COMPILER
    status = report(snippets, results)
    return EXIT_MARKER if errors else status


def jobs(text: str) -> int:
    count = int(text)
    if count < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return count


def parser_for() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, epilog=EXIT_HELP,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    add = parser.add_argument
    add("--source-root", type=Path, default=Path("."), help="repository root (default: the current directory)")
    add("--compiler", default=os.environ.get("CXX") or "c++", help="C++ compiler (default: $CXX, else c++)")
    add("--eigen-include", type=Path, default=Path("/usr/include/eigen3"), help="Eigen include directory")
    add("-j", "--jobs", type=jobs, default=os.cpu_count() or 1,
        help="compilers run at once (default: the number of processors)")
    add("--list", action="store_true", help="print each snippet as PAGE:FENCE-LINE KIND and compile nothing")
    add("paths", nargs="*", metavar="PATH",
        help="pages or directories to check, relative to the root (default: every tracked page)")
    return parser


if __name__ == "__main__":
    PARSER = parser_for()
    ARGS = PARSER.parse_args()
    ARGS.extra_args = shlex.split(os.environ.get(EXTRA_ARG_VARIABLE, ""))
    raise SystemExit(run(ARGS, PARSER))
