"""The include-order rule: sections, folder groups, the file's own folder first, and order within a group."""

import re
import posixpath
from fnmatch import fnmatchcase
from typing import NamedTuple

from style_source import Source
from style_config import INCLUDE_ROOTS, ORDER_DEPENDENCIES, STANDARD_HEADERS

INCLUDE = re.compile(r'^[ \t]*#[ \t]*include[ \t]*([<"])([^>"]*)[>"]')
SECTIONS = ("project", "third-party", "standard")
PROJECT, THIRD_PARTY, STANDARD = range(3)
LIBRARY_ROOT = INCLUDE_ROOTS[0] + "/"


class Include(NamedTuple):
    line: int
    shown: str
    spelled: str
    section: int
    folder: str


class Finding(NamedTuple):
    anchor: Include
    offender: Include
    message: str
    ordering: bool = True


def resolve(source: Source, spelled: str) -> str | None:
    for root in (posixpath.dirname(source.path), *INCLUDE_ROOTS):
        candidate = posixpath.normpath(posixpath.join(root, spelled))
        if candidate in source.tracked:
            return candidate
    return None


def classify(source: Source, line: int, delimiter: str, spelled: str) -> Include:
    shown = f'"{spelled}"' if delimiter == '"' else f"<{spelled}>"
    resolved = resolve(source, spelled) if delimiter == '"' else None
    if resolved is not None:
        return Include(line, shown, spelled, PROJECT, posixpath.dirname(resolved))
    if delimiter == "<" and spelled in STANDARD_HEADERS:
        return Include(line, shown, spelled, STANDARD, "")
    return Include(line, shown, spelled, THIRD_PARTY, posixpath.dirname(spelled))


def folder_name(section: int, folder: str) -> str:
    if section == STANDARD:
        return "the standard library"
    if section == PROJECT:
        return folder.removeprefix(LIBRARY_ROOT)
    return f"third-party {folder or 'top-level'}"


def key(include: Include) -> tuple[int, str]:
    return include.section, include.folder


def name_order(include: Include) -> tuple[int, str]:
    name = posixpath.basename(include.spelled)
    return len(name), name


def block_tokens(source: Source) -> list[tuple[str, int, Include | None]]:
    """The include block, every line up to the first line of code, as (kind, line, include) tokens."""
    tokens, continued = [], False
    for number, (raw, blanked) in enumerate(zip(source.raw.split("\n"), source.blanked.split("\n")), 1):
        code = blanked.strip()
        if continued or (code.startswith("#") and not INCLUDE.match(raw)):
            tokens.append(("break", number, None))
        elif code.startswith("#"):
            match = INCLUDE.match(raw)
            tokens.append(("include", number, classify(source, number, match.group(1), match.group(2))))
        elif code:
            break
        else:
            tokens.append(("comment" if raw.strip() else "blank", number, None))
        continued = code.endswith("\\")
    return tokens


def regions(tokens: list) -> list[list[list[Include]]]:
    """Split the block at every non-include directive, and each region into groups at blank lines.

    A comment line neither ends a group nor starts one, so a comment above an include stays with its group.
    """
    found, groups, gap = [], [], True
    for kind, _, include in tokens + [("break", 0, None)]:
        if kind == "include":
            if gap:
                groups.append([])
            groups[-1].append(include)
        if kind == "break" and groups:
            found.append(groups)
            groups = []
        if kind != "comment":
            gap = kind != "include"
    return found


def blank_runs(tokens: list) -> list[tuple[int, str]]:
    places = [index for index, (kind, _, _) in enumerate(tokens) if kind == "include"]
    findings, run = [], 0
    for index, (kind, line, _) in enumerate(tokens):
        if kind == "blank":
            run += 1
            continue
        if run > 1 and places and places[0] < index <= places[-1]:
            findings.append((line, f"{run} blank lines before this line; include groups are separated by exactly one"))
        run = 0
    return findings


def group_findings(group: list[Include]) -> list[Finding]:
    findings = []
    for previous, include in zip(group, group[1:]):
        if key(previous) != key(include):
            findings.append(Finding(previous, include, f"{include.shown} ({folder_name(*key(include))}) shares a group "
                                    f"with {previous.shown} ({folder_name(*key(previous))}); give each folder its "
                                    "own group", ordering=False))
        elif name_order(include) < name_order(previous):
            findings.append(Finding(previous, include, f"{include.shown} sorts before {previous.shown}: shorter file "
                                    "names first, equal lengths alphabetically"))
    return findings


def section_findings(groups: list[list[Include]]) -> list[Finding]:
    findings, latest = [], groups[0][0]
    for group in groups[1:]:
        if group[0].section < latest.section:
            findings.append(Finding(latest, group[0], f"{group[0].shown} is a {SECTIONS[group[0].section]} include "
                                    f"after the {SECTIONS[latest.section]} group of {latest.shown}; order the "
                                    "sections project, third-party, standard"))
        elif group[0].section > latest.section:
            latest = group[0]
    return findings


def folder_findings(source: Source, groups: list[list[Include]]) -> list[Finding]:
    findings, seen, foreign = [], {}, None
    own = (PROJECT, posixpath.dirname(source.path))
    for group in groups:
        lead = group[0]
        if key(lead) in seen:
            findings.append(Finding(seen[key(lead)], lead, f"{folder_name(*key(lead))} is split into more than one "
                                    f"group at {lead.shown}; keep a folder's includes together"))
        elif key(lead) == own and foreign is not None:
            findings.append(Finding(foreign, lead, f"the file's own folder {folder_name(*own)} must lead the project "
                                    f"includes, ahead of {folder_name(*key(foreign))}"))
        elif lead.section == PROJECT and foreign is None and key(lead) != own:
            foreign = lead
        seen.setdefault(key(lead), lead)
    return findings


def exempt(path: str, finding: Finding) -> bool:
    return finding.ordering and any(fnmatchcase(path, pattern) and finding.anchor.spelled == earlier
               and fnmatchcase(finding.offender.spelled, later) for pattern, earlier, later, _ in ORDER_DEPENDENCIES)


def check_include_order(source: Source) -> list[tuple[int, str]]:
    tokens = block_tokens(source)
    findings = []
    for groups in regions(tokens):
        findings += [finding for group in groups for finding in group_findings(group)]
        findings += section_findings(groups) + folder_findings(source, groups)
    kept = {(finding.offender.line, finding.message) for finding in findings if not exempt(source.path, finding)}
    return sorted(kept | set(blank_runs(tokens)))
