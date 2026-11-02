"""Style rules over the data members of every record: member naming and in-declaration member initializers."""

import re

from style_source import Source
from style_scopes import Member, Record, data_members, records
from style_config import CALLBACK_TYPE_PATTERNS, CARVE_OUT_NESTED, CARVE_OUT_SUFFIXES

Finding = tuple[int, str]

LEADING_SPECIFIERS = r"(?:(?:static|inline|constexpr|mutable|const)\s+)*"
CALLBACK_TYPES = tuple(re.compile(LEADING_SPECIFIERS + pattern) for pattern in CALLBACK_TYPE_PATTERNS)


def is_callback(member: Member) -> bool:
    return any(pattern.match(member.type) for pattern in CALLBACK_TYPES)


def naming_problems(member: Member) -> list[str]:
    problems = []
    if member.static and member.name.startswith("m_"):
        problems.append("a static member takes no m_ prefix")
    if not member.static and member.access != "public" and not member.name.startswith("m_"):
        problems.append(f"a {member.access} member takes the m_ prefix")
    if member.name.endswith("_"):
        problems.append("no member name ends in an underscore")
    if is_callback(member) and not member.name.endswith("_cb"):
        problems.append("a callback member takes the _cb suffix")
    return problems


def check_member_naming(source: Source) -> list[Finding]:
    return [(member.line, f"{record.qualified}::{member.name}: {'; '.join(problems)}")
            for record in records(source.blanked) for member in data_members(record)
            if (problems := naming_problems(member))]


def is_carve_out(record: Record, members: list[Member]) -> bool:
    named = record.name is not None and (record.name.endswith(CARVE_OUT_SUFFIXES) or
                                         (record.nested and record.name == CARVE_OUT_NESTED))
    return named and record.kind == "struct" and not record.has_constructor and \
        all(member.access == "public" for member in members if not member.static)


def initialized(source: Source) -> list[tuple[Record, Member]]:
    found = []
    for record in records(source.blanked):
        members = data_members(record)
        if not is_carve_out(record, members):
            found += [(record, member) for member in members if member.init is not None and not member.static]
    return found


def check_member_init(source: Source) -> list[Finding]:
    return [(member.line, f"{record.qualified}::{member.name} is initialized in its declaration; initialize it in "
                          "a constructor's init list") for record, member in initialized(source)]


def list_members(source: Source) -> list[str]:
    return [f"{source.path}:{member.line}: {record.qualified}: {member.name} "
            f"{' '.join(source.raw[member.init[0]:member.init[1]].split())}" for record, member in initialized(source)]
