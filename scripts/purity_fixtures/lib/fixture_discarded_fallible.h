#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_DISCARDED_FALLIBLE_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_DISCARDED_FALLIBLE_H

/// Deliberately non-conforming. Never compiled, part of no build target.
///
/// A bare-statement call to create(...), which the gate's callable list names
/// as fallible on every declaration in the library. The typed failure it
/// returns is dropped on the floor.
///
/// Rule 4c anchors its pattern at the start of the statement, so a mention
/// inside a comment can never match it and the comment above is documentation
/// rather than a filter exercise. The rules whose filters the fixtures do
/// exercise are 1, 2, 3 and 4b.

struct fixture_config
{
};

struct fixture_thing
{
    static auto create(const fixture_config& cfg) -> int;
};

inline void fixture_drop()
{
    fixture_thing::create(fixture_config{});
}

#endif
