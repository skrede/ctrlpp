#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_UNWRAP_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_UNWRAP_H

/// Deliberately non-conforming. Never compiled, part of no build target.
///
/// This comment repeats the banned pattern -- a call to .value() -- so the
/// comment filter is exercised. Remove that filter and this file reports two
/// hits instead of one, which the self-test's exact-count assertion catches.

struct fixture_optional
{
    int value() const;
};

inline int fixture_unwrap(const fixture_optional& r)
{
    return r.value();
}

#endif
