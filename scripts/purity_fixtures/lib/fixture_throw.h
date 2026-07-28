#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_THROW_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_THROW_H

/// Deliberately non-conforming. Never compiled, part of no build target.
///
/// This comment contains a bare throw statement so the comment filter is
/// exercised. Remove that filter and this file reports two hits instead of
/// one, which the self-test's exact-count assertion catches. Twelve comments
/// in the real library match the same pattern, which is why the filter exists.

struct fixture_error
{
};

inline void fixture_raise()
{
    throw fixture_error{};
}

#endif
