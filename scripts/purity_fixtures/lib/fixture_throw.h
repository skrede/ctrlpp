#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_THROW_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_THROW_H

/// Deliberately non-conforming. Never compiled, part of no build target.
///
/// This comment contains a bare throw statement so the comment filter is
/// exercised. Remove that filter and this file reports three hits instead of
/// two, which the self-test's exact-count assertion catches. Twelve comments
/// in the real library match the same pattern, which is why the filter exists.
///
/// The two rethrows below prove the sanctioned construct in both directions:
/// the one inside the forced-unwind handler must not count, and the one inside
/// the catch-all that follows it must.

#include <cxxabi.h>

struct fixture_error
{
};

inline void fixture_raise()
{
    throw fixture_error{};
}

inline void fixture_rethrow()
{
    try
    {
        fixture_raise();
    }
    catch(const abi::__forced_unwind&)
    {
        throw;
    }
    catch(...)
    {
        throw;
    }
}

#endif
