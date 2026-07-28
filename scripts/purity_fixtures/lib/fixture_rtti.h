#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_RTTI_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_RTTI_H

/// Deliberately non-conforming. Never compiled, part of no build target.
///
/// This comment mentions dynamic_cast<T> so the comment filter is exercised.
/// Remove that filter and this file reports two hits instead of one.

struct fixture_base
{
    virtual ~fixture_base();
};

struct fixture_derived : fixture_base
{
};

inline const fixture_derived* fixture_downcast(const fixture_base* b)
{
    return dynamic_cast<const fixture_derived*>(b);
}

#endif
