#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_SITE_ANNOTATION_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_SITE_ANNOTATION_H

/// Deliberately non-conforming. Never compiled, part of no build target.
///
/// A site-level [[nodiscard]] on a trivial size query. The getter shape is
/// deliberate: it is exactly what the scope fence exists to reject, and a
/// fixture using a fallible return would blur the distinction rule 4b draws
/// between a discard that silently breaks correctness and a discard that is
/// merely pointless.
///
/// The occurrence in this comment exercises the comment filter; remove the
/// filter and this file reports two hits instead of one.

struct fixture_widget
{
    [[nodiscard]] int size() const { return count_; }

    int count_{0};
};

#endif
