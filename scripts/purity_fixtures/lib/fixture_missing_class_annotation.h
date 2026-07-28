#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_MISSING_CLASS_ANNOTATION_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_MISSING_CLASS_ANNOTATION_H

/// Deliberately non-conforming. Never compiled, part of no build target.
///
/// The self-test's manifest carries two rows for this file. The first names
/// `fixture_present_disposition`, whose attribute IS below, so rule 4a's
/// positive branch and rule 4b's manifest exemption are both exercised. The
/// second names `struct [[nodiscard]] fixture_disposition`, whose attribute is
/// deliberately absent, and that is the one violation this fixture produces.
///
/// The pairing is what lets the self-test tell a deleted fixture apart from a
/// deleted annotation: deleting this file breaks both rows and the count goes
/// to two instead of one.
///
/// The second row's exact text appears in THIS comment on purpose: a rule
/// written as "the attribute appears somewhere in the file" would pass against
/// this fixture, and rule 4a must not. That is what makes 4a a presence
/// assertion on an exact declaration rather than a text search, and the
/// commented occurrence is also what proves rule 4b's comment filter is live
/// here.

struct [[nodiscard]] fixture_present_disposition
{
    double commanded{};
};

struct fixture_disposition
{
    double commanded{};
    double realized{};
};

#endif
