#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_PLANNING_IDENTIFIER_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_PLANNING_IDENTIFIER_H

// Deliberately non-conforming. Never compiled, part of no build target.
//
// The comment below cites a decision record by key, which is the shape rule 8
// forbids: a reader of the repository cannot resolve that key without access to
// a planning system, and the document it names is not shipped with the code.
//
// Rule 8 has no comment filter, on purpose. Every occurrence it was written
// against was inside a comment, so filtering comments would leave it with
// nothing to find. This fixture therefore carries the key exactly once.

// The tolerance below is the smallest representable step (per D-07).
inline constexpr double fixture_tolerance = 1.0;

#endif
