#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_MILESTONE_IDENTIFIER_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_MILESTONE_IDENTIFIER_H

// Deliberately non-conforming. Never compiled, part of no build target.
//
// The branch reference below uses the versioned development-branch shape that
// the shipped-material gate rejects. It occurs exactly once so the self-test
// can distinguish this family from the decision-record fixture.

// This behavior first shipped on milestone/v7.8.9.
inline constexpr bool fixture_milestone_behavior = true;

#endif
