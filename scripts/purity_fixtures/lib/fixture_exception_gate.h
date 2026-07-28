#ifndef HPP_GUARD_PURITY_FIXTURES_FIXTURE_EXCEPTION_GATE_H
#define HPP_GUARD_PURITY_FIXTURES_FIXTURE_EXCEPTION_GATE_H

/// Deliberately non-conforming. Never compiled, part of no build target.
///
/// The exception-mode macro appears in a file that is neither the
/// fallible-result backport header nor the configuration header that defines
/// it. Rule 5 is a file-presence rule with no comment filter, on purpose: a
/// mention of that macro anywhere else in the library is itself the signal, so
/// this comment and the directive below are one violation rather than two.

#if CTRLPP_HAS_EXCEPTIONS
inline constexpr bool fixture_exceptions_enabled = true;
#else
inline constexpr bool fixture_exceptions_enabled = false;
#endif

#endif
