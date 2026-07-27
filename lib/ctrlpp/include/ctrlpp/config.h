#ifndef HPP_GUARD_CTRLPP_CONFIG_H
#define HPP_GUARD_CTRLPP_CONFIG_H

/// @brief Library-wide build configuration macros. This header is dependency-free
/// by design (zero includes) so embedded builds can consume it before anything else.
///
/// CTRLPP_NO_EXCEPTIONS: define at build time to opt out of the one place the
/// library is contractually obliged to throw. No construction path is affected:
/// every type is built through a fallible factory returning ctrlpp::expected,
/// and those are available in every build; so is the fallible `setup(problem)`
/// on the optional backend adapters, which is their only setup shape. What the
/// macro gates is the throw that ctrlpp::expected::value() is contractually
/// required to have, which becomes std::abort() instead.
///
/// CTRLPP_HAS_EXCEPTIONS: expands to 1 exactly when the toolchain has exception
/// support enabled and CTRLPP_NO_EXCEPTIONS is not defined; expands to 0
/// otherwise. Two toolchain spellings are accepted: the standard feature-test
/// macro __cpp_exceptions, and _CPPUNWIND, which is what the Microsoft toolchain
/// has always set from its /EH switch. Accepting both keeps a supported compiler
/// from silently landing on the exception-free contract while its runtime still
/// unwinds; neither macro is defined by a compiler whose exceptions are off, so
/// the disjunction cannot turn them on anywhere.

#if (defined(__cpp_exceptions) || defined(_CPPUNWIND)) && !defined(CTRLPP_NO_EXCEPTIONS)
    #define CTRLPP_HAS_EXCEPTIONS 1
#else
    #define CTRLPP_HAS_EXCEPTIONS 0
#endif

#endif
