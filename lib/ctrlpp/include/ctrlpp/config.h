#ifndef HPP_GUARD_CTRLPP_CONFIG_H
#define HPP_GUARD_CTRLPP_CONFIG_H

/// @brief Library-wide build configuration macros. This header is dependency-free
/// by design (zero includes) so embedded builds can consume it before anything else.
///
/// CTRLPP_NO_EXCEPTIONS: define at build time to opt out of every throwing
/// convenience wrapper in the library. The expected-based factories and solvers
/// remain available unconditionally.
///
/// CTRLPP_HAS_EXCEPTIONS: expands to 1 exactly when the compiler has exception
/// support enabled (__cpp_exceptions) and CTRLPP_NO_EXCEPTIONS is not defined;
/// expands to 0 otherwise.

#if defined(__cpp_exceptions) && !defined(CTRLPP_NO_EXCEPTIONS)
    #define CTRLPP_HAS_EXCEPTIONS 1
#else
    #define CTRLPP_HAS_EXCEPTIONS 0
#endif

#endif
