#ifndef HPP_GUARD_CTRLPP_CONFIG_H
#define HPP_GUARD_CTRLPP_CONFIG_H

// Includes nothing, so an embedded build can include it before anything else. Defining CTRLPP_NO_EXCEPTIONS
// turns the one throw the library is required to have, in ctrlpp::expected::value(), into std::abort().
// _CPPUNWIND is the Microsoft toolchain's counterpart of __cpp_exceptions, set from its /EH switch.
#if(defined(__cpp_exceptions) || defined(_CPPUNWIND)) && !defined(CTRLPP_NO_EXCEPTIONS)
    #define CTRLPP_HAS_EXCEPTIONS 1
#else
    #define CTRLPP_HAS_EXCEPTIONS 0
#endif

#endif
