#ifndef HPP_GUARD_BENCHMARKS_BENCH_CONSTRUCT_H
#define HPP_GUARD_BENCHMARKS_BENCH_CONSTRUCT_H

/// @brief Unwrap a fallible ctrlpp factory result inside a benchmark harness.
///
/// Every ctrlpp type is built through a factory returning ctrlpp::expected.
/// Benchmarks build from fixed, known-good configurations, so a rejection means
/// the harness itself is misconfigured; it is reported and the run stops. A
/// benchmark must never measure a stand-in object it silently substituted, and
/// it must never report a number for a configuration the library refused.

#include <cstdio>
#include <cstdlib>
#include <utility>
#include <type_traits>

namespace ctrlpp::bench
{

template <typename Built>
auto built_or_exit(Built&& built, char const* what) -> std::remove_cvref_t<decltype(*built)>
{
    if(!built.has_value())
    {
        std::fprintf(stderr, "benchmark configuration rejected: %s\n", what);
        std::exit(EXIT_FAILURE);
    }
    return *std::forward<Built>(built);
}

}

#endif
