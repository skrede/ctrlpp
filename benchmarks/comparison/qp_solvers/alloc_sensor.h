#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ALLOC_SENSOR_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ALLOC_SENSOR_H

#include <cstdint>

namespace ctrlpp::bench::qp_solvers
{

/// Counts the C allocator's allocation calls (malloc, calloc, realloc, posix_memalign, aligned_alloc,
/// memalign) while armed, from every library in the process, after argmin's
/// tests/unit/alloc_sensor_b_linux.cpp. The sensor exists only where the build compiled its glibc
/// translation unit; elsewhere the count is reported unavailable.
auto alloc_sensor_present() -> bool;

void arm_alloc_sensor();

/// Disarms and returns the count since arming, or -1 without the sensor.
auto disarm_alloc_sensor() -> std::int64_t;

}

#endif
