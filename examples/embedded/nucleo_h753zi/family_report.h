#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_FAMILY_REPORT_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_FAMILY_REPORT_H

#include "alloc_sensor.h"

#include "golden_verdict.h"
#include "golden_reference.h"

#include <cmath>
#include <cstdio>
#include <cstddef>
#include <cstdint>

namespace ctrlpp {

struct window_figures
{
    std::uint32_t allocations;
    std::uint32_t eigen_trips;
};

struct family_result
{
    const char *name;
    double value;
    double golden;
    double departure;
    double bound;
    bool pass;
};

window_figures read_window() noexcept;

void report_allocations(const window_figures &figures);

// One line per family, so a family that ran and printed nothing is a missing
// tag to the capture script rather than an omission it cannot see.
void report_family(const family_result &result, const window_figures &figures);

// Armed around the steps alone and reset per family, so a nonzero count belongs
// to this family and not to the program.
template<template<class> class Demo>
void drive_family(const char *name, double golden, double bound)
{
    auto demo = Demo<double>::make();
    if(!demo.has_value())
    {
        std::printf("[family] %s verdict=REFUSED -- %s\n", name, describe(demo.error()));
        return;
    }
    alloc_sensor_reset();
    alloc_sensor_arm();
    double value = 0.0;
    for(std::size_t k = 0; k < kGoldenSteps; ++k)
        value = demo->step();
    alloc_sensor_disarm();
    const double departure = std::abs(value - golden);
    report_family({name, value, golden, departure, bound, within_bound(departure, bound)}, read_window());
}

}

#endif
