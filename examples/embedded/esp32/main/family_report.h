#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_FAMILY_REPORT_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_FAMILY_REPORT_H

#include "alloc_sensor.h"

#include "golden_verdict.h"
#include "golden_reference.h"
#include "control_loop_demo.h"

#include <cmath>
#include <cstddef>
#include <cstdint>

namespace ctrlpp {

struct window_figures
{
    bool measured;
    std::uint32_t allocations;
    std::uint32_t eigen_trips;
    heap_view heap;
    std::uint32_t heap_self;
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
// line to a reader rather than an omission it cannot see.
void report_family(const family_result &result, const window_figures &figures);

void report_refusal(const char *name, const char *reason);

void report_absent_predictive();

void report_golden(const control_loop_demo<float> &demo, const golden_verdict &verdict, const window_figures &figures);

// Armed around the steps alone and reset per family, so a nonzero count belongs
// to this family and not to the program. The value is widened to double before
// the subtraction, which then rounds once in double, far below the
// single-precision bound it is judged against.
template<template<class> class Demo>
void drive_family(const char *name, double golden, double bound)
{
    auto demo = Demo<float>::make();
    if(!demo.has_value())
    {
        report_refusal(name, describe(demo.error()));
        return;
    }
    alloc_sensor_reset();
    alloc_sensor_arm();
    float value = 0.0f;
    for(std::size_t k = 0; k < kGoldenSteps; ++k)
        value = demo->step();
    alloc_sensor_disarm();
    const double departure = std::abs(static_cast<double>(value) - golden);
    report_family({name, static_cast<double>(value), golden, departure, bound, within_bound(departure, bound)}, read_window());
}

}

#endif
