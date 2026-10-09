#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_FAMILY_REPORT_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_FAMILY_REPORT_H

#include "golden_verdict.h"
#include "golden_reference.h"
#include "control_loop_demo.h"

#include <cmath>
#include <cstddef>

namespace ctrlpp {

struct family_result
{
    const char *name;
    double value;
    double golden;
    double departure;
    double bound;
    bool pass;
};

// One line per family, so a family that ran and printed nothing is a missing
// line to a reader rather than an omission it cannot see.
void report_family(const family_result &result);

void report_refusal(const char *name, const char *reason);

void report_unmeasured_allocations();

void report_absent_predictive();

void report_golden(const control_loop_demo<float> &demo, const golden_verdict &verdict);

// The value is widened to double before the subtraction, which then rounds once
// in double, far below the single-precision bound it is judged against.
template<template<class> class Demo>
void drive_family(const char *name, double golden, double bound)
{
    auto demo = Demo<float>::make();
    if(!demo.has_value())
    {
        report_refusal(name, describe(demo.error()));
        return;
    }
    float value = 0.0f;
    for(std::size_t k = 0; k < kGoldenSteps; ++k)
        value = demo->step();
    const double departure = std::abs(static_cast<double>(value) - golden);
    report_family({name, static_cast<double>(value), golden, departure, bound, within_bound(departure, bound)});
}

}

#endif
