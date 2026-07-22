#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H

namespace ctrlpp
{

constexpr double kDt    = 0.02;
constexpr int    kSteps = 200;

// Host-double reference for the shared double-integrator problem, produced by an
// independent host build (generate_golden.cpp) so on-target results are diffed
// against an external reference rather than self-reported.
constexpr double kHostK0        = 9.467093672;
constexpr double kHostK1        = 5.281739638;
constexpr double kHostFinalNorm = 2.020686726e-05;

}

#endif
