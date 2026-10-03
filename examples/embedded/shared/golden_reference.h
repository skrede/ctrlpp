#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_GOLDEN_REFERENCE_H

#include <cstddef>

namespace ctrlpp
{

constexpr double kDt    = 0.02;
constexpr int    kSteps = 200;

// The reference below is an independent host BUILD of the same algorithm at
// double precision (generate_golden.cpp), not an independent implementation. It
// establishes that a board's arithmetic and its FPU agree with a desktop double
// build on this problem, which is a portability claim; it cannot detect a wrong
// Riccati solve, because both sides would be wrong identically. The
// independent-implementation check is the cross-validation suite under
// validation/. The constants carry 17 significant digits so each is the host's
// double exactly, which is what lets the verdict charge the host run nothing
// for its gain.
constexpr double kHostK0        = 9.4670936721545083;
constexpr double kHostK1        = 5.281739638040154;
constexpr double kHostFinalNorm = 2.0206867258019932e-05;

// The kernel runs one step per index from 0 through kSteps inclusive.
constexpr std::size_t kGoldenSteps = static_cast<std::size_t>(kSteps) + 1;
constexpr std::size_t kGoldenNx    = 2;
constexpr std::size_t kGoldenNu    = 1;

}

#endif
