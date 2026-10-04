#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_TIMING_RUN_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_TIMING_RUN_H

#include "cycle_counter.h"

#include "golden_reference.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace ctrlpp {

// Enables and proves the counter, then times the instrument's own overhead, the
// gain design and a golden-length block of each family's steps. Each region's
// first execution after reset is its cold figure, taken before anything else
// runs it, so this must run ahead of the families' golden runs.
counter_state run_timing();

// Brackets each solve of the predictive run, whose controller cannot be built a
// second time beside the one the golden run judges.
class solve_timer
{
public:
    solve_timer();

    void enter(std::size_t k) noexcept;

    void leave(std::size_t k, std::int32_t iterations) noexcept;

    // Prints nothing unless the counter still proves live.
    void report() const;

private:
    std::array<cycle_span, kPredictiveRunSteps> spans_;
    std::array<std::int32_t, kPredictiveRunSteps> iterations_;
    std::uint32_t start_;

    void report_work(std::uint32_t slowest_at) const;
};

}

#endif
