#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_TIMING_RUN_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_TIMING_RUN_H

#include "cycle_counter.h"

#include "golden_reference.h"

#include <array>
#include <cstddef>
#include <cstdint>

namespace ctrlpp
{

// Each region's first execution after reset is its cold figure, so this runs ahead of the families'
// golden runs.
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
    std::array<cycle_span, kPredictiveRunSteps> m_spans;
    std::array<std::int32_t, kPredictiveRunSteps> m_iterations;
    std::uint32_t m_start;

    void report_work(std::uint32_t slowest_at) const;
};

}

#endif
