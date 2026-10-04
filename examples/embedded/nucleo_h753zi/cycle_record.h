#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_CYCLE_RECORD_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_CYCLE_RECORD_H

#include "cycle_counter.h"

#include "golden_reference.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ctrlpp {

// The largest sample set any region records: every solve of the predictive run.
constexpr std::size_t kTimingCapacity = kPredictiveRunSteps;

struct cycle_distribution
{
    std::uint32_t n;
    std::uint32_t wraps;
    std::uint32_t min;
    std::uint32_t median;
    std::uint32_t max;
    std::uint32_t max_at;
};

struct timing_region
{
    const char *name;
    const char *regime;
    const char *basis;
    std::uint32_t block;
};

// A wrapped sample is counted and kept out of the figures. The worst sample's
// position is kept, because a worst that is always the first repeat says the
// region was re-warmed by the code that ran before it.
class cycle_record
{
public:
    cycle_record()
            : cycles_{}
            , n_(0)
            , wraps_(0)
            , max_at_(0)
    {
    }

    void add(cycle_span span) noexcept
    {
        if(span.wrapped)
        {
            ++wraps_;
            return;
        }
        if(n_ > 0 && span.cycles > cycles_[max_at_])
            max_at_ = n_;
        cycles_[n_++] = span.cycles;
    }

    cycle_distribution summarize() noexcept
    {
        if(n_ == 0)
            return {0, wraps_, 0, 0, 0, 0};
        const std::uint32_t max_at = max_at_;
        std::sort(cycles_.begin(), cycles_.begin() + n_);
        return {n_, wraps_, cycles_[0], cycles_[n_ / 2], cycles_[n_ - 1], max_at};
    }

private:
    std::array<std::uint32_t, kTimingCapacity> cycles_;
    std::uint32_t n_;
    std::uint32_t wraps_;
    std::uint32_t max_at_;
};

// One tagged line per region and regime. The figures are raw: each carries the
// instrument overhead, which the report publishes on its own lines.
void report_timing(const timing_region &region, const cycle_distribution &d);

}

#endif
