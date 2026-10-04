#include "cycle_record.h"

#include <cstdio>
#include <cstdint>
#include <cinttypes>

namespace ctrlpp {

namespace {

constexpr std::uint64_t kTenthsPerCycle = 10;

}

void report_timing(const timing_region &region, const cycle_distribution &d)
{
    if(d.n == 0)
    {
        std::printf("[timing] region=%s regime=%s verdict=WRAPPED wraps=%" PRIu32 "\n", region.name, region.regime, d.wraps);
        return;
    }
    const std::uint64_t tenths = (std::uint64_t{d.max} * kTenthsPerCycle + region.block / 2U) / region.block;
    std::printf("[timing] region=%s regime=%s n=%" PRIu32 " n_basis=%s block=%" PRIu32 " raw_min=%" PRIu32 " raw_median=%" PRIu32 " raw_max=%" PRIu32, region.name, region.regime, d.n,
                region.basis, region.block, d.min, d.median, d.max);
    std::printf(" max_at=%" PRIu32 " per_iteration_raw_max=%" PRIu32 ".%" PRIu32 " wraps=%" PRIu32 "\n", d.max_at, static_cast<std::uint32_t>(tenths / kTenthsPerCycle),
                static_cast<std::uint32_t>(tenths % kTenthsPerCycle), d.wraps);
}

}
