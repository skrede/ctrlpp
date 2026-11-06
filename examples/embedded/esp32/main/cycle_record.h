#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_CYCLE_RECORD_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_ESP32_MAIN_CYCLE_RECORD_H

#include "cycle_counter.h"

#include "golden_reference.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <algorithm>

namespace ctrlpp
{

// The largest sample set any region records: one warm repeat per golden step.
inline constexpr std::size_t timing_capacity = golden::sample_count;

struct cycle_distribution
{
    std::uint32_t n;
    std::uint32_t wraps;
    std::uint32_t min;
    std::uint32_t median;
    std::uint32_t max;
    std::uint32_t max_at;
    std::uint32_t tick_crossed;
};

struct timing_region
{
    const char* name;
    const char* regime;
    const char* basis;
    std::uint32_t block;
};

// A wrapped sample is counted and kept out of the figures. The worst sample's
// position is kept, because a worst that is always the first repeat says the
// region was re-warmed by the code that ran before it. A tick crossing is
// counted over every sample, wrapped or not.
class cycle_record
{
public:
    cycle_record()
        : m_cycles{}
        , m_n(0)
        , m_wraps(0)
        , m_max_at(0)
        , m_tick_crossed(0)
    {
    }

    void clear() noexcept
    {
        m_n = 0;
        m_wraps = 0;
        m_max_at = 0;
        m_tick_crossed = 0;
    }

    void add(cycle_span span) noexcept
    {
        m_tick_crossed += span.tick_crossed ? 1U : 0U;
        if(span.wrapped)
        {
            ++m_wraps;
            return;
        }
        if(m_n > 0 && span.cycles > m_cycles[m_max_at])
            m_max_at = m_n;
        m_cycles[m_n++] = span.cycles;
    }

    cycle_distribution summarize() noexcept
    {
        if(m_n == 0)
            return {0, m_wraps, 0, 0, 0, 0, m_tick_crossed};
        const std::uint32_t max_at = m_max_at;
        std::sort(m_cycles.begin(), m_cycles.begin() + m_n);
        return {m_n, m_wraps, m_cycles[0], m_cycles[m_n / 2], m_cycles[m_n - 1], max_at, m_tick_crossed};
    }

private:
    std::array<std::uint32_t, timing_capacity> m_cycles;
    std::uint32_t m_n;
    std::uint32_t m_wraps;
    std::uint32_t m_max_at;
    std::uint32_t m_tick_crossed;
};

}

#endif
