#include "stack_watermark_board.h"

#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cinttypes>

extern "C" std::uint32_t _ebss[];
extern "C" std::uint32_t _estack[];

namespace ctrlpp {

namespace {

// The harness's own reach below the origin: the walk's caller pushes three
// words below it. A build that deepens the harness reads a nonzero floor, which
// withholds the figures rather than misstating them.
#ifndef CTRLPP_MCU_STACK_GAP_BYTES
#define CTRLPP_MCU_STACK_GAP_BYTES 12
#endif

constexpr std::uint32_t kGapBytes = CTRLPP_MCU_STACK_GAP_BYTES;

static_assert(kGapBytes % sizeof(std::uint32_t) == 0, "the window's top must sit on the word grid");

// The host instrument's seed cut to this part's 32-bit word. The address
// exclusive-or works at the pointer width on either part, so no single stored
// value can pass as undisturbed at more than one address.
constexpr std::uint32_t kPaintSeed = static_cast<std::uint32_t>(0x5FC3A96E1D47B208ULL);

struct program_record
{
    std::uint32_t *top;
    std::uint32_t *deepest;
    std::uint32_t floor;
};

constinit program_record program{nullptr, nullptr, 0};

// The heap lives in AXI SRAM, so nothing in DTCM sits between the end of bss
// and the stack. The stack descends through that memory long before it reaches
// the nominal reserve, so the reserve is not the window's honest lower bound and
// the end of bss is.
std::uint32_t *window_bottom() noexcept
{
    return _ebss;
}

std::uintptr_t stack_top() noexcept
{
    return reinterpret_cast<std::uintptr_t>(_estack);
}

std::uint32_t paint_word(const std::uint32_t *address) noexcept
{
    return kPaintSeed ^ static_cast<std::uint32_t>(reinterpret_cast<std::uintptr_t>(address));
}

[[gnu::noinline]] void paint_below(std::uint32_t *top) noexcept
{
    for(std::uint32_t *word = window_bottom(); word != top; ++word)
        *word = paint_word(word);
}

// The scan runs upward from the bottom of the window and stops at the first
// word that no longer holds its pattern. That word is by construction the
// deepest disturbed one, which is stronger than walking down from the origin
// and stopping at the first undisturbed word: the latter assumes the disturbed
// region is contiguous, and this does not.
[[gnu::noinline]] std::uint32_t *deepest_below(std::uint32_t *top) noexcept
{
    for(std::uint32_t *word = window_bottom(); word != top; ++word)
    {
        if(*word != paint_word(word))
            return word;
    }
    return top;
}

// Every repaint is preceded by this walk, so the program's record keeps the
// deepest word disturbed since the startup paint across the regions' repaints.
void fold_program() noexcept
{
    if(program.top == nullptr)
        return;
    std::uint32_t *const deepest = deepest_below(program.top);
    if(deepest < program.deepest)
        program.deepest = deepest;
}

stack_reading reading_of(const stack_window &window, const std::uint32_t *deepest, std::uint32_t floor) noexcept
{
    const auto address = reinterpret_cast<std::uintptr_t>(deepest);
    const auto top     = reinterpret_cast<std::uintptr_t>(window.top);
    const auto bottom  = reinterpret_cast<std::uintptr_t>(window_bottom());
    const bool touched = address != top;
    const auto room    = static_cast<std::uint32_t>(stack_top() - bottom);
    const auto from    = touched ? static_cast<std::uint32_t>(stack_top() - address) : 0U;
    const auto used    = touched ? static_cast<std::uint32_t>(window.reference - address) : 0U;
    return {used, from, room - from, room, static_cast<std::uint32_t>(window.reference - top), static_cast<std::uint32_t>(top - bottom), floor, touched && address == bottom};
}

}

void paint_stack_at_startup() noexcept
{
    if constexpr(kStackWatermarkPresent)
    {
        const std::uint32_t floor = close_stack_window(open_stack_window(), 0).used;
        const stack_window window = open_stack_window();
        program                   = {window.top, window.top, floor};
    }
}

[[gnu::noinline]] stack_window open_stack_window() noexcept
{
    fold_program();
    std::uint32_t anchor        = 0;
    const std::uintptr_t origin = reinterpret_cast<std::uintptr_t>(&anchor) & ~std::uintptr_t{sizeof(std::uint32_t) - 1};
    std::uint32_t *const top    = reinterpret_cast<std::uint32_t *>(origin - kGapBytes);
    paint_below(top);
    return {origin, top};
}

stack_reading close_stack_window(const stack_window &window, std::uint32_t floor) noexcept
{
    return reading_of(window, deepest_below(window.top), floor);
}

stack_reading read_program_stack() noexcept
{
    fold_program();
    return reading_of({stack_top(), program.top}, program.deepest, program.floor);
}

void report_stack(const char *region, const stack_reading &r)
{
    if constexpr(!kStackWatermarkPresent)
        return;
    std::printf("[stack] region=%s ", region);
    if(r.floor != 0 || r.saturated)
        std::printf("used=withheld from_top=withheld free=withheld");
    else
        std::printf("used=%" PRIu32 " from_top=%" PRIu32 " free=%" PRIu32, r.used, r.from_top, r.free);
    std::printf(" of=%" PRIu32 " bytes gap=%" PRIu32 " window=%" PRIu32 " floor=%" PRIu32 " saturated=%s\n", r.room, r.gap, r.window, r.floor, r.saturated ? "yes" : "no");
}

stack_step_region::stack_step_region(std::size_t first, std::size_t last)
        : reading_{}
        , window_{0, nullptr}
        , floor_(0)
        , first_(first)
        , last_(last)
{
}

void stack_step_region::enter(std::size_t k) noexcept
{
    if(k != first_)
        return;
    floor_  = close_stack_window(open_stack_window(), 0).used;
    window_ = open_stack_window();
}

void stack_step_region::leave(std::size_t k) noexcept
{
    if(k == last_)
        reading_ = close_stack_window(window_, floor_);
}

const stack_reading &stack_step_region::reading() const noexcept
{
    return reading_;
}

}
