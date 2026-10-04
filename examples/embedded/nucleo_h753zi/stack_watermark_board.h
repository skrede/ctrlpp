#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_STACK_WATERMARK_BOARD_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_NUCLEO_H753ZI_STACK_WATERMARK_BOARD_H

#include <atomic>
#include <cstddef>
#include <cstdint>

namespace ctrlpp {

// The evidence image measures its stack. The timing image paints nothing, so
// the code it times runs without the instrument beside it.
inline constexpr bool kStackWatermarkPresent = CTRLPP_MCU_ALLOC_SENTINEL != 0;

struct stack_window
{
    std::uintptr_t reference;
    std::uint32_t *top;
};

// used is the depth below the region's reference, the host tables' quantity;
// from_top is the depth below the top of the stack, the figure a stack reserve
// must cover; free and room restate it in the other board leg's phrasing. A
// nonzero floor or a saturated window withholds all three.
struct stack_reading
{
    std::uint32_t used;
    std::uint32_t from_top;
    std::uint32_t free;
    std::uint32_t room;
    std::uint32_t gap;
    std::uint32_t window;
    std::uint32_t floor;
    bool saturated;
};

// Called from the C library's initialization hook, which runs before the
// static constructors, so the program's figure covers them.
void paint_stack_at_startup() noexcept;

// Paints below the caller's frame. The reference is the painting call's own
// origin, so a chain measured between open and close starts a few bytes above
// it and used understates that chain by those bytes; from_top does not.
stack_window open_stack_window() noexcept;

stack_reading close_stack_window(const stack_window &window, std::uint32_t floor) noexcept;

// Everything since the startup paint, as deep as any region has reached.
stack_reading read_program_stack() noexcept;

void report_stack(const char *region, const stack_reading &reading);

// A floor pass, then the measured pass: the same paint and the same walk, with
// nothing between them the first time. The call must not be inlined: a chain
// merged into this frame sits above the origin, where the paint cannot see it,
// and the figure would then report the harness.
template<class Call>
stack_reading painted_pass(Call &&call) noexcept
{
    const std::uint32_t floor = close_stack_window(open_stack_window(), 0).used;
    const stack_window window = open_stack_window();
    std::atomic_signal_fence(std::memory_order_seq_cst);
    call();
    std::atomic_signal_fence(std::memory_order_seq_cst);
    return close_stack_window(window, floor);
}

// A region that spans steps first through last of a run whose loop the caller
// cannot split, opened from the run's entry hook and closed from its exit hook.
class stack_step_region
{
public:
    stack_step_region(std::size_t first, std::size_t last);

    void enter(std::size_t k) noexcept;

    void leave(std::size_t k) noexcept;

    const stack_reading &reading() const noexcept;

private:
    stack_reading reading_;
    stack_window window_;
    std::uint32_t floor_;
    std::size_t first_;
    std::size_t last_;
};

}

#endif
