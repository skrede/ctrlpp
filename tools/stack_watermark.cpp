/// @file
/// @brief Instrument: the whole-chain runtime stack watermark of a library hot
///        path, measured by painting a region below the current frame and
///        reading back the deepest disturbed word.
///
/// The `allocation-free?` column of the real-time safety matrix is a HEAP
/// statement. A no-allocation test counts the heap and is blind to the stack,
/// so a path proved allocation-free still carries an unmeasured frame cost, and
/// on a four-to-sixteen-kibibyte task stack the frame cost is what breaks a
/// caller. This program measures that cost directly. It is driven across state
/// dimensions by `stack_watermark.sh` beside it, which also takes the
/// complementary per-function frame measurement and prints both side by side.
///
/// **Two measurements exist and neither substitutes for the other.** The
/// compiler's per-function frame report attributes nothing to callees, so a sum
/// over named functions is a LOWER BOUND that excludes every library frame
/// beneath them. The runtime watermark below covers the whole chain, including
/// every Eigen and libm frame the report cannot see.
///
/// The procedure, which is the one the published Riccati figures were taken
/// with:
///
///   1. Spawn a thread whose stack is set to 64 MiB through the thread
///      attribute rather than inherited and hoped for. The painted region must
///      fit below the current frame without reaching the guard page, which is
///      why the stack is sixteen times the painted region.
///   2. Inside the thread, take the current frame's address as the origin.
///   3. Paint 4 MiB below the origin with a 64-bit pattern.
///   4. Run the chain under measurement. Construction and one warm-up call
///      happen OUTSIDE the painted window, because construction is offline and
///      only the hot path is being measured -- the same split the library's
///      no-allocation tests already make.
///   5. Walk the painted region and report the deepest word that no longer
///      holds its pattern. That offset in bytes from the origin is the
///      watermark.
///   6. Repeat with NO call in between to obtain the harness floor, and report
///      it beside every figure. A nonzero floor means the figure is not
///      attributable to the call, and every row derived from it would inherit
///      that.
///
/// **The per-word pattern includes the word's own address, and that is a
/// refinement over the original procedure rather than the original.** The
/// written record of the first measurement says a 64-bit pattern was used and
/// does not say which; a fixed constant can be matched by a word of legitimate
/// data, which reads as undisturbed and truncates the measurement downward.
/// Exclusive-oring the word's address into the pattern makes a collision
/// require the datum to equal a different value at every address, and costs
/// nothing. Whether the original harness did this is not known, so a figure
/// reproduced here is reproduced under a pattern that is at least as strict.
///
/// **The corpus is named because a figure that fails to reproduce has the
/// corpus as one of its candidate causes.** The discrete Riccati row is
/// measured over the discrete damped chain the internal Riccati benchmark
/// builds: a forward-Euler step of the continuous damped chain, `-0.5` on the
/// diagonal and `1.0` on the superdiagonal at a step of `0.01`, one input per
/// group of states, `Q = I` and `R = 0.1 I`. That weighting puts the weight
/// scale at exactly one, so equilibration is inactive and the entry point runs
/// its acceptance check once.
///
/// Building: standalone, no test framework and no build system. Eigen enters as
/// a system include, the way the library's own build treats it, so the only
/// diagnostics a build reports are this file's own.
///
///     g++ -std=c++20 -O2 -fno-exceptions -fno-rtti -pthread
///         -I lib/ctrlpp/include -isystem /usr/include/eigen3
///         -DWATERMARK_STATE_DIMENSION=4
///         tools/stack_watermark.cpp -o /tmp/stack_watermark
///
/// The measured row and its dimensions are compile-time, so the same instrument
/// serves further rows without being rewritten and so the optimizer sees the
/// chain exactly as a caller at that dimension would.
///
/// Output is one line carrying the row, the dimensions, the scalar type, the
/// whole-chain watermark in bytes and the harness floor in bytes. The floor is
/// on the same line as the figure by construction: a watermark reported without
/// its floor is not a measurement of anything.

#include "ctrlpp/control/dare.h"

#include <Eigen/Core>

#include <atomic>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <pthread.h>

/// The row under measurement. One selector per hot path the matrix publishes.
#define WATERMARK_ROW_RICCATI_DISCRETE 1

#ifndef WATERMARK_ROW
#define WATERMARK_ROW WATERMARK_ROW_RICCATI_DISCRETE
#endif

#ifndef WATERMARK_STATE_DIMENSION
#define WATERMARK_STATE_DIMENSION 2
#endif

#ifndef WATERMARK_INPUT_DIMENSION
#define WATERMARK_INPUT_DIMENSION 1
#endif

#if WATERMARK_ROW != WATERMARK_ROW_RICCATI_DISCRETE
#error "WATERMARK_ROW names no row this instrument implements"
#endif

namespace {

/// Sixteen times the painted region. The paint must fit below the current frame
/// without reaching the guard page, and a thread stack is the only stack whose
/// size this program gets to choose.
constexpr std::size_t thread_stack_bytes = std::size_t{64} * 1024 * 1024;

/// The painted window. Every published figure is three orders of magnitude
/// below this, so a watermark that reached the bottom would be reported as a
/// saturation rather than as a number.
constexpr std::size_t painted_bytes = std::size_t{4} * 1024 * 1024;

/// The top of the painted window is held this far below the origin so that the
/// painting and walking routines' own frames are never painted over. Anything
/// shallower than the gap is therefore invisible to the measurement, which is
/// harmless because the reported figure is the DEEPEST disturbed word, and it
/// is what makes the harness floor read zero. A floor that does not read zero
/// says the gap is too small for this build's harness frames.
constexpr std::size_t harness_gap_bytes = 1024;

/// The pattern's fixed half. Exclusive-ored with each word's own address, so no
/// single value written by the chain can pass as undisturbed at more than one
/// address.
constexpr std::uint64_t paint_seed = 0x5FC3A96E1D47B208ULL;

/// The chain's result is stored here so the optimizer cannot delete the call
/// whose stack cost is the whole point of the program.
std::atomic<std::uint64_t> chain_sink{0};

auto paint_word(const std::uint64_t *address) noexcept -> std::uint64_t
{
    return paint_seed ^ static_cast<std::uint64_t>(reinterpret_cast<std::uintptr_t>(address));
}

/// The origin, rounded down to the word grid so that every painted address is
/// aligned and the gap and window sizes stay exact multiples of the word.
auto aligned_origin(unsigned char *frame) noexcept -> unsigned char *
{
    const auto address = reinterpret_cast<std::uintptr_t>(frame);
    return reinterpret_cast<unsigned char *>(address & ~std::uintptr_t{7});
}

[[gnu::noinline]] void paint_region(unsigned char *origin) noexcept
{
    auto *const top    = reinterpret_cast<std::uint64_t *>(origin - harness_gap_bytes);
    auto *const bottom = reinterpret_cast<std::uint64_t *>(origin - painted_bytes);
    for(std::uint64_t *word = bottom; word != top; ++word)
        *word = paint_word(word);
}

/// The scan runs upward from the bottom of the window and stops at the first
/// word that no longer holds its pattern. That word is by construction the
/// deepest disturbed one, which is stronger than walking down from the origin
/// and stopping at the first undisturbed word: the latter assumes the disturbed
/// region is contiguous, and this does not.
[[gnu::noinline]] auto deepest_disturbed(unsigned char *origin) noexcept -> std::size_t
{
    const auto *const top    = reinterpret_cast<const std::uint64_t *>(origin - harness_gap_bytes);
    const auto *const bottom = reinterpret_cast<const std::uint64_t *>(origin - painted_bytes);
    for(const std::uint64_t *word = bottom; word != top; ++word)
    {
        if(*word != paint_word(word))
            return static_cast<std::size_t>(origin - reinterpret_cast<const unsigned char *>(word));
    }
    return 0;
}

/// One painted pass. With `invoke` false this is the harness floor: the same
/// paint and the same walk with no call in between.
///
/// The routine is kept out of line, and so is the chain it calls, because an
/// inlined chain would have its frame merged into this one -- above the origin,
/// where the paint cannot see it -- and the measurement would silently report
/// the harness instead of the chain.
template<typename Invocable>
[[gnu::noinline]] auto painted_pass(const Invocable &call, bool invoke) noexcept -> std::size_t
{
    unsigned char anchor       = 0;
    unsigned char *const origin = aligned_origin(&anchor);

    paint_region(origin);
    std::atomic_signal_fence(std::memory_order_seq_cst);

    if(invoke)
        call();

    std::atomic_signal_fence(std::memory_order_seq_cst);
    return deepest_disturbed(origin);
}

template<std::size_t NX, std::size_t NU>
struct discrete_damped_chain
{
    Eigen::Matrix<double, int(NX), int(NX)> A;
    Eigen::Matrix<double, int(NX), int(NU)> B;
    Eigen::Matrix<double, int(NX), int(NX)> Q;
    Eigen::Matrix<double, int(NU), int(NU)> R;
};

/// Forward-Euler step of the continuous damped chain, identical to the corpus
/// the internal Riccati benchmark sweeps. Every eigenvalue of `A` sits at
/// `1 - 0.5 dt`, inside the unit disk, so the discrete equation is well posed at
/// every size and the whole chain including the acceptance check executes.
template<std::size_t NX, std::size_t NU>
auto build_discrete_damped_chain() -> discrete_damped_chain<NX, NU>
{
    constexpr int    n  = int(NX);
    constexpr int    nu = int(NU);
    constexpr double dt = 0.01;

    discrete_damped_chain<NX, NU> corpus{};

    corpus.A = Eigen::Matrix<double, n, n>::Identity();
    for(std::size_t i = 0; i < NX; ++i)
        corpus.A(int(i), int(i)) += dt * -0.5;
    for(std::size_t i = 0; i + 1 < NX; ++i)
        corpus.A(int(i), int(i + 1)) = dt * 1.0;

    corpus.B                = Eigen::Matrix<double, n, nu>::Zero();
    const std::size_t group = NX / NU;
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t last_row = ((j + 1) * group < NX ? (j + 1) * group : NX) - 1;
        corpus.B(int(last_row), int(j)) = dt;
    }

    corpus.Q = Eigen::Matrix<double, n, n>::Identity();
    corpus.R = 0.1 * Eigen::Matrix<double, nu, nu>::Identity();

    return corpus;
}

/// The measured chain: the discrete Riccati entry point, called exactly as a
/// caller at this dimension calls it. Kept out of line so its frame lies below
/// the origin where the paint can see it.
template<std::size_t NX, std::size_t NU>
[[gnu::noinline]] void run_discrete_riccati(const discrete_damped_chain<NX, NU> &corpus) noexcept
{
    auto solved = ctrlpp::dare<double, NX, NU>(corpus.A, corpus.B, corpus.Q, corpus.R);

    std::uint64_t bits = 1;
    if(solved)
    {
        const double reduction = solved->P.sum() + solved->K.sum();
        std::memcpy(&bits, &reduction, sizeof(bits));
    }
    chain_sink.store(bits, std::memory_order_relaxed);
}

struct measurement
{
    std::size_t watermark_bytes;
    std::size_t floor_bytes;
    bool        solved;
};

auto watermark_thread(void *argument) noexcept -> void *
{
    constexpr std::size_t nx = WATERMARK_STATE_DIMENSION;
    constexpr std::size_t nu = WATERMARK_INPUT_DIMENSION;

    auto *const result = static_cast<measurement *>(argument);

    // Construction and warm-up sit outside every painted window: construction is
    // offline and only the hot path is being measured.
    const auto corpus = build_discrete_damped_chain<nx, nu>();
    const auto call   = [&corpus]() noexcept { run_discrete_riccati<nx, nu>(corpus); };
    call();

    result->solved      = chain_sink.load(std::memory_order_relaxed) != 1;
    result->floor_bytes = painted_pass(call, false);
    result->watermark_bytes = painted_pass(call, true);

    return nullptr;
}

}

auto main() -> int
{
    pthread_attr_t attributes;
    if(pthread_attr_init(&attributes) != 0)
    {
        std::fprintf(stderr, "stack_watermark: thread attributes unavailable\n");
        return 1;
    }

    if(pthread_attr_setstacksize(&attributes, thread_stack_bytes) != 0)
    {
        std::fprintf(stderr, "stack_watermark: a %zu-byte thread stack was refused\n", thread_stack_bytes);
        static_cast<void>(pthread_attr_destroy(&attributes));
        return 1;
    }

    measurement result{0, 0, false};
    pthread_t   thread{};
    if(pthread_create(&thread, &attributes, &watermark_thread, &result) != 0)
    {
        std::fprintf(stderr, "stack_watermark: the measurement thread could not be created\n");
        static_cast<void>(pthread_attr_destroy(&attributes));
        return 1;
    }

    static_cast<void>(pthread_join(thread, nullptr));
    static_cast<void>(pthread_attr_destroy(&attributes));

    std::printf("row=riccati-discrete nx=%d nu=%d scalar=double solved=%s watermark_bytes=%zu floor_bytes=%zu\n",
                int{WATERMARK_STATE_DIMENSION},
                int{WATERMARK_INPUT_DIMENSION},
                result.solved ? "yes" : "no",
                result.watermark_bytes,
                result.floor_bytes);

    return 0;
}
