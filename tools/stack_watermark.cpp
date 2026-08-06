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
///   1. Spawn a thread whose stack is set through the thread attribute rather
///      than inherited and hoped for, at sixteen times the painted region. The
///      painted region must fit below the current frame without reaching the
///      guard page, which is what that relation buys.
///   2. Inside the thread, take the current frame's address as the origin.
///   3. Paint the window below the origin with a 64-bit pattern. The window is
///      4 MiB by default and is a compile-time selector, because it is a
///      ceiling on what can be reported and the deepest chains this instrument
///      is driven over reach it.
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
/// **The estimator rows are measured over the same plant**, so a filter figure
/// and a Riccati figure at the same state dimension are figures about one
/// system rather than two. The dynamics are that forward-Euler chain; the
/// measurement map takes output `i` from state `i mod NX`, which is a defined
/// map at every pair of dimensions including the ones where there are more
/// outputs than states; the attitude rows instead propagate a constant body
/// rate on SO(3) and measure the gravity direction, repeated to whatever output
/// dimension is asked for, because their state is a rotation and not a vector.
///
/// **The predictive controller row is measured over that same forward-Euler
/// chain as its prediction model**, with a state weight of `10 I`, an input
/// weight of `0.1 I` and no path or terminal constraint, so its only
/// constraints are the dynamics continuity and the initial state. It is the
/// only row whose build needs the optional nonlinear-programming backend, and
/// the backend's revision is part of its provenance: the frames are an output
/// of inlining that backend's templates, so figures taken against two revisions
/// of it are not comparable and the driver prints the revision it resolved.
///
/// **This row grids over what a CALLER CHOOSES and reports what that induces.**
/// The controller's decision dimension `NV = (NH + 1) * NX + NH * NU` and its
/// constraint bound `MaxM = NX * (NH + 1)` are derived from the state, input
/// and horizon dimensions rather than selected, so a grid over the derived pair
/// would describe configurations no caller can reach. Both are computed here
/// from the selectors and printed beside every figure.
///
/// **The input dimension is part of the provenance of every row here.** It is a
/// separate compile-time selector rather than a constant, it is printed on
/// every output line, and it moves the discrete Riccati watermark by more than
/// a thousand bytes at two states. A stack figure quoted without it is not
/// reproducible.
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
/// whole-chain watermark in bytes, the harness floor in bytes, the painted
/// window and whether the figure saturated it. The floor is on the same line as
/// the figure by construction: a watermark reported without its floor is not a
/// measurement of anything, and one that reached the bottom of its window is a
/// lower bound rather than a figure.

/// The row under measurement. One selector per hot path the matrix publishes.
///
/// The unchecked discrete Riccati row is the same chain with the acceptance
/// check taken off it, composed from the parts the entry point composes so that
/// the two differ in the check and in nothing else. It exists because the cost
/// of the check is only separable from the cost of the solve if both are
/// measured by one instrument in one run.
#define WATERMARK_ROW_RICCATI_DISCRETE           1
#define WATERMARK_ROW_RICCATI_DISCRETE_UNCHECKED 2
#define WATERMARK_ROW_KALMAN                     3
#define WATERMARK_ROW_EKF                        4
#define WATERMARK_ROW_UKF                        5
#define WATERMARK_ROW_MANIFOLD_UKF               6
#define WATERMARK_ROW_MEKF                       7
#define WATERMARK_ROW_NMPC_STATIC                8

#ifndef WATERMARK_ROW
#define WATERMARK_ROW WATERMARK_ROW_RICCATI_DISCRETE
#endif

#include "ctrlpp/control/dare.h"

#include "ctrlpp/estimation/ekf.h"
#include "ctrlpp/estimation/ukf.h"
#include "ctrlpp/estimation/mekf.h"
#include "ctrlpp/estimation/kalman.h"
#include "ctrlpp/estimation/manifold_ukf.h"

#include "ctrlpp/lie/so3.h"

#include "ctrlpp/model/state_space.h"

// The predictive controller row alone needs the optional nonlinear-programming
// backend, and pulling its headers into every row would put that backend's
// instantiation cost on rows that never call it.
#if WATERMARK_ROW == WATERMARK_ROW_NMPC_STATIC
#include "ctrlpp/nmpc.h"

#include "ctrlpp/mpc/argmin_solver.h"
#endif

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <atomic>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <pthread.h>

/// The state dimension the caller chooses. It is `NX` for the Riccati rows and
/// for the three vector-state filters, and the BIAS dimension `NB` for the
/// multiplicative error-state filter, whose error state is `3 + NB`. The
/// manifold filter has no such selector at all: its state is a rotation, so its
/// state dimension is three by construction and is not the caller's to choose.
#ifndef WATERMARK_STATE_DIMENSION
#define WATERMARK_STATE_DIMENSION 2
#endif

#ifndef WATERMARK_MEASUREMENT_DIMENSION
#define WATERMARK_MEASUREMENT_DIMENSION 1
#endif

#ifndef WATERMARK_INPUT_DIMENSION
#define WATERMARK_INPUT_DIMENSION 1
#endif

/// The prediction horizon, which is the predictive controller row's third
/// caller-chosen dimension and is read by no other row.
#ifndef WATERMARK_HORIZON
#define WATERMARK_HORIZON 5
#endif

#if WATERMARK_ROW < WATERMARK_ROW_RICCATI_DISCRETE || WATERMARK_ROW > WATERMARK_ROW_NMPC_STATIC
#error "WATERMARK_ROW names no row this instrument implements"
#endif

namespace {

/// The painted window and the thread stack that has to hold it, both in
/// mebibytes.
///
/// **BOTH ARE SELECTORS RATHER THAN CONSTANTS, BECAUSE THE WINDOW IS A CEILING
/// ON WHAT THIS INSTRUMENT CAN REPORT AND THE INTERIOR OF THE GRID REACHES IT.**
/// A chain deeper than the window disturbs the window's bottom word, and the
/// walk then returns the window's own size rather than the chain's depth: a
/// number that looks exactly like a measurement and is not one. That case is
/// reported below as a saturation rather than left for a reader to notice, and
/// the window is raised for the dimensions that need it rather than the figure
/// being silently truncated.
///
/// **Raising the window cannot move a figure that was already resolved.** The
/// walk runs upward from the bottom and returns the DEEPEST disturbed word, so a
/// deeper bottom can only add words the chain never touched, which the walk
/// passes over. That is a claim about the instrument rather than about the
/// library, so the driver checks it by re-measuring points the published tables
/// already carry instead of asserting it here.
#ifndef WATERMARK_PAINTED_MIB
#define WATERMARK_PAINTED_MIB 4
#endif

#ifndef WATERMARK_THREAD_STACK_MIB
#define WATERMARK_THREAD_STACK_MIB 64
#endif

constexpr std::size_t painted_bytes = std::size_t{WATERMARK_PAINTED_MIB} * 1024 * 1024;

/// Sixteen times the painted region. The paint must fit below the current frame
/// without reaching the guard page, and a thread stack is the only stack whose
/// size this program gets to choose. The relation is the one the published
/// figures were taken under, and it is asserted rather than left as a default a
/// selector could quietly drop below.
constexpr std::size_t thread_stack_bytes = std::size_t{WATERMARK_THREAD_STACK_MIB} * 1024 * 1024;

static_assert(thread_stack_bytes >= 16 * painted_bytes,
              "the thread stack must be at least sixteen times the painted window");

/// The top of the painted window is held this far below the origin so that the
/// painting and walking routines' own frames are never painted over. Anything
/// shallower than the gap is therefore invisible to the measurement, and it is
/// what makes the harness floor read zero. A floor that does not read zero says
/// the gap is too small for this build's harness frames.
///
/// **The gap is a RESOLUTION FLOOR, and it binds on the cheap rows rather than
/// on the expensive ones.** A chain whose whole frame cost is shallower than
/// the gap disturbs nothing the walk can see and reports zero, which is the
/// same value the harness floor reports. The Riccati figures sit three
/// thousand bytes and more below the gap so it never bound there; the estimator
/// chains at their smallest dimensions do not, so the gap is exposed as a
/// compile-time selector and swept rather than assumed.
///
/// Shrinking it cannot move a figure that was already resolved: the walk runs
/// upward from the bottom of the window and returns the DEEPEST disturbed word,
/// so extending the window's top can only add shallower words the return value
/// ignores. That invariance is a claim about the instrument, so the driver
/// checks it by re-running the Riccati rows at both gaps rather than asserting
/// it here.
#ifndef WATERMARK_HARNESS_GAP_BYTES
#define WATERMARK_HARNESS_GAP_BYTES 1024
#endif

constexpr std::size_t harness_gap_bytes = WATERMARK_HARNESS_GAP_BYTES;

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

/// The step of the plant every vector-state row is measured over, shared so that
/// a filter figure and a Riccati figure at the same state dimension are figures
/// about ONE system. Forward-Euler step of the continuous damped chain,
/// identical to the corpus the internal Riccati benchmark sweeps. Every
/// eigenvalue sits at `1 - 0.5 dt`, inside the unit disk, so the discrete
/// Riccati equation is well posed at every size and the whole chain including
/// the acceptance check executes.
template<std::size_t NX>
auto chain_state_matrix() -> Eigen::Matrix<double, int(NX), int(NX)>
{
    constexpr int    n  = int(NX);
    constexpr double dt = 0.01;

    Eigen::Matrix<double, n, n> A = Eigen::Matrix<double, n, n>::Identity();
    for(std::size_t i = 0; i < NX; ++i)
        A(int(i), int(i)) += dt * -0.5;
    for(std::size_t i = 0; i + 1 < NX; ++i)
        A(int(i), int(i + 1)) = dt * 1.0;
    return A;
}

/// One input per group of states, acting on the last state of its group. The
/// group size is clamped at one so the map stays defined when there are more
/// inputs than states; at the input dimensions this instrument is driven at,
/// which are never above the state dimension, the clamp is inactive and the
/// matrix is the one the Riccati benchmark builds.
template<std::size_t NX, std::size_t NU>
auto chain_input_matrix() -> Eigen::Matrix<double, int(NX), int(NU)>
{
    constexpr int    n  = int(NX);
    constexpr int    nu = int(NU);
    constexpr double dt = 0.01;

    Eigen::Matrix<double, n, nu> B = Eigen::Matrix<double, n, nu>::Zero();
    const std::size_t           group = NU <= NX ? NX / NU : std::size_t{1};
    for(std::size_t j = 0; j < NU; ++j)
    {
        const std::size_t reach    = (j + 1) * group;
        const std::size_t last_row = (reach < NX ? reach : NX) - 1;
        B(int(last_row), int(j)) = dt;
    }
    return B;
}

/// Output `i` reads state `i mod NX`. The modulus is what makes the map defined
/// at every pair of dimensions, including the held-dimension sweeps where the
/// output dimension runs past the state dimension; there it duplicates rows,
/// which leaves the innovation covariance nonsingular because the measurement
/// noise is the identity.
template<std::size_t NX, std::size_t NY>
auto chain_output_matrix() -> Eigen::Matrix<double, int(NY), int(NX)>
{
    Eigen::Matrix<double, int(NY), int(NX)> H = Eigen::Matrix<double, int(NY), int(NX)>::Zero();
    for(std::size_t i = 0; i < NY; ++i)
        H(int(i), int(i % NX)) = 1.0;
    return H;
}

template<std::size_t NX, std::size_t NU>
struct discrete_damped_chain
{
    Eigen::Matrix<double, int(NX), int(NX)> A;
    Eigen::Matrix<double, int(NX), int(NU)> B;
    Eigen::Matrix<double, int(NX), int(NX)> Q;
    Eigen::Matrix<double, int(NU), int(NU)> R;
};

template<std::size_t NX, std::size_t NU>
auto build_discrete_damped_chain() -> discrete_damped_chain<NX, NU>
{
    constexpr int nu = int(NU);

    discrete_damped_chain<NX, NU> corpus{};
    corpus.A = chain_state_matrix<NX>();
    corpus.B = chain_input_matrix<NX, NU>();
    corpus.Q = Eigen::Matrix<double, int(NX), int(NX)>::Identity();
    corpus.R = 0.1 * Eigen::Matrix<double, nu, nu>::Identity();
    return corpus;
}

/// The measured chain: the discrete Riccati entry point, called exactly as a
/// caller at this dimension calls it. Kept out of line so its frame lies below
/// the origin where the paint can see it.
template<std::size_t NX, std::size_t NU>
[[gnu::noinline]] void run_discrete_riccati(discrete_damped_chain<NX, NU> &corpus) noexcept
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

/// The same solve with the acceptance check taken off it, composed from the
/// parts the entry point composes and stopping where the entry point begins
/// verifying. The two rows therefore differ in the check and in nothing else,
/// which is the only way the cost of the check is separable from the cost of the
/// solve.
///
/// This composition mirrors the entry point rather than calling it, so it is the
/// one place in this file that can drift from the library without failing to
/// compile. The corpus puts the weight scale at exactly one, so the entry
/// point's equilibration branch is inactive and there is nothing on that branch
/// to mirror.
template<std::size_t NX, std::size_t NU>
[[gnu::noinline]] void run_discrete_riccati_unchecked(discrete_damped_chain<NX, NU> &corpus) noexcept
{
    std::uint64_t bits = 1;

    auto operands = ctrlpp::detail::factor_dare_symplectic_operands<double, NX, NU>(corpus.A, corpus.B, corpus.R);
    if(operands)
    {
        auto symplectic = ctrlpp::detail::build_dare_symplectic<double, NX>(corpus.A, corpus.Q, *operands);
        if(symplectic)
        {
            auto solved = ctrlpp::detail::dare_solve_from_symplectic<double, NX, NU>(*symplectic);
            if(solved)
            {
                const double reduction = solved->P.sum();
                std::memcpy(&bits, &reduction, sizeof(bits));
            }
        }
    }
    chain_sink.store(bits, std::memory_order_relaxed);
}

/// The dynamics the two Jacobian-based filters linearize and the sigma-point
/// filter propagates. Linear, so the linearization is exact and no filter is
/// measured on a chain a modeling error made longer or shorter than another's.
template<std::size_t NX, std::size_t NU>
struct chain_dynamics
{
    Eigen::Matrix<double, int(NX), int(NX)> F;
    Eigen::Matrix<double, int(NX), int(NU)> G;

    auto operator()(const ctrlpp::Vector<double, NX> &x, const ctrlpp::Vector<double, NU> &u) const
        -> ctrlpp::Vector<double, NX>
    {
        return ctrlpp::Vector<double, NX>{F * x + G * u};
    }

    auto jacobian_x(const ctrlpp::Vector<double, NX> &, const ctrlpp::Vector<double, NU> &) const
        -> Eigen::Matrix<double, int(NX), int(NX)>
    {
        return F;
    }

    auto jacobian_u(const ctrlpp::Vector<double, NX> &, const ctrlpp::Vector<double, NU> &) const
        -> Eigen::Matrix<double, int(NX), int(NU)>
    {
        return G;
    }
};

template<std::size_t NX, std::size_t NY>
struct chain_measurement
{
    Eigen::Matrix<double, int(NY), int(NX)> H;

    auto operator()(const ctrlpp::Vector<double, NX> &x) const -> ctrlpp::Vector<double, NY>
    {
        return ctrlpp::Vector<double, NY>{H * x};
    }

    auto jacobian(const ctrlpp::Vector<double, NX> &) const -> Eigen::Matrix<double, int(NY), int(NX)>
    {
        return H;
    }
};

/// A constant body rate integrated on SO(3). The attitude rows have no vector
/// state to give the damped chain to, so their plant is this instead, and the
/// two families of figures are not comparable at equal state dimension for that
/// reason.
struct body_rate_dynamics
{
    double dt = 0.01;

    auto operator()(const Eigen::Quaternion<double> &q, const ctrlpp::Vector<double, 3> &omega) const
        -> Eigen::Quaternion<double>
    {
        const ctrlpp::Vector<double, 3> increment = (omega * dt).eval();
        return (q * ctrlpp::so3::exp(increment)).normalized();
    }
};

/// The gravity direction in body axes, repeated to whatever output dimension is
/// asked for. Repetition rather than truncation, so the map is defined for an
/// output dimension above three as well as below it.
template<std::size_t NY>
auto gravity_direction(const Eigen::Quaternion<double> &q) -> ctrlpp::Vector<double, NY>
{
    const Eigen::Matrix<double, 3, 3> rotation = q.toRotationMatrix();
    const ctrlpp::Vector<double, 3>   down     = rotation.transpose().col(2);

    ctrlpp::Vector<double, NY> z;
    for(std::size_t i = 0; i < NY; ++i)
        z(int(i)) = down(int(i % 3));
    return z;
}

template<std::size_t NY>
struct attitude_measurement
{
    auto operator()(const Eigen::Quaternion<double> &q) const -> ctrlpp::Vector<double, NY>
    {
        return gravity_direction<NY>(q);
    }
};

template<std::size_t NB, std::size_t NY>
struct attitude_bias_measurement
{
    auto operator()(const Eigen::Quaternion<double> &q, const ctrlpp::Vector<double, NB> &) const
        -> ctrlpp::Vector<double, NY>
    {
        return gravity_direction<NY>(q);
    }
};

/// A filter row's carried objects. The filter is held as the factory returned
/// it, so a configuration the factory refuses is reported as an unsolved row
/// rather than measured as if it had run.
template<typename Filter, typename Error, std::size_t NI, std::size_t NY>
struct filter_chain
{
    ctrlpp::expected<Filter, Error> filter;
    ctrlpp::Vector<double, NI>      input;
    ctrlpp::Vector<double, NY>      measurement_sample;
};

/// The painted window covers the prediction and the update as ONE block rather
/// than each separately with a maximum taken afterwards. That is what the
/// allocation guards already do, and the deeper of the two is what a task stack
/// has to hold in either case.
template<typename Chain>
[[gnu::noinline]] void run_filter(Chain &chain) noexcept
{
    std::uint64_t bits = 1;
    if(chain.filter)
    {
        auto &filter = *chain.filter;
        filter.predict(chain.input);
        const auto stepped = filter.update(chain.measurement_sample);
        if(stepped)
        {
            const auto   state     = filter.state();
            const double reduction = state.sum();
            std::memcpy(&bits, &reduction, sizeof(bits));
        }
    }
    chain_sink.store(bits, std::memory_order_relaxed);
}

template<std::size_t NX, std::size_t NU, std::size_t NY>
auto build_kalman_chain()
    -> filter_chain<ctrlpp::kalman_filter<double, NX, NU, NY>, ctrlpp::filter_error, NU, NY>
{
    ctrlpp::discrete_state_space<double, NX, NU, NY> system;
    system.A = chain_state_matrix<NX>();
    system.B = chain_input_matrix<NX, NU>();
    system.C = chain_output_matrix<NX, NY>();
    system.D = Eigen::Matrix<double, int(NY), int(NU)>::Zero();

    const ctrlpp::kalman_config<double, NX, NU, NY> config{
        .Q  = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01,
        .R  = ctrlpp::Matrix<double, NY, NY>::Identity(),
        .x0 = ctrlpp::Vector<double, NX>::Zero(),
        .P0 = ctrlpp::Matrix<double, NX, NX>::Identity() * 10.0};

    return {ctrlpp::kalman_filter<double, NX, NU, NY>::create(system, config),
            ctrlpp::Vector<double, NU>::Zero(),
            ctrlpp::Vector<double, NY>::Constant(0.1)};
}

template<std::size_t NX, std::size_t NU, std::size_t NY>
auto build_ekf_chain()
    -> filter_chain<ctrlpp::ekf<double, NX, NU, NY, chain_dynamics<NX, NU>, chain_measurement<NX, NY>>,
                    ctrlpp::filter_error, NU, NY>
{
    const chain_dynamics<NX, NU>    dynamics{chain_state_matrix<NX>(), chain_input_matrix<NX, NU>()};
    const chain_measurement<NX, NY> measurement_map{chain_output_matrix<NX, NY>()};

    const ctrlpp::ekf_config<double, NX, NU, NY> config{
        .Q  = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01,
        .R  = ctrlpp::Matrix<double, NY, NY>::Identity(),
        .x0 = ctrlpp::Vector<double, NX>::Zero(),
        .P0 = ctrlpp::Matrix<double, NX, NX>::Identity() * 10.0};

    return {ctrlpp::ekf<double, NX, NU, NY, chain_dynamics<NX, NU>, chain_measurement<NX, NY>>::create(
                dynamics, measurement_map, config),
            ctrlpp::Vector<double, NU>::Zero(),
            ctrlpp::Vector<double, NY>::Constant(0.1)};
}

template<std::size_t NX, std::size_t NU, std::size_t NY>
auto build_ukf_chain()
    -> filter_chain<ctrlpp::ukf<double, NX, NU, NY, chain_dynamics<NX, NU>, chain_measurement<NX, NY>>,
                    ctrlpp::filter_error, NU, NY>
{
    const chain_dynamics<NX, NU>    dynamics{chain_state_matrix<NX>(), chain_input_matrix<NX, NU>()};
    const chain_measurement<NX, NY> measurement_map{chain_output_matrix<NX, NY>()};

    const ctrlpp::ukf_config<double, NX, NU, NY> config{
        .Q  = ctrlpp::Matrix<double, NX, NX>::Identity() * 0.01,
        .R  = ctrlpp::Matrix<double, NY, NY>::Identity(),
        .x0 = ctrlpp::Vector<double, NX>::Zero(),
        .P0 = ctrlpp::Matrix<double, NX, NX>::Identity() * 10.0};

    return {ctrlpp::ukf<double, NX, NU, NY, chain_dynamics<NX, NU>, chain_measurement<NX, NY>>::create(
                dynamics, measurement_map, config),
            ctrlpp::Vector<double, NU>::Zero(),
            ctrlpp::Vector<double, NY>::Constant(0.1)};
}

template<std::size_t NY>
auto build_manifold_ukf_chain()
    -> filter_chain<ctrlpp::manifold_ukf<double, NY, body_rate_dynamics, attitude_measurement<NY>>,
                    ctrlpp::filter_error, 3, NY>
{
    ctrlpp::manifold_ukf_config<double, NY> config;
    config.Q *= 1e-6;

    ctrlpp::Vector<double, 3> rate;
    rate << 0.01, -0.02, 0.03;

    return {ctrlpp::manifold_ukf<double, NY, body_rate_dynamics, attitude_measurement<NY>>::create(
                body_rate_dynamics{}, attitude_measurement<NY>{}, config),
            rate,
            gravity_direction<NY>(Eigen::Quaternion<double>::Identity())};
}

template<std::size_t NB, std::size_t NY>
auto build_mekf_chain()
    -> filter_chain<ctrlpp::mekf<double, NB, NY, attitude_bias_measurement<NB, NY>>,
                    ctrlpp::filter_error, 3, NY>
{
    ctrlpp::mekf_config<double, NB, NY> config;
    config.Q *= 1e-6;

    ctrlpp::Vector<double, 3> rate;
    rate << 0.01, -0.02, 0.03;

    return {ctrlpp::mekf<double, NB, NY, attitude_bias_measurement<NB, NY>>::create(
                attitude_bias_measurement<NB, NY>{}, config),
            rate,
            gravity_direction<NY>(Eigen::Quaternion<double>::Identity())};
}

#if WATERMARK_ROW == WATERMARK_ROW_NMPC_STATIC

/// The predictive controller's carried objects, and the two DERIVED dimensions
/// computed here from the three the caller chose.
///
/// `NV = (NH + 1) * NX + NH * NU` counts one state block per horizon node
/// including the initial one, plus one input block per interval;
/// `MaxM = NX * (NH + 1)` counts the dynamics-continuity rows plus the
/// initial-state rows, which are the only equalities this pose carries. There
/// are no path or terminal constraints and no slack, so the bound is met with
/// equality and the solver's per-call multiplier storage stays inline.
template<std::size_t NX, std::size_t NU, std::size_t NH>
struct controller_chain
{
    static constexpr int decision_dimension   = static_cast<int>((NH + 1) * NX + NH * NU);
    static constexpr int constraint_dimension = static_cast<int>(NX * (NH + 1));

    using solver_type =
        ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, decision_dimension, constraint_dimension>;
    using controller_type =
        ctrlpp::nmpc_static<double, NX, NU, NH, solver_type, chain_dynamics<NX, NU>>;

    chain_dynamics<NX, NU>     dynamics;
    controller_type            controller;
    ctrlpp::Vector<double, NX> state;
};

/// The measured call: one steady-state solve followed by the caller's own plant
/// step. That pair is the loop body the shipped allocation test arms, so the
/// two pieces of evidence for this row cover the same block rather than two
/// different ones. The plant step is a fixed-size linear map and the solve is
/// everything else.
template<std::size_t NX, std::size_t NU, std::size_t NH>
[[gnu::noinline]] void run_controller(controller_chain<NX, NU, NH> &chain) noexcept
{
    std::uint64_t bits = 1;

    const auto stepped = chain.controller.solve(chain.state);
    if(stepped)
    {
        const double reduction = stepped->input.sum();
        std::memcpy(&bits, &reduction, sizeof(bits));
        chain.state = chain.dynamics(chain.state, stepped->input);
    }
    chain_sink.store(bits, std::memory_order_relaxed);
}

/// Construction and the walk into steady state both sit outside every painted
/// window. The first solves construct the solver's state and flush its one-time
/// lazy instantiation, and each of them advances the plant, so the call the
/// paint sees is a warm-started solve at a state the previous solve did not
/// answer -- a steady-state step rather than a converged re-solve at an
/// unchanged pose, which would report a shorter chain than a caller runs.
template<std::size_t NX, std::size_t NU, std::size_t NH>
auto build_controller_chain(int warm_up_solves) -> controller_chain<NX, NU, NH>
{
    const chain_dynamics<NX, NU> dynamics{chain_state_matrix<NX>(), chain_input_matrix<NX, NU>()};

    ctrlpp::nmpc_config<double, NX, NU> config;
    config.horizon = static_cast<int>(NH);
    config.Q       = ctrlpp::Matrix<double, NX, NX>::Identity() * 10.0;
    config.R       = ctrlpp::Matrix<double, NU, NU>::Identity() * 0.1;

    ctrlpp::Vector<double, NX> initial_state = ctrlpp::Vector<double, NX>::Zero();
    initial_state(0)                         = 1.0;

    controller_chain<NX, NU, NH> chain{
        dynamics,
        typename controller_chain<NX, NU, NH>::controller_type{dynamics, config},
        initial_state};

    for(int solve = 0; solve < warm_up_solves; ++solve)
        run_controller(chain);
    return chain;
}

#endif

constexpr std::size_t selected_state_dimension       = WATERMARK_STATE_DIMENSION;
constexpr std::size_t selected_measurement_dimension = WATERMARK_MEASUREMENT_DIMENSION;
constexpr std::size_t selected_input_dimension       = WATERMARK_INPUT_DIMENSION;
constexpr std::size_t selected_horizon               = WATERMARK_HORIZON;

/// Each row names itself, states the dimensions it was ACTUALLY instantiated at
/// -- which for the attitude rows is not always the dimension asked for, since
/// their state is a rotation -- and supplies a builder and a runner. Everything
/// below the selection is common, so no row carries a measurement procedure of
/// its own.
#if WATERMARK_ROW == WATERMARK_ROW_RICCATI_DISCRETE

constexpr const char *row_name                        = "riccati-discrete";
constexpr std::size_t instantiated_state_dimension    = selected_state_dimension;
constexpr std::size_t instantiated_input_dimension    = selected_input_dimension;
constexpr std::size_t instantiated_output_dimension   = 0;

auto build_measured_chain()
{
    return build_discrete_damped_chain<selected_state_dimension, selected_input_dimension>();
}

template<typename Chain>
void run_measured_chain(Chain &chain) noexcept
{
    run_discrete_riccati<selected_state_dimension, selected_input_dimension>(chain);
}

#elif WATERMARK_ROW == WATERMARK_ROW_RICCATI_DISCRETE_UNCHECKED

constexpr const char *row_name                        = "riccati-discrete-unchecked";
constexpr std::size_t instantiated_state_dimension    = selected_state_dimension;
constexpr std::size_t instantiated_input_dimension    = selected_input_dimension;
constexpr std::size_t instantiated_output_dimension   = 0;

auto build_measured_chain()
{
    return build_discrete_damped_chain<selected_state_dimension, selected_input_dimension>();
}

template<typename Chain>
void run_measured_chain(Chain &chain) noexcept
{
    run_discrete_riccati_unchecked<selected_state_dimension, selected_input_dimension>(chain);
}

#elif WATERMARK_ROW == WATERMARK_ROW_KALMAN

constexpr const char *row_name                        = "kalman-filter";
constexpr std::size_t instantiated_state_dimension    = selected_state_dimension;
constexpr std::size_t instantiated_input_dimension    = selected_input_dimension;
constexpr std::size_t instantiated_output_dimension   = selected_measurement_dimension;

auto build_measured_chain()
{
    return build_kalman_chain<selected_state_dimension, selected_input_dimension, selected_measurement_dimension>();
}

template<typename Chain>
void run_measured_chain(Chain &chain) noexcept
{
    run_filter(chain);
}

#elif WATERMARK_ROW == WATERMARK_ROW_EKF

constexpr const char *row_name                        = "ekf";
constexpr std::size_t instantiated_state_dimension    = selected_state_dimension;
constexpr std::size_t instantiated_input_dimension    = selected_input_dimension;
constexpr std::size_t instantiated_output_dimension   = selected_measurement_dimension;

auto build_measured_chain()
{
    return build_ekf_chain<selected_state_dimension, selected_input_dimension, selected_measurement_dimension>();
}

template<typename Chain>
void run_measured_chain(Chain &chain) noexcept
{
    run_filter(chain);
}

#elif WATERMARK_ROW == WATERMARK_ROW_UKF

constexpr const char *row_name                        = "ukf";
constexpr std::size_t instantiated_state_dimension    = selected_state_dimension;
constexpr std::size_t instantiated_input_dimension    = selected_input_dimension;
constexpr std::size_t instantiated_output_dimension   = selected_measurement_dimension;

auto build_measured_chain()
{
    return build_ukf_chain<selected_state_dimension, selected_input_dimension, selected_measurement_dimension>();
}

template<typename Chain>
void run_measured_chain(Chain &chain) noexcept
{
    run_filter(chain);
}

#elif WATERMARK_ROW == WATERMARK_ROW_MANIFOLD_UKF

/// The state dimension selector is IGNORED here and the row reports three,
/// because the manifold filter's state is a rotation. Reporting the requested
/// value would publish a dimension the caller cannot choose.
constexpr const char *row_name                        = "manifold-ukf";
constexpr std::size_t instantiated_state_dimension    = 3;
constexpr std::size_t instantiated_input_dimension    = 3;
constexpr std::size_t instantiated_output_dimension   = selected_measurement_dimension;

auto build_measured_chain()
{
    return build_manifold_ukf_chain<selected_measurement_dimension>();
}

template<typename Chain>
void run_measured_chain(Chain &chain) noexcept
{
    run_filter(chain);
}

#elif WATERMARK_ROW == WATERMARK_ROW_MEKF

/// The selector carries the BIAS dimension, and the error state the covariance
/// recursion runs at is three larger. Both are reported: the first is what the
/// caller picks and the second is what the frames are sized by.
constexpr const char *row_name                        = "mekf";
constexpr std::size_t instantiated_state_dimension    = selected_state_dimension;
constexpr std::size_t instantiated_input_dimension    = 3;
constexpr std::size_t instantiated_output_dimension   = selected_measurement_dimension;

auto build_measured_chain()
{
    return build_mekf_chain<selected_state_dimension, selected_measurement_dimension>();
}

template<typename Chain>
void run_measured_chain(Chain &chain) noexcept
{
    run_filter(chain);
}

#else

/// The three selectors this row reads are the state, input and horizon
/// dimensions, which is what a caller picks. There is no measurement dimension
/// to report and the two dimensions the frames are sized by are derived, so
/// they are printed as their own fields rather than folded into one of the
/// three above.
constexpr const char *row_name                        = "nmpc-static";
constexpr std::size_t instantiated_state_dimension    = selected_state_dimension;
constexpr std::size_t instantiated_input_dimension    = selected_input_dimension;
constexpr std::size_t instantiated_output_dimension   = 0;

using selected_controller_chain =
    controller_chain<selected_state_dimension, selected_input_dimension, selected_horizon>;

constexpr std::size_t instantiated_decision_dimension =
    static_cast<std::size_t>(selected_controller_chain::decision_dimension);
constexpr std::size_t instantiated_constraint_dimension =
    static_cast<std::size_t>(selected_controller_chain::constraint_dimension);

/// The steady-state chain has walked into steady state before it is returned;
/// the first-call chain has not, so the call the paint sees on it is the solve
/// that emplaces the solver.
constexpr int steady_state_warm_up_solves = 20;

auto build_measured_chain()
{
    return build_controller_chain<selected_state_dimension, selected_input_dimension, selected_horizon>(
        steady_state_warm_up_solves);
}

auto build_first_call_chain()
{
    return build_controller_chain<selected_state_dimension, selected_input_dimension, selected_horizon>(0);
}

/// Construction, measured as its own chain rather than reported as an absence.
///
/// It is offline and it is excluded from both solve figures above, exactly as
/// the allocation guards exclude it -- but the frame report puts the deepest
/// frame in the whole translation unit inside this constructor, larger than
/// either solve's whole chain, and a caller who constructs the controller on
/// the task's own stack pays it. Publishing the solve figures alone would leave
/// the largest of the three unmeasured.
///
/// Kept out of line for the same reason every other measured chain is: an
/// inlined body would put its locals in the painting routine's own frame, above
/// the origin, where the paint cannot see them.
[[gnu::noinline]] void run_construction() noexcept
{
    auto built = build_first_call_chain();

    std::uint64_t bits      = 0;
    const double  reduction = built.state.sum();
    std::memcpy(&bits, &reduction, sizeof(bits));
    chain_sink.store(bits, std::memory_order_relaxed);
}

template<typename Chain>
void run_measured_chain(Chain &chain) noexcept
{
    run_controller(chain);
}

#endif

struct measurement
{
    std::size_t watermark_bytes;
    std::size_t first_call_bytes;
    std::size_t construction_bytes;
    std::size_t floor_bytes;
    bool        solved;
};

/// A peak equal to the window's own size means the chain reached the bottom
/// word, so the walk returned where it stopped looking rather than where the
/// chain stopped writing. Such a figure is a LOWER BOUND and not a measurement,
/// and it is reported as one.
auto is_saturated(const measurement &result) noexcept -> bool
{
    return result.watermark_bytes == painted_bytes || result.first_call_bytes == painted_bytes
           || result.construction_bytes == painted_bytes;
}

auto watermark_thread(void *argument) noexcept -> void *
{
    auto *const result = static_cast<measurement *>(argument);

    // Construction and warm-up sit outside every painted window: construction is
    // offline and only the hot path is being measured.
    auto       chain = build_measured_chain();
    const auto call  = [&chain]() noexcept { run_measured_chain(chain); };
    call();

    result->solved             = chain_sink.load(std::memory_order_relaxed) != 1;
    result->first_call_bytes   = 0;
    result->construction_bytes = 0;
    result->floor_bytes        = painted_pass(call, false);
    result->watermark_bytes    = painted_pass(call, true);

#if WATERMARK_ROW == WATERMARK_ROW_NMPC_STATIC
    // ON THIS ROW THE FIRST CALL IS A DIFFERENT CHAIN FROM THE STEADY-STATE
    // ONE, AND IT IS THE DEEPER OF THE TWO. The solver instance is emplaced
    // lazily on the first solve and only reset on every later one, so the
    // first solve reaches a setup path that no steady-state solve enters. The
    // per-function frame report already shows that path's frame exceeding the
    // steady-state whole-chain peak at several dimensions, which would make a
    // supported maximum derived from the steady-state figure alone an
    // UNDERSTATEMENT -- and an understated stack figure overflows a task
    // rather than returning a wrong answer.
    //
    // It is measured on a freshly built chain, so construction stays outside
    // the painted window exactly as it does above and the difference between
    // the two figures is the first solve and nothing else.
    auto       fresh      = build_first_call_chain();
    const auto first_call = [&fresh]() noexcept { run_measured_chain(fresh); };
    result->first_call_bytes = painted_pass(first_call, true);

    const auto construct      = []() noexcept { run_construction(); };
    result->construction_bytes = painted_pass(construct, true);
#endif

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

    measurement result{0, 0, 0, 0, false};
    pthread_t   thread{};
    if(pthread_create(&thread, &attributes, &watermark_thread, &result) != 0)
    {
        std::fprintf(stderr, "stack_watermark: the measurement thread could not be created\n");
        static_cast<void>(pthread_attr_destroy(&attributes));
        return 1;
    }

    static_cast<void>(pthread_join(thread, nullptr));
    static_cast<void>(pthread_attr_destroy(&attributes));

    std::printf("row=%s nx=%zu ny=%zu nu=%zu",
                row_name,
                instantiated_state_dimension,
                instantiated_output_dimension,
                instantiated_input_dimension);

    // The horizon and the two dimensions it induces are printed by the one row
    // that has them and by no other. Printing a zero for a row that carries no
    // horizon would put a number where there is no quantity.
#if WATERMARK_ROW == WATERMARK_ROW_NMPC_STATIC
    std::printf(" nh=%zu nv=%zu maxm=%zu first_call_bytes=%zu construction_bytes=%zu",
                selected_horizon,
                instantiated_decision_dimension,
                instantiated_constraint_dimension,
                result.first_call_bytes,
                result.construction_bytes);
#endif

    std::printf(" gap_bytes=%zu painted_bytes=%zu saturated=%s scalar=double solved=%s"
                " watermark_bytes=%zu floor_bytes=%zu\n",
                harness_gap_bytes,
                painted_bytes,
                is_saturated(result) ? "yes" : "no",
                result.solved ? "yes" : "no",
                result.watermark_bytes,
                result.floor_bytes);

    return 0;
}
