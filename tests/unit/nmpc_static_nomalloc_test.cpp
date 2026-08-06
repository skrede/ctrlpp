// Prove the compile-time-horizon nmpc_static solve hot path does zero heap
// allocation in steady state, using the belt-and-suspenders harness from
// nomalloc_harness.h: an eigen_assert that stores into a pollable sentinel and
// survives -DNDEBUG, plus a global operator-new counter that catches heap
// traffic outside Eigen's own bookkeeping. The harness header must stay the
// FIRST include of this file.
//
// The first case is the shipped pin: with the decision dimension NV fixed at
// compile time, the argmin bridge sizes its decision-vector storage with
// fixed-size Eigen types and argmin's fixed-N NW-SQP substrate solves without
// heap growth, so after a warm-up solve the steady-state step/solve loop
// allocates nothing.
//
// The second case walks the boundary below, over the same corpus the published
// stack figures for this row were taken on.
//
// THIS ROW ARMS TWO CONFIGURATIONS BESIDE THE PIN, NOT THE 25 IT PUBLISHES A
// FIGURE FOR, AND THE REASON IS MEASURED RATHER THAN ASSERTED. Every point here
// instantiates the whole nonlinear-programming substrate and then solves at that
// dimension. A unit holding all 25 peaked at 21.9 GB of compiler resident set
// and was killed by the kernel before it finished; a unit holding eight -- the
// endpoint of every published line plus this bracket -- built, ran green at all
// eight, and cost 17.4 GB to compile and 232 seconds to run, against roughly one
// second for the pin alone. Neither is a cost this suite can carry on an
// ordinary machine. What is kept is the pair that WALKS THE BOUNDARY, which is
// the fact this row did not previously have; the published endpoints above it
// are NOT armed and nothing here is a claim about them.
//
// THE HEAP CLAIM HAS A BOUNDARY ON THIS ROW TOO, AND IT IS AT A DECISION
// DIMENSION OF 48. Eigen applies a Householder sequence BY BLOCK once the
// sequence is at least BlockSize = 48 long and the destination has more than one
// column (HouseholderSequence.h), and the blocked application declares its block
// reflector factor with the destination's compile-time column count, which for a
// dynamically sized destination block is Dynamic -- so the factor is HEAP
// ALLOCATED, through Eigen's own aligned allocator. The five estimator rows
// cross the same boundary at a measurement dimension of 48, walked from both
// sides in estimation_<row>_nomalloc_test.cpp; here it is walked in the decision
// dimension, which is the dimension the published table for this row is indexed
// by. Every point below states which side it is on, and the two bracket points
// exist only to put the boundary between two measured configurations rather than
// between a measured one and an assertion.
//
// THE GLOBAL COUNTER CANNOT SEE THAT ALLOCATION -- it goes through Eigen's
// aligned_malloc, which calls std::malloc directly -- so the allocating points
// assert both that the sentinel fires AND that the counter stayed at zero, which
// is the measured demonstration that a guard carrying only the counter would
// report every one of them allocation-free.

#include "nomalloc_harness.h"

#include "ctrlpp/types.h"

#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/nmpc.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cstddef>
#include <utility>


namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NH = 5;
constexpr double dt = 0.1;

// NV = (NH+1)*NX + NH*NU = 6*2 + 5*1 = 17 for this slack-free double_integrator
// cell; the caller pins the static solver on the same compile-time dimension.
constexpr int NV = static_cast<int>((NH + 1) * NX + NH * NU);

// MaxM caps the CONSTRAINT axis: this cell carries only the
// dynamics-continuity + initial-state equalities, M_eq = NX*(NH+1) = 2*6 = 12,
// and no general inequalities (box bounds, if any, are free via argmin's +2N
// slack). Binding MaxM = M_eq makes argmin's per-call result-multiplier storage
// inline, closing the last per-solve heap on argmin's side.
constexpr int MaxM = static_cast<int>(NX * (NH + 1));

auto double_integrator = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{ return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)}; };

using StaticSolver = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, NV, MaxM>;
using NmpcStaticDI = ctrlpp::nmpc_static<double, NX, NU, NH, StaticSolver, decltype(double_integrator)>;


// The decision dimension at which the linear-algebra library switches to its
// blocked Householder application. Read out of that library rather than fitted,
// and walked from both sides by the two bracket points in the grid below.
constexpr int blocked_householder_boundary = 48;

// The corpus, in the shape the stack-watermark instrument builds it, so an armed
// point and a published stack figure are the same configuration of the same
// plant. Forward-Euler damped chain, one input per group of states.
template <std::size_t N>
auto chain_state_matrix() -> ctrlpp::Matrix<double, N, N>
{
    constexpr double step = 0.01;

    ctrlpp::Matrix<double, N, N> A = ctrlpp::Matrix<double, N, N>::Identity();
    for(std::size_t i = 0; i < N; ++i)
        A(int(i), int(i)) += step * -0.5;
    for(std::size_t i = 0; i + 1 < N; ++i)
        A(int(i), int(i + 1)) = step * 1.0;
    return A;
}

template <std::size_t N, std::size_t M>
auto chain_input_matrix() -> ctrlpp::Matrix<double, N, M>
{
    constexpr double step = 0.01;

    ctrlpp::Matrix<double, N, M> B = ctrlpp::Matrix<double, N, M>::Zero();
    const std::size_t            group = M <= N ? N / M : std::size_t{1};
    for(std::size_t j = 0; j < M; ++j)
    {
        const std::size_t reach = (j + 1) * group;
        B(int((reach < N ? reach : N) - 1), int(j)) = step;
    }
    return B;
}

template <std::size_t N, std::size_t M>
struct chain_dynamics
{
    ctrlpp::Matrix<double, N, N> F;
    ctrlpp::Matrix<double, N, M> G;

    auto operator()(const ctrlpp::Vector<double, N>& x, const ctrlpp::Vector<double, M>& u) const
        -> ctrlpp::Vector<double, N>
    {
        return ctrlpp::Vector<double, N>{F * x + G * u};
    }
};

// One armed grid point: the three dimensions a caller of this row chooses. The
// two the frames are sized by are DERIVED from them and are computed here rather
// than restated.
template <std::size_t States, std::size_t Inputs, std::size_t Horizon>
struct configuration
{
    static constexpr std::size_t states  = States;
    static constexpr std::size_t inputs  = Inputs;
    static constexpr std::size_t horizon = Horizon;

    static constexpr int decision_dimension   = static_cast<int>((Horizon + 1) * States + Horizon * Inputs);
    static constexpr int constraint_dimension = static_cast<int>(States * (Horizon + 1));
};

// The traversal the estimator grids use, over this row's three-dimensional
// points: arm() applies the callable to every configuration and returns how many
// it armed, so the count a case asserts on is produced by the traversal rather
// than restated beside it.
template <typename... Configurations>
struct over
{
    static constexpr std::size_t size = sizeof...(Configurations);

    template <typename Arm>
    static std::size_t arm(Arm&& arm_point)
    {
        (arm_point(Configurations{}), ...);
        return size;
    }
};

// Warm-up-then-arm for one configuration. Construction and the walk into steady
// state happen OUTSIDE the window; the count is sampled inside it and before any
// Catch2 macro runs, because framework macros may themselves allocate. Which
// side of the boundary the point is on decides what is asserted, and the counter
// is required to read zero on BOTH sides -- above the boundary because that is
// the measured demonstration that it cannot see the library's own allocator.
template <typename Configuration>
void arm_configuration(Configuration)
{
    constexpr std::size_t states  = Configuration::states;
    constexpr std::size_t inputs  = Configuration::inputs;
    constexpr std::size_t horizon = Configuration::horizon;
    constexpr int         nv      = Configuration::decision_dimension;
    constexpr int         maxm    = Configuration::constraint_dimension;

    using solver_type     = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, nv, maxm>;
    using controller_type = ctrlpp::nmpc_static<double, states, inputs, horizon, solver_type,
                                                chain_dynamics<states, inputs>>;

    static_assert(controller_type::problem_dimension == nv,
        "nmpc_static must derive the compile-time decision dimension NV");
    static_assert(solver_type::bridge_type::strict_allocation_free,
        "the fixed-(NV,MaxM) argmin bridge backing nmpc_static must be strict_allocation_free");

    const chain_dynamics<states, inputs> dynamics{chain_state_matrix<states>(),
                                                  chain_input_matrix<states, inputs>()};

    ctrlpp::nmpc_config<double, states, inputs> config;
    config.horizon = static_cast<int>(horizon);
    config.Q       = ctrlpp::Matrix<double, states, states>::Identity() * 10.0;
    config.R       = ctrlpp::Matrix<double, inputs, inputs>::Identity() * 0.1;

    controller_type controller{dynamics, config};

    ctrlpp::Vector<double, states> x = ctrlpp::Vector<double, states>::Zero();
    x(0)                             = 1.0;

    bool warm_ok = true;
    for(int step = 0; step < 20; ++step)
    {
        auto warm = controller.solve(x);
        warm_ok   = warm_ok && warm.has_value();
        if(warm)
            x = dynamics(x, warm->input);
    }
    REQUIRE(warm_ok);

    std::size_t allocations = 0;
    bool        all_ok      = true;
    bool        eigen_fired = false;
    {
        ctrlpp_test::scoped_no_malloc guard;
        for(int step = 0; step < 4; ++step)
        {
            auto u = controller.solve(x);
            all_ok = all_ok && u.has_value();
            if(u)
                x = dynamics(x, u->input);
        }
        allocations = guard.allocations();
        eigen_fired = guard.eigen_violation();
    }

    REQUIRE(all_ok);
    if constexpr(nv >= blocked_householder_boundary)
        REQUIRE(eigen_fired);
    else
        REQUIRE_FALSE(eigen_fired);

    REQUIRE(allocations == 0);
    REQUIRE(x.allFinite());
}

// The armed set: both sides of the boundary, and nothing else. Two
// configurations chosen so that the decision dimension steps across 48 while
// staying small enough that the pair costs a fraction of what the published
// endpoints cost to compile and to run.
constexpr std::size_t armed_points = 2;

using controller_grid = over<configuration<3, 1, 11>, configuration<6, 1, 6>>;

static_assert(controller_grid::size == armed_points,
              "the armed count must be produced by the list rather than restated beside it");

// The bracket, spelled out so the two configurations that carry it cannot be
// dropped as duplicates of a ladder rung.
static_assert(configuration<3, 1, 11>::decision_dimension == blocked_householder_boundary - 1);
static_assert(configuration<6, 1, 6>::decision_dimension == blocked_householder_boundary);

}


TEST_CASE("nmpc_static steady-state solve performs zero heap allocation",
          "[nmpc][argmin][static][hardening][nomalloc]")
{
    // The static bridge must carry the compile-time decision dimension and be
    // marked allocation-free by type before we ever measure the runtime path.
    static_assert(NmpcStaticDI::problem_dimension == NV,
        "nmpc_static must derive the compile-time decision dimension NV");
    static_assert(StaticSolver::bridge_type::strict_allocation_free,
        "the fixed-(NV,MaxM) argmin bridge backing nmpc_static must be "
        "strict_allocation_free (both decision and constraint axes bound)");

    auto config = ctrlpp::nmpc_config<double, NX, NU>{
        .horizon = static_cast<int>(NH),
        .Q = 10.0 * Eigen::Matrix2d::Identity(),
        .R = 0.1 * Eigen::Matrix<double, 1, 1>::Identity(),
    };

    NmpcStaticDI controller{double_integrator, config};

    Eigen::Vector2d x{1.0, 0.0};

    // Warm-up OUTSIDE the armed window: the first solves construct the argmin
    // solver state and flush any one-time lazy instantiation.
    for(int step = 0; step < 20; ++step)
    {
        auto warm = controller.solve(x);
        REQUIRE(warm.has_value());
        x = double_integrator(x, warm->input);
    }

    // Steady-state window: arm the guard, run a step/solve loop, and sample the
    // allocation count immediately, before any Catch2 macro runs inside the
    // window (framework macros may themselves allocate).
    std::size_t allocations = 0;
    bool all_ok = true;
    bool eigen_ok = true;
    {
        ctrlpp_test::scoped_no_malloc guard;
        for(int step = 0; step < 32; ++step)
        {
            auto u = controller.solve(x);
            all_ok = all_ok && u.has_value();
            x = double_integrator(x, u->input);
        }
        allocations = guard.allocations();
        eigen_ok = !guard.eigen_violation();
    }

    REQUIRE(eigen_ok);
    REQUIRE(all_ok);
    REQUIRE(allocations == 0);
    REQUIRE(x.allFinite());
}

TEST_CASE("nmpc_static is allocation-free below a decision dimension of 48 and allocates in the solve "
          "at and above it",
          "[nmpc][argmin][static][hardening][nomalloc][grid]")
{
    const std::size_t armed = controller_grid::arm([](auto point) { arm_configuration(point); });

    REQUIRE(armed == armed_points);
}
