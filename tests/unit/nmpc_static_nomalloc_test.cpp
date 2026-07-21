// Prove the compile-time-horizon nmpc_static solve hot path does zero heap
// allocation in steady state, using the belt-and-suspenders harness from
// nomalloc_harness.h: a throwing eigen_assert that survives -DNDEBUG plus a
// global operator-new counter that catches heap traffic outside Eigen's own
// bookkeeping. The harness header must stay the FIRST include of this file.
//
// This is the EMB-02 proof for the static NMPC path (Route A / SEED-002): with
// the decision dimension NV pinned at compile time, the argmin bridge sizes its
// decision-vector storage with fixed-size Eigen types and argmin's fixed-N
// NW-SQP substrate solves without heap growth, so after a warm-up solve the
// steady-state step/solve loop allocates nothing. Strict-zero sign-off is a
// joint ctrlpp+argmin measurement; this test pins the ctrlpp-side contract at
// argmin's fixed-N floor.

#include "nomalloc_harness.h"

#include "ctrlpp/mpc/argmin_solver.h"
#include "ctrlpp/nmpc.h"

#include <catch2/catch_test_macros.hpp>

#include <Eigen/Dense>

#include <cstddef>


namespace
{

constexpr std::size_t NX = 2;
constexpr std::size_t NU = 1;
constexpr std::size_t NH = 5;
constexpr double dt = 0.1;

// NV = (NH+1)*NX + NH*NU = 6*2 + 5*1 = 17 for this slack-free double_integrator
// cell; the caller pins the static solver on the same compile-time dimension.
constexpr int NV = static_cast<int>((NH + 1) * NX + NH * NU);

// MaxM caps the CONSTRAINT axis (argmin SEED-044): this cell carries only the
// dynamics-continuity + initial-state equalities, M_eq = NX*(NH+1) = 2*6 = 12,
// and no general inequalities (box bounds, if any, are free via argmin's +2N
// slack). Binding MaxM = M_eq makes argmin's per-call result-multiplier storage
// inline, closing the last per-solve heap on argmin's side.
constexpr int MaxM = static_cast<int>(NX * (NH + 1));

auto double_integrator = [](const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{ return Eigen::Vector2d{x(0) + dt * x(1), x(1) + dt * u(0)}; };

using StaticSolver = ctrlpp::argmin_solver<double, ctrlpp::argmin_nw_sqp, true, NV, MaxM>;
using NmpcStaticDI = ctrlpp::nmpc_static<double, NX, NU, NH, StaticSolver, decltype(double_integrator)>;

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
