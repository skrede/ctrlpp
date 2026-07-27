// Regression coverage for the argmin QP backend policy.
//
// argmin_qp_solver.h is __has_include-gated on argmin/qp/sparse_admm_qp.h. The
// default argmin pin (the milestone/v0.3.5 tip) ships that header, so with
// CTRLPP_HAS_ARGMIN the policy is defined by default and the real checks below run:
// concept conformance and a functional solve on a small bound-constrained QP with
// a closed-form optimum. Only an override to a pre-argmin/qp/ pin leaves the policy
// compiled out, in which case this file is an empty (trivially passing) TU.

#include "ctrlpp/mpc/argmin_qp_solver.h"

#include <catch2/catch_test_macros.hpp>

// Always-present sentinel: if the pin is overridden back to a pre-argmin/qp/ SHA
// the QP header is absent and every gated case below compiles out, which would
// leave a Catch2 binary with zero registered tests -- and Catch2 exits non-zero
// ("No tests ran"), reddening CI. This keeps at least one test case in every
// configuration and records which path was taken.
TEST_CASE("argmin_qp_solver policy gate", "[mpc][argmin][qp]")
{
#if defined(CTRLPP_HAS_ARGMIN_QP)
    SUCCEED("argmin QP header present: policy active, functional cases below run");
#else
    SUCCEED("argmin QP header absent (pin override to a pre-QP SHA): policy compiled out");
#endif
}

#if defined(CTRLPP_HAS_ARGMIN_QP)

#include "ctrlpp/mpc/qp_solver.h"
#include "ctrlpp/mpc/qp_types.h"

#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <Eigen/Dense>
#include <Eigen/Sparse>

using Catch::Matchers::WithinAbs;

static_assert(ctrlpp::qp_solver<ctrlpp::argmin_qp_solver>,
              "argmin_qp_solver must model the qp_solver concept");

namespace
{

// Identity Hessian, zero linear term, box l <= x <= u. The unconstrained
// optimum is x = -q = 0, which lies inside the box, so x* = 0 with the box
// inactive -- an oracle the ADMM iterate must reach independently of scaling.
auto make_box_qp() -> ctrlpp::qp_problem<double>
{
    Eigen::SparseMatrix<double> P(2, 2);
    P.insert(0, 0) = 1.0;
    P.insert(1, 1) = 1.0;
    P.makeCompressed();

    Eigen::SparseMatrix<double> A(2, 2);
    A.insert(0, 0) = 1.0;
    A.insert(1, 1) = 1.0;
    A.makeCompressed();

    return {.P = P,
            .q = Eigen::Vector2d::Zero(),
            .A = A,
            .l = Eigen::Vector2d::Constant(-1.0),
            .u = Eigen::Vector2d::Constant(1.0)};
}

} // namespace

TEST_CASE("argmin_qp_solver setup then resolve reaches the analytic optimum", "[mpc][argmin][qp]")
{
    auto problem = make_box_qp();

    ctrlpp::argmin_qp_solver solver;
    REQUIRE(solver.setup(problem).has_value());

    ctrlpp::qp_update<double> update{
        .q = problem.q,
        .l = problem.l,
        .u = problem.u,
        .warm_x = {},
        .warm_y = {}};

    auto result = solver.solve(update);

    CHECK(result.status == ctrlpp::solve_status::optimal);
    REQUIRE(result.x.size() == 2);
    CHECK_THAT(result.x(0), WithinAbs(0.0, 1e-4));
    CHECK_THAT(result.x(1), WithinAbs(0.0, 1e-4));
}

TEST_CASE("argmin_qp_solver preset constructors both solve the box QP", "[mpc][argmin][qp]")
{
    auto problem = make_box_qp();
    ctrlpp::qp_update<double> update{
        .q = problem.q,
        .l = problem.l,
        .u = problem.u,
        .warm_x = {},
        .warm_y = {}};

    for(auto preset : {ctrlpp::qp_preset::accuracy, ctrlpp::qp_preset::speed})
    {
        ctrlpp::argmin_qp_solver solver{preset};
        REQUIRE(solver.setup(problem).has_value());
        auto result = solver.solve(update);
        CHECK(result.status == ctrlpp::solve_status::optimal);
        // Both presets reach the optimum; accuracy is tighter, so check speed
        // only to the (looser) stopping tolerance.
        REQUIRE(result.x.size() == 2);
        CHECK_THAT(result.x(0), WithinAbs(0.0, 1e-2));
        CHECK_THAT(result.x(1), WithinAbs(0.0, 1e-2));
    }
}

TEST_CASE("argmin_qp_solver reports error before setup", "[mpc][argmin][qp]")
{
    ctrlpp::argmin_qp_solver solver;
    ctrlpp::qp_update<double> update{
        .q = Eigen::Vector2d::Zero(),
        .l = Eigen::Vector2d::Constant(-1.0),
        .u = Eigen::Vector2d::Constant(1.0),
        .warm_x = {},
        .warm_y = {}};

    // No setup() has run, so there is no factorization to resolve against.
    CHECK(solver.solve(update).status == ctrlpp::solve_status::error);
}

#endif // CTRLPP_HAS_ARGMIN_QP
