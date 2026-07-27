#ifndef HPP_GUARD_CTRLPP_TESTS_STUB_NLP_SOLVER_H
#define HPP_GUARD_CTRLPP_TESTS_STUB_NLP_SOLVER_H

// Backend-free NLP solver satisfying the ctrlpp::nlp_solver concept, the
// nonlinear twin of stub_qp_solver.h. It exists for the same reason: every
// nonlinear controller test translation unit that links a real NLP backend
// lives in the exceptions carve-out tree, so a default-tree-only run would
// exercise no nonlinear controller construction at all.
//
// Usage contract:
//   1. Construct the stub (optionally with a report_lengths knob) and hand it
//      to a controller factory's solver-taking overload. The stub is moved in.
//   2. setup() records the decision dimension of the problem, so a later
//      solve() can size its result from the problem rather than from a
//      hard-coded number.
//   3. solve() reports solve_status::optimal and returns a zero decision vector
//      of `primal_length()` entries.
//
// The reported result length reuses ctrlpp_test::report_lengths from
// stub_qp_solver.h and defaults to conforming. `short_dual` has no distinct
// meaning here (an NLP result carries no dual block), so it behaves as
// conforming. Every variant still reports solve_status::optimal, so a consumer
// that trusts the reported status alone is handed a decision vector it cannot
// legally read.
//
// As on the QP stub, the knob is settable as a template argument as well as a
// runtime field, because the nonlinear moving-horizon estimator
// default-constructs its solver member and offers no injection seam.

#include "ctrlpp/mpc/nlp_solver.h"

#include "stub_qp_solver.h"

#include <Eigen/Core>

#include <algorithm>

namespace ctrlpp_test
{

template <typename Scalar, report_lengths Reported = report_lengths::conforming>
struct stub_nlp_solver
{
    using scalar_type = Scalar;

    report_lengths lengths{Reported};
    int n_vars{0};
    int solve_count{0};

    stub_nlp_solver() = default;

    explicit stub_nlp_solver(report_lengths reported)
        : lengths{reported}
    {
    }

    void setup(const ctrlpp::nlp_problem<Scalar>& problem) { n_vars = problem.n_vars; }

    // The compile-time-dimension contract is a distinct type, not a conversion
    // of the runtime-erased one, so the static controller path needs its own
    // overload here.
    template <int NV>
    void setup(const ctrlpp::nlp_problem_static<Scalar, NV>& problem)
    {
        n_vars = problem.n_vars;
    }

    [[nodiscard]] auto primal_length() const -> int
    {
        switch(lengths)
        {
        case report_lengths::short_primal:
            return std::max(n_vars - 1, 0);
        case report_lengths::empty:
            return 0;
        case report_lengths::conforming:
        case report_lengths::short_dual:
        default:
            return n_vars;
        }
    }

    auto solve(const ctrlpp::nlp_update<Scalar>&) -> ctrlpp::nlp_result<Scalar>
    {
        ++solve_count;
        ctrlpp::nlp_result<Scalar> result{};
        result.status = ctrlpp::solve_status::optimal;
        result.x = Eigen::VectorX<Scalar>::Zero(primal_length());
        result.objective = Scalar{0};
        result.solve_time = Scalar{0};
        result.iterations = 1;
        result.primal_residual = Scalar{0};
        return result;
    }
};

}

#endif
