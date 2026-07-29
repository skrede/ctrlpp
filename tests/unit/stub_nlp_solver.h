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
// The reported result length and value controls reuse their definitions from
// stub_qp_solver.h. `short_dual` has no distinct meaning here (an NLP result
// carries no dual block), so it behaves as conforming. The status is independently
// configurable, so consumers can be tested against malformed results under
// every status they accept.
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

template <typename Scalar,
          report_lengths Reported = report_lengths::conforming,
          report_values Values = report_values::finite,
          ctrlpp::solve_status Status = ctrlpp::solve_status::optimal>
struct stub_nlp_solver
{
    using scalar_type = Scalar;

    report_lengths lengths{Reported};
    report_values values{Values};
    ctrlpp::solve_status status{Status};
    int n_vars{0};
    int solve_count{0};

    stub_nlp_solver() = default;

    explicit stub_nlp_solver(report_lengths reported)
        : lengths{reported}
    {
    }

    stub_nlp_solver(report_lengths reported_lengths,
                    report_values reported_values,
                    ctrlpp::solve_status reported_status)
        : lengths{reported_lengths}
        , values{reported_values}
        , status{reported_status}
    {
    }

    auto setup(const ctrlpp::nlp_problem<Scalar>& problem) -> ctrlpp::expected<void, stub_setup_error>
    {
        n_vars = problem.n_vars;
        return {};
    }

    // The compile-time-dimension contract is a distinct type, not a conversion
    // of the runtime-erased one, so the static controller path needs its own
    // overload here.
    template <int NV>
    auto setup(const ctrlpp::nlp_problem_static<Scalar, NV>& problem) -> ctrlpp::expected<void, stub_setup_error>
    {
        n_vars = problem.n_vars;
        return {};
    }

    auto primal_length() const -> int
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
        result.status = status;
        result.x = Eigen::VectorX<Scalar>::Zero(primal_length());
        result.objective = Scalar{0};
        result.solve_time = Scalar{0};
        result.iterations = 1;
        result.primal_residual = Scalar{0};
        if(values != report_values::finite)
        {
            Scalar bad{};
            switch(values)
            {
            case report_values::nan:
                bad = std::numeric_limits<Scalar>::quiet_NaN();
                break;
            case report_values::positive_infinity:
                bad = std::numeric_limits<Scalar>::infinity();
                break;
            case report_values::negative_infinity:
                bad = -std::numeric_limits<Scalar>::infinity();
                break;
            case report_values::finite:
                break;
            }
            result.x.setConstant(bad);
            result.objective = bad;
            result.solve_time = bad;
            result.primal_residual = bad;
        }
        return result;
    }
};

}

#endif
