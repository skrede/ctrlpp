#ifndef HPP_GUARD_CTRLPP_TESTS_STUB_QP_SOLVER_H
#define HPP_GUARD_CTRLPP_TESTS_STUB_QP_SOLVER_H

// Backend-free QP solver satisfying the ctrlpp::qp_solver concept.
//
// Why this exists: every controller test translation unit that links a real QP
// backend lives in the exceptions carve-out tree, so a default-tree-only run
// exercises no controller construction at all. This stub has no external
// dependency, so the construction-time guards can be covered in the DEFAULT
// (-fno-exceptions) tree, which is where the shipping discipline is dogfooded.
//
// Usage contract:
//   1. Construct the stub (optionally with a report_lengths knob) and hand it
//      to a controller factory's solver-taking overload. The stub is moved in.
//   2. setup() records the problem dimensions, so a later solve() can size its
//      result from the problem rather than from a hard-coded number.
//   3. solve() reports solve_status::optimal and returns a zero primal of
//      `primal_length()` entries and a zero dual of `dual_length()` entries.
//
// The reported result lengths are a documented knob, defaulting to conforming:
//   * report_lengths::conforming  - primal is exactly the decision dimension of
//                                   the problem passed to setup(), dual is
//                                   exactly its constraint-row count. This is
//                                   the default and models a well-behaved
//                                   backend.
//   * report_lengths::short_primal - primal is one entry shorter than the
//                                   decision dimension; dual conforms.
//   * report_lengths::short_dual  - dual is one entry shorter than the
//                                   constraint-row count; primal conforms.
//   * report_lengths::empty       - both are empty, the degenerate case of a
//                                   backend that reported success and returned
//                                   nothing.
// The short and empty variants exist so result-shape validation can be tested
// against a backend that satisfies the concept, reports optimality, and still
// returns storage the extraction cannot legally read. Lengths are clamped at
// zero, so the short variants degrade to empty on a zero-sized problem.

#include "ctrlpp/mpc/qp_types.h"

#include <Eigen/Core>

#include <algorithm>

namespace ctrlpp_test
{

enum class report_lengths
{
    conforming,
    short_primal,
    short_dual,
    empty
};

template <typename Scalar>
struct stub_qp_solver
{
    using scalar_type = Scalar;

    report_lengths lengths{report_lengths::conforming};
    int n_dec{0};
    int n_con{0};
    int solve_count{0};

    stub_qp_solver() = default;

    explicit stub_qp_solver(report_lengths reported)
        : lengths{reported}
    {
    }

    void setup(const ctrlpp::qp_problem<Scalar>& problem)
    {
        n_dec = static_cast<int>(problem.P.cols());
        n_con = static_cast<int>(problem.A.rows());
    }

    [[nodiscard]] auto primal_length() const -> int
    {
        switch(lengths)
        {
        case report_lengths::short_primal:
            return std::max(n_dec - 1, 0);
        case report_lengths::empty:
            return 0;
        case report_lengths::conforming:
        case report_lengths::short_dual:
        default:
            return n_dec;
        }
    }

    [[nodiscard]] auto dual_length() const -> int
    {
        switch(lengths)
        {
        case report_lengths::short_dual:
            return std::max(n_con - 1, 0);
        case report_lengths::empty:
            return 0;
        case report_lengths::conforming:
        case report_lengths::short_primal:
        default:
            return n_con;
        }
    }

    auto solve(const ctrlpp::qp_update<Scalar>&) -> ctrlpp::qp_result<Scalar>
    {
        ++solve_count;
        ctrlpp::qp_result<Scalar> result;
        result.status = ctrlpp::solve_status::optimal;
        result.x = Eigen::VectorX<Scalar>::Zero(primal_length());
        result.y = Eigen::VectorX<Scalar>::Zero(dual_length());
        result.objective = Scalar{0};
        result.solve_time = Scalar{0};
        result.iterations = 1;
        result.primal_residual = Scalar{0};
        result.dual_residual = Scalar{0};
        return result;
    }
};

}

#endif
