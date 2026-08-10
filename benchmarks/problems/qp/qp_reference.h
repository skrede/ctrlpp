#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_QP_QP_REFERENCE_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_QP_QP_REFERENCE_H

#include "qp/qp_accuracy.h"

#include <Eigen/Dense>

#include <vector>
#include <cstdint>

namespace ctrlpp::bench::problems::qp
{

enum class row_state : std::int8_t
{
    at_lower = -1,
    inactive = 0,
    at_upper = 1
};

struct reference_answer
{
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    std::vector<row_state> rows;
};

inline auto row_states_of(std::int32_t code, Eigen::Index count) -> std::vector<row_state>
{
    std::vector<row_state> states(static_cast<std::size_t>(count), row_state::inactive);
    for(Eigen::Index i = 0; i < count; ++i)
    {
        states[static_cast<std::size_t>(i)] = static_cast<row_state>(code % 3 - 1);
        code /= 3;
    }
    return states;
}

inline auto active_rows_of(const std::vector<row_state>& states) -> std::vector<Eigen::Index>
{
    std::vector<Eigen::Index> active;
    for(std::size_t i = 0; i < states.size(); ++i)
        if(states[i] != row_state::inactive)
            active.push_back(static_cast<Eigen::Index>(i));
    return active;
}

inline auto bound_held_by(const dense_program& program, const std::vector<row_state>& states, Eigen::Index row)
    -> double
{
    return states[static_cast<std::size_t>(row)] == row_state::at_upper ? program.u(row) : program.l(row);
}

/// Stationarity of the Lagrangian, together with equality on the rows the
/// assignment holds at a bound. Every other multiplier is zero by
/// complementarity, so the assignment fixes the system completely.
inline auto answer_for(const dense_program& program, const std::vector<row_state>& states) -> reference_answer
{
    const std::vector<Eigen::Index> active = active_rows_of(states);
    const Eigen::Index n = program.P.cols();
    const Eigen::Index m = static_cast<Eigen::Index>(active.size());
    Eigen::MatrixXd system = Eigen::MatrixXd::Zero(n + m, n + m);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + m);
    system.topLeftCorner(n, n) = program.P;
    rhs.head(n) = -program.q;
    for(Eigen::Index j = 0; j < m; ++j)
    {
        system.row(n + j).head(n) = program.A.row(active[j]);
        system.col(n + j).head(n) = program.A.row(active[j]).transpose();
        rhs(n + j) = bound_held_by(program, states, active[j]);
    }
    const Eigen::VectorXd solution = system.fullPivLu().solve(rhs);
    Eigen::VectorXd y = Eigen::VectorXd::Zero(program.A.rows());
    for(Eigen::Index j = 0; j < m; ++j)
        y(active[j]) = solution(n + j);
    return {solution.head(n), y, states};
}

inline bool is_optimal(const dense_program& program, const reference_answer& answer)
{
    const Eigen::VectorXd row = program.A * answer.x;
    for(Eigen::Index i = 0; i < row.size(); ++i)
    {
        const row_state state = answer.rows[static_cast<std::size_t>(i)];
        if(state == row_state::inactive && (row(i) > program.u(i) || row(i) < program.l(i)))
            return false;
        if(state == row_state::at_upper && answer.y(i) < 0.0)
            return false;
        if(state == row_state::at_lower && answer.y(i) > 0.0)
            return false;
    }
    return true;
}

/// A strictly convex program has one Karush-Kuhn-Tucker point, so the single
/// assignment of the rows to {inside its bounds, at its lower bound, at its
/// upper bound} that satisfies the conditions names the active set exactly, with
/// no tolerance entering anywhere. The cost is 3^rows small solves, which is why
/// this serves a corpus of a handful of rows and is not a solver.
inline auto solve_by_enumeration(const dense_program& program) -> reference_answer
{
    std::int32_t assignments = 1;
    for(Eigen::Index i = 0; i < program.A.rows(); ++i)
        assignments *= 3;
    for(std::int32_t code = 0; code < assignments; ++code)
    {
        const reference_answer candidate = answer_for(program, row_states_of(code, program.A.rows()));
        if(is_optimal(program, candidate))
            return candidate;
    }
    return {};
}

inline auto active_row_count(const reference_answer& answer) -> std::int32_t
{
    return static_cast<std::int32_t>(active_rows_of(answer.rows).size());
}

/// A comparison whose rows are all slack at the solution races the solvers on an
/// effectively unconstrained program, whatever the solvers are.
inline bool poses_active_constraint(const dense_program& program)
{
    return active_row_count(solve_by_enumeration(program)) > 0;
}

}

#endif
