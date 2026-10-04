#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_JUDGE_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_JUDGE_H

#include "arm.h"
#include "registry.h"

#include "qp/qp_accuracy.h"

#include <Eigen/Dense>

#include <cmath>
#include <limits>
#include <string_view>

namespace ctrlpp::bench::qp_solvers
{

/// The three parts of qp_accuracy.h's kkt_relative_residual and their maximum, which is that function's
/// value. An answer of the wrong shape or with a non-finite entry scores NaN throughout.
struct kkt_terms
{
    double stationarity;
    double primal;
    double complementarity;
    double worst;
};

inline auto judge(const problems::qp::dense_program& program, const Eigen::VectorXd& x, const Eigen::VectorXd& y)
    -> kkt_terms
{
    constexpr double missing = std::numeric_limits<double>::quiet_NaN();
    if(x.size() != program.P.cols() || y.size() != program.A.rows() || !x.allFinite() || !y.allFinite())
        return {missing, missing, missing, missing};
    kkt_terms terms{problems::qp::stationarity_residual(program, x, y), problems::qp::primal_violation(program, x),
                    problems::qp::complementarity_residual(program, x, y), 0.0};
    const bool finite = std::isfinite(terms.stationarity) && std::isfinite(terms.primal) && std::isfinite(terms.complementarity);
    terms.worst = finite ? std::max({terms.stationarity, terms.primal, terms.complementarity}) : missing;
    return terms;
}

/// NaN compares false, so a non-finite answer fails.
inline bool passes(const kkt_terms& terms, tier level)
{
    return terms.worst <= threshold(level);
}

inline auto fail_class(bool passed, solver_claim claim) -> std::string_view
{
    if(passed)
        return "none";
    switch(claim)
    {
        case solver_claim::solved:
            return "inaccurate";
        case solver_claim::iteration_cap:
            return "max_iter";
        case solver_claim::infeasible:
            return "infeasible_claim";
        default:
            return "error";
    }
}

}

#endif
