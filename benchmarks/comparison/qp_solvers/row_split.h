#ifndef HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ROW_SPLIT_H
#define HPP_GUARD_BENCHMARKS_COMPARISON_QP_SOLVERS_ROW_SPLIT_H

#include <Eigen/Dense>

#include <vector>

namespace ctrlpp::bench::qp_solvers
{

/// Separates the rows with coinciding bounds (equalities) from the two-sided rows, for solvers that
/// pose the two kinds through different arguments, and puts their multipliers back in row order.
struct row_split
{
    std::vector<Eigen::Index> equalities;
    std::vector<Eigen::Index> inequalities;

    auto equality_rows(const Eigen::MatrixXd& A) const -> Eigen::MatrixXd
    {
        return A(equalities, Eigen::all);
    }

    auto inequality_rows(const Eigen::MatrixXd& A) const -> Eigen::MatrixXd
    {
        return A(inequalities, Eigen::all);
    }

    void gather(const Eigen::VectorXd& l, const Eigen::VectorXd& u, Eigen::VectorXd& b, Eigen::VectorXd& lower,
                Eigen::VectorXd& upper) const
    {
        b = l(equalities);
        lower = l(inequalities);
        upper = u(inequalities);
    }

    void scatter(const Eigen::VectorXd& y, const Eigen::VectorXd& z, Eigen::VectorXd& out) const
    {
        out.resize(static_cast<Eigen::Index>(equalities.size() + inequalities.size()));
        out(equalities) = y;
        out(inequalities) = z;
    }
};

inline auto split_rows(const Eigen::VectorXd& l, const Eigen::VectorXd& u) -> row_split
{
    row_split split;
    for(Eigen::Index i = 0; i < l.size(); ++i)
        (l(i) == u(i) ? split.equalities : split.inequalities).push_back(i);
    return split;
}

}

#endif
