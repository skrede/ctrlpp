#include "dual_recovery.h"

#include <Eigen/SparseQR>
#include <Eigen/OrderingMethods>

#include <cmath>
#include <vector>

namespace ctrlpp::bench::qp_solvers
{

namespace
{

auto active_rows(const problems::qp::dense_program& program, const Eigen::VectorXd& row, double threshold)
    -> std::vector<Eigen::Index>
{
    std::vector<Eigen::Index> active;
    for(Eigen::Index i = 0; i < row.size(); ++i)
    {
        const double l = program.l(i);
        const double u = program.u(i);
        const bool equality = l == u;
        const bool at_lower = row(i) - l <= threshold * (1.0 + std::abs(l));
        const bool at_upper = u - row(i) <= threshold * (1.0 + std::abs(u));
        if(equality || at_lower || at_upper)
            active.push_back(i);
    }
    return active;
}

auto active_columns(const Eigen::SparseMatrix<double>& A, const std::vector<Eigen::Index>& active)
    -> Eigen::SparseMatrix<double>
{
    const Eigen::SparseMatrix<double> At = A.transpose();
    Eigen::SparseMatrix<double> columns(At.rows(), static_cast<Eigen::Index>(active.size()));
    std::vector<Eigen::Triplet<double>> entries;
    for(std::size_t j = 0; j < active.size(); ++j)
        for(Eigen::SparseMatrix<double>::InnerIterator it(At, active[j]); it; ++it)
            entries.emplace_back(it.row(), static_cast<Eigen::Index>(j), it.value());
    columns.setFromTriplets(entries.begin(), entries.end());
    columns.makeCompressed();
    return columns;
}

}

auto recover_duals(const problems::qp::dense_program& program, const Eigen::SparseMatrix<double>& A,
                   const Eigen::VectorXd& x, double threshold) -> Eigen::VectorXd
{
    Eigen::VectorXd y = Eigen::VectorXd::Zero(program.A.rows());
    if(!x.allFinite() || x.size() != program.P.cols())
        return Eigen::VectorXd::Constant(program.A.rows(), std::nan(""));
    const std::vector<Eigen::Index> active = active_rows(program, A * x, threshold);
    if(active.empty())
        return y;
    const Eigen::SparseMatrix<double> columns = active_columns(A, active);
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr(columns);
    const Eigen::VectorXd gradient = program.P * x + program.q;
    const Eigen::VectorXd solved = qr.solve(-gradient);
    for(std::size_t j = 0; j < active.size(); ++j)
        y(active[j]) = solved(static_cast<Eigen::Index>(j));
    return y;
}

}
