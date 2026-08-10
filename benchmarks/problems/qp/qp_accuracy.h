#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_QP_QP_ACCURACY_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_QP_QP_ACCURACY_H

#include "qp/dense_mpc_shaped.h"

#include <Eigen/Dense>

#include <algorithm>

namespace ctrlpp::bench::problems::qp
{

/// The dense form of the posed program, held once so every arm's answer is
/// scored against the same operands rather than against the shape its own
/// backend happened to receive them in.
struct dense_program
{
    Eigen::MatrixXd P;
    Eigen::VectorXd q;
    Eigen::MatrixXd A;
    Eigen::VectorXd l;
    Eigen::VectorXd u;
};

/// Each part of the optimality condition is divided by the sum of the norms of
/// the terms that cancel to form it, so the figure is dimensionless and no
/// tolerance is chosen anywhere. A part whose terms all vanish is reported
/// absolutely rather than divided by zero.
inline double relative_part(double residual, double scale)
{
    return scale > 0.0 ? residual / scale : residual;
}

inline double stationarity_residual(const dense_program& program, const Eigen::VectorXd& x,
                                    const Eigen::VectorXd& y)
{
    const Eigen::VectorXd curvature = program.P * x;
    const Eigen::VectorXd carried = program.A.transpose() * y;
    return relative_part((curvature + program.q + carried).lpNorm<Eigen::Infinity>(),
                         curvature.lpNorm<Eigen::Infinity>() + program.q.lpNorm<Eigen::Infinity>()
                             + carried.lpNorm<Eigen::Infinity>());
}

inline double primal_violation(const dense_program& program, const Eigen::VectorXd& x)
{
    const Eigen::VectorXd row = program.A * x;
    const double excess = (row - program.u).cwiseMax(0.0).lpNorm<Eigen::Infinity>();
    const double shortfall = (program.l - row).cwiseMax(0.0).lpNorm<Eigen::Infinity>();
    return relative_part(std::max(excess, shortfall),
                         row.lpNorm<Eigen::Infinity>() + program.l.lpNorm<Eigen::Infinity>()
                             + program.u.lpNorm<Eigen::Infinity>());
}

/// A multiplier is entitled to be nonzero only on the bound its own row sits on,
/// so the product of the multiplier with the slack against that bound must
/// vanish. The two signs are separated because a two-sided row's multiplier
/// names which of its bounds is active by its sign.
inline double complementarity_residual(const dense_program& program, const Eigen::VectorXd& x,
                                       const Eigen::VectorXd& y)
{
    const Eigen::VectorXd row = program.A * x;
    const double upper = y.cwiseMax(0.0).cwiseProduct(program.u - row).lpNorm<Eigen::Infinity>();
    const double lower = y.cwiseMin(0.0).cwiseProduct(program.l - row).lpNorm<Eigen::Infinity>();
    return relative_part(std::max(upper, lower),
                         y.lpNorm<Eigen::Infinity>()
                             * (row.lpNorm<Eigen::Infinity>() + program.l.lpNorm<Eigen::Infinity>()
                                + program.u.lpNorm<Eigen::Infinity>()));
}

/// The worst of the three conditions a primal-dual pair must satisfy to be the
/// program's answer, taken on that arm's OWN pair. It is the quadratic
/// program's counterpart of a Riccati residual: it says the answer solves the
/// problem, which a cross-arm deviation cannot say because two wrong answers can
/// agree.
inline double kkt_relative_residual(const dense_program& program, const Eigen::VectorXd& x,
                                    const Eigen::VectorXd& y)
{
    return std::max({stationarity_residual(program, x, y), primal_violation(program, x),
                     complementarity_residual(program, x, y)});
}

inline dense_program make_dense_mpc_program()
{
    return {Eigen::MatrixXd(make_dense_mpc_hessian()), make_dense_mpc_gradient(),
            Eigen::MatrixXd(make_dense_mpc_constraint_matrix()), make_dense_mpc_lower_bound(),
            make_dense_mpc_upper_bound()};
}

}

#endif
