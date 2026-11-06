#ifndef HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_REFERENCE_H
#define HPP_GUARD_CTRLPP_EXAMPLES_EMBEDDED_SHARED_PREDICTIVE_REFERENCE_H

#include "predictive_demo.h"
#include "golden_reference.h"
#include "predictive_verdict.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstddef>

namespace ctrlpp
{

struct predictive_constants
{
    double k0;
    double k1;
    double input_row;
    double kkt_condition;
    double jacobian_sigma;
    double jacobian_norm;
    double cost;
};

namespace detail
{

constexpr int predictive_kkt_size = golden::predictive_nv + golden::predictive_max_m;

using kkt_matrix = Eigen::Matrix<double, predictive_kkt_size, predictive_kkt_size>;
using jacobian_matrix = Eigen::Matrix<double, golden::predictive_max_m, golden::predictive_nv>;

inline Eigen::Matrix2d predictive_state_matrix()
{
    Eigen::Matrix2d A;
    A << 1.0, golden::predictive_dt, 0.0, 1.0;
    return A;
}

// Backward Riccati recursion over the horizon, with the terminal weight the formulation substitutes when
// none is configured: the state weight.
inline Eigen::RowVector2d finite_horizon_gain()
{
    const Eigen::Matrix2d A = predictive_state_matrix();
    const Eigen::Vector2d B(0.0, golden::predictive_dt);
    const Eigen::Matrix2d Q = golden::predictive_state_weight * Eigen::Matrix2d::Identity();
    Eigen::Matrix2d P = Q;
    Eigen::RowVector2d K = Eigen::RowVector2d::Zero();
    for(std::size_t k = 0; k < golden::predictive_nh; ++k)
    {
        const double s = golden::predictive_input_weight + B.dot(P * B);
        K = (B.transpose() * P * A) / s;
        P = Q + A.transpose() * P * (A - B * K);
    }
    return K;
}

// Rows: the initial state, then x_{k+1} - A*x_k - B*u_k, in the formulation's layout of every state first
// and every input after.
inline jacobian_matrix predictive_jacobian()
{
    const Eigen::Index nx = static_cast<Eigen::Index>(golden::predictive_nx);
    const Eigen::Index inputs = (static_cast<Eigen::Index>(golden::predictive_nh) + 1) * nx;
    jacobian_matrix J = jacobian_matrix::Zero();
    J.block<2, 2>(0, 0) = Eigen::Matrix2d::Identity();
    for(Eigen::Index k = 0; k < static_cast<Eigen::Index>(golden::predictive_nh); ++k)
    {
        J.block<2, 2>((k + 1) * nx, (k + 1) * nx) = Eigen::Matrix2d::Identity();
        J.block<2, 2>((k + 1) * nx, k * nx) = -predictive_state_matrix();
        J((k + 1) * nx + 1, inputs + k) = -golden::predictive_dt;
    }
    return J;
}

inline kkt_matrix predictive_kkt()
{
    const Eigen::Index inputs = (static_cast<Eigen::Index>(golden::predictive_nh) + 1) * static_cast<Eigen::Index>(golden::predictive_nx);
    kkt_matrix M = kkt_matrix::Zero();
    for(Eigen::Index i = 0; i < golden::predictive_nv; ++i)
        M(i, i) = i < inputs ? golden::predictive_state_weight : golden::predictive_input_weight;
    const jacobian_matrix J = predictive_jacobian();
    M.block<golden::predictive_max_m, golden::predictive_nv>(golden::predictive_nv, 0) = J;
    M.block<golden::predictive_nv, golden::predictive_max_m>(0, golden::predictive_nv) = J.transpose();
    return M;
}

inline double reference_cost(const Eigen::RowVector2d& K)
{
    Eigen::Vector2d x(1.0, 0.0);
    double cost = 0.0;
    for(std::size_t k = 0; k < golden::predictive_run_steps; ++k)
    {
        const double u = -(K(0) * x[0] + K(1) * x[1]);
        cost += stage_cost(x, u);
        x = predictive_plant<double>{}(x, Eigen::Matrix<double, 1, 1>(u));
    }
    return cost;
}

}

// The first-input row of M^-1 is y^T for M^T*y = e_u0; the condition number and the Jacobian's norms are
// what the bounds in predictive_bound.h charge the host computation of these figures at.
inline predictive_constants compute_predictive_constants()
{
    const detail::kkt_matrix M = detail::predictive_kkt();
    const detail::jacobian_matrix J = detail::predictive_jacobian();
    const detail::kkt_matrix inverse = M.fullPivLu().inverse();
    const Eigen::Index u0 = (static_cast<Eigen::Index>(golden::predictive_nh) + 1) * static_cast<Eigen::Index>(golden::predictive_nx);
    const Eigen::RowVector2d K = detail::finite_horizon_gain();
    const auto sigma = J.jacobiSvd().singularValues();
    const double one_norm = J.cwiseAbs().colwise().sum().maxCoeff();
    const double inf_norm = J.cwiseAbs().rowwise().sum().maxCoeff();
    const double condition = M.cwiseAbs().rowwise().sum().maxCoeff() * inverse.cwiseAbs().rowwise().sum().maxCoeff();
    return {K(0), K(1), inverse.row(u0).cwiseAbs().sum(), condition, sigma.minCoeff(), std::sqrt(one_norm * inf_norm), detail::reference_cost(K)};
}

}

#endif
