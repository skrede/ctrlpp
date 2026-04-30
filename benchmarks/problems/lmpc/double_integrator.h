#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_LMPC_DOUBLE_INTEGRATOR_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_LMPC_DOUBLE_INTEGRATOR_H

#include "ctrlpp/mpc.h"
#include "ctrlpp/model/state_space.h"

#include <Eigen/Core>

#include <cstddef>

namespace ctrlpp::bench::problems::lmpc
{

inline constexpr double double_integrator_dt = 0.1;

inline auto make_double_integrator_4_2_state_space()
    -> ctrlpp::discrete_state_space<double, 4, 2, 4>
{
    Eigen::Matrix4d A = Eigen::Matrix4d::Identity();
    A(0, 1) = double_integrator_dt;
    A(2, 3) = double_integrator_dt;

    Eigen::Matrix<double, 4, 2> B = Eigen::Matrix<double, 4, 2>::Zero();
    B(0, 0) = 0.5 * double_integrator_dt * double_integrator_dt;
    B(1, 0) = double_integrator_dt;
    B(2, 1) = 0.5 * double_integrator_dt * double_integrator_dt;
    B(3, 1) = double_integrator_dt;

    return {
        .A = A,
        .B = B,
        .C = Eigen::Matrix4d::Identity(),
        .D = Eigen::Matrix<double, 4, 2>::Zero(),
    };
}

inline auto make_double_integrator_4_2_config(int horizon)
    -> ctrlpp::mpc_config<double, 4, 2>
{
    ctrlpp::mpc_config<double, 4, 2> cfg{};
    cfg.horizon = horizon;
    cfg.Q = Eigen::Matrix4d::Identity();
    cfg.R = 0.1 * Eigen::Matrix2d::Identity();
    cfg.u_min = Eigen::Vector2d::Constant(-1.0);
    cfg.u_max = Eigen::Vector2d::Constant(1.0);
    cfg.x_min = Eigen::Vector4d::Constant(-5.0);
    cfg.x_max = Eigen::Vector4d::Constant(5.0);
    return cfg;
}

inline auto double_integrator_4_2_x0_default() -> Eigen::Vector4d
{
    Eigen::Vector4d x0;
    x0 << 1.0, 0.0, -0.5, 0.0;
    return x0;
}

}

#endif
