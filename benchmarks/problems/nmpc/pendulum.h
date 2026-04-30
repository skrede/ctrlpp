#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_NMPC_PENDULUM_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_NMPC_PENDULUM_H

#include "ctrlpp/mpc/nmpc_config.h"

#include <Eigen/Core>

#include <cmath>
#include <cstddef>

namespace ctrlpp::bench::problems::nmpc
{

inline constexpr double pendulum_dt = 0.05;
inline constexpr double pendulum_gravity = 9.81;
inline constexpr double pendulum_length = 1.0;

inline auto pendulum_2(const Eigen::Vector2d& x,
                       const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    const double theta = x(0);
    const double omega = x(1);
    const double alpha = -pendulum_gravity / pendulum_length * std::sin(theta) + u(0);
    return Eigen::Vector2d{theta + pendulum_dt * omega,
                           omega + pendulum_dt * alpha};
}

inline auto make_pendulum_config(int horizon) -> ctrlpp::nmpc_config<double, 2, 1>
{
    return {
        .horizon = horizon,
        .Q = Eigen::Matrix2d::Identity(),
        .R = Eigen::Matrix<double, 1, 1>::Identity() * 0.1,
    };
}

inline auto pendulum_x0_default() -> Eigen::Vector2d
{
    return {0.6, 0.0};
}

}

#endif
