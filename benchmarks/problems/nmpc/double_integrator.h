#ifndef HPP_GUARD_BENCHMARKS_PROBLEMS_NMPC_DOUBLE_INTEGRATOR_H
#define HPP_GUARD_BENCHMARKS_PROBLEMS_NMPC_DOUBLE_INTEGRATOR_H

#include "ctrlpp/mpc/nmpc_config.h"

#include <Eigen/Core>

#include <cstddef>

namespace ctrlpp::bench::problems::nmpc
{

inline constexpr double double_integrator_dt = 0.1;

inline auto double_integrator_2(const Eigen::Vector2d& x,
                                const Eigen::Matrix<double, 1, 1>& u) -> Eigen::Vector2d
{
    return Eigen::Vector2d{x(0) + double_integrator_dt * x(1),
                           x(1) + double_integrator_dt * u(0)};
}

inline auto double_integrator_4(const Eigen::Vector4d& x,
                                const Eigen::Vector2d& u) -> Eigen::Vector4d
{
    return Eigen::Vector4d{
        x(0) + double_integrator_dt * x(1),
        x(1) + double_integrator_dt * u(0),
        x(2) + double_integrator_dt * x(3),
        x(3) + double_integrator_dt * u(1)};
}

inline auto double_integrator_8(const Eigen::Matrix<double, 8, 1>& x,
                                const Eigen::Vector4d& u) -> Eigen::Matrix<double, 8, 1>
{
    Eigen::Matrix<double, 8, 1> xn;
    xn(0) = x(0) + double_integrator_dt * x(1);
    xn(1) = x(1) + double_integrator_dt * u(0);
    xn(2) = x(2) + double_integrator_dt * x(3);
    xn(3) = x(3) + double_integrator_dt * u(1);
    xn(4) = x(4) + double_integrator_dt * x(5);
    xn(5) = x(5) + double_integrator_dt * u(2);
    xn(6) = x(6) + double_integrator_dt * x(7);
    xn(7) = x(7) + double_integrator_dt * u(3);
    return xn;
}

/// The same plant, carrying the analytic partials a solver can be handed
/// instead of differencing the map itself.
struct differentiable_double_integrator_2
{
    auto operator()(const Eigen::Vector2d& x, const Eigen::Matrix<double, 1, 1>& u) const -> Eigen::Vector2d
    {
        return double_integrator_2(x, u);
    }

    auto jacobian_x(const Eigen::Vector2d&, const Eigen::Matrix<double, 1, 1>&) const -> Eigen::Matrix2d
    {
        Eigen::Matrix2d partials;
        partials << 1.0, double_integrator_dt, 0.0, 1.0;
        return partials;
    }

    auto jacobian_u(const Eigen::Vector2d&, const Eigen::Matrix<double, 1, 1>&) const
        -> Eigen::Matrix<double, 2, 1>
    {
        return Eigen::Matrix<double, 2, 1>{0.0, double_integrator_dt};
    }
};

template <std::size_t NX, std::size_t NU>
auto make_nmpc_quadratic_config(int horizon) -> ctrlpp::nmpc_config<double, NX, NU>
{
    return {
        .horizon = horizon,
        .Q = Eigen::Matrix<double, NX, NX>::Identity(),
        .R = Eigen::Matrix<double, NU, NU>::Identity() * 0.1,
    };
}

template <std::size_t NX>
auto unit_first_axis_x0() -> Eigen::Matrix<double, NX, 1>
{
    Eigen::Matrix<double, NX, 1> x0 = Eigen::Matrix<double, NX, 1>::Zero();
    x0(0) = 1.0;
    return x0;
}

}

#endif
