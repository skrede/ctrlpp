// Usage: ./ctrlpp_nmpc_01_pendulum | gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '-' using 1:2 with lines"
// Redirect: ./ctrlpp_nmpc_01_pendulum > output.csv

#include "ctrlpp/mpc/nlopt_solver.h"
#include "ctrlpp/mpc/nmpc_config.h"
#include "ctrlpp/nmpc.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numbers>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr double g = 9.81;
    constexpr double l = 0.5;
    constexpr double m = 0.5;
    constexpr double b = 0.1;
    constexpr double dt = 0.05;

    auto dynamics = [](const ctrlpp::Vector<double, NX>& x, const ctrlpp::Vector<double, NU>& u) -> ctrlpp::Vector<double, NX>
    {
        double theta = x[0];
        double theta_dot = x[1];
        double torque = u[0];
        double theta_ddot = (g / l) * std::sin(theta) - b * theta_dot + torque / (m * l * l);
        return (ctrlpp::Vector<double, NX>() << theta + theta_dot * dt, theta_dot + theta_ddot * dt).finished().eval();
    };

    ctrlpp::nmpc_config<double, NX, NU> cfg{.horizon = 30,
                                            .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
                                            .R = Eigen::Matrix<double, 1, 1>::Constant(0.01),
                                            .Qf = Eigen::Vector2d(100.0, 10.0).asDiagonal(),
                                            .u_min = Eigen::Matrix<double, 1, 1>::Constant(-5.0),
                                            .u_max = Eigen::Matrix<double, 1, 1>::Constant(5.0)};

    ctrlpp::nmpc<double, NX, NU, ctrlpp::nlopt_solver<double>, decltype(dynamics)> controller(dynamics, cfg);

    Eigen::Vector2d x(std::numbers::pi - 0.3, 0.0);
    Eigen::Vector2d x_ref(std::numbers::pi, 0.0);
    constexpr double duration = 10.0;

    std::cout << "time,theta,theta_dot,torque\n";

    for(double t = 0.0; t < duration; t += dt)
    {
        auto u_opt = controller.solve(x, x_ref);
        if(!u_opt)
        {
            std::cerr << "NMPC solve failed at t=" << t << "\n";
            return EXIT_FAILURE;
        }

        Eigen::Matrix<double, 1, 1> u = *u_opt;

        std::cout << std::fixed << std::setprecision(4) << t << "," << x[0] << "," << x[1] << "," << u[0] << "\n";

        x = dynamics(x, u);
    }
}
