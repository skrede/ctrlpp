// Usage: ./ctrlpp_mpc_01_regulation | gnuplot -p -e "plot '-' using 1:2 with lines"
// Redirect: ./ctrlpp_mpc_01_regulation > output.csv

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/propagate.h"

#include <Eigen/Dense>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr double dt = 0.1;

    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        .A = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        .B = (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(),
        .C = Eigen::Matrix2d::Identity(),
        .D = Eigen::Matrix<double, 2, 1>::Zero()
    };

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 20,
        .Q       = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R       = Eigen::Matrix<double, 1, 1>::Identity(),
        .Qf      = std::nullopt,
        .u_min   = Eigen::Matrix<double, 1, 1>::Constant(-1.0),
        .u_max   = Eigen::Matrix<double, 1, 1>::Constant(1.0),
        .x_min   = Eigen::Vector2d(-std::numeric_limits<double>::infinity(), -2.0),
        .x_max   = Eigen::Vector2d(std::numeric_limits<double>::infinity(), 2.0)
    };

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> controller(sys, cfg);

    Eigen::Vector2d x(5.0, 0.0);
    constexpr double duration = 10.0;

    std::cout << "time,x_0,x_1,control,cost\n";

    for(double t = 0.0; t < duration; t += dt)
    {
        auto u_opt = controller.solve(x);
        if(!u_opt)
        {
            std::cerr << "MPC solve failed at t=" << t << "\n";
            return EXIT_FAILURE;
        }

        auto diag = controller.diagnostics();
        Eigen::Matrix<double, 1, 1> u = *u_opt;

        std::cout << std::fixed << std::setprecision(4)
            << t << "," << x[0] << "," << x[1] << ","
            << u[0] << "," << diag.cost << "\n";

        x = ctrlpp::propagate(sys, x, u);
    }
}
