// Usage: ./ctrlpp_mpc_04_mobile_robot | gnuplot -p -e "plot '-' using 2:3 with lines"
// Redirect: ./ctrlpp_mpc_04_mobile_robot > output.csv

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/propagate.h"

#include <Eigen/Dense>

#include <cstdlib>
#include <iomanip>
#include <iostream>

int main()
{
    constexpr std::size_t NX = 3;
    constexpr std::size_t NU = 2;
    constexpr double dt = 0.1;
    constexpr double v0 = 1.0;

    // Linearized unicycle about straight-line motion (theta=0, v=v0)
    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        .A = (Eigen::Matrix3d() <<
            1.0, 0.0, 0.0,
            0.0, 1.0, v0 * dt,
            0.0, 0.0, 1.0).finished(),
        .B = (Eigen::Matrix<double, 3, 2>() <<
            dt, 0.0,
            0.0, 0.0,
            0.0, dt).finished(),
        .C = Eigen::Matrix3d::Identity(),
        .D = Eigen::Matrix<double, 3, 2>::Zero()
    };

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = 15,
        .Q = Eigen::Vector3d(1.0, 10.0, 5.0).asDiagonal(),
        .R = Eigen::Vector2d(0.1, 0.1).asDiagonal(),
        .Qf = std::nullopt,
        .u_min = Eigen::Vector2d(0.0, -1.0),
        .u_max = Eigen::Vector2d(2.0, 1.0),
        .x_min = std::nullopt,
        .x_max = std::nullopt,
        .du_max = std::nullopt,
        .soft_penalty = 1e4,
        .hard_state_constraints = false
    };

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> controller(sys, cfg);

    // Initial state: offset from straight-line path
    Eigen::Vector3d x(0.0, 0.5, 0.1);
    constexpr double duration = 10.0;
    constexpr int sim_steps = static_cast<int>(duration / dt);

    std::cout << "time,x,y,theta,v,omega\n";

    for (int k = 0; k < sim_steps; ++k) {
        // Reference: straight line at nominal speed, y=0, theta=0
        Eigen::Vector3d x_ref(v0 * k * dt, 0.0, 0.0);

        auto u_opt = controller.solve(x, x_ref);
        if (!u_opt) {
            std::cerr << "MPC solve failed at step " << k << "\n";
            return EXIT_FAILURE;
        }

        Eigen::Vector2d u = *u_opt;
        double t = k * dt;

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << x[0] << "," << x[1] << "," << x[2] << ","
                  << u[0] << "," << u[1] << "\n";

        x = ctrlpp::propagate(sys, x, u);
    }
}
