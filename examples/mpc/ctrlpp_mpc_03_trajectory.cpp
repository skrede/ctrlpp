// Usage: ./ctrlpp_mpc_03_trajectory | gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '-' using 1:2 with lines"
// Redirect: ./ctrlpp_mpc_03_trajectory > output.csv

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/model/propagate.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <span>
#include <vector>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr double dt = 0.1;
    constexpr int horizon = 20;
    constexpr double duration = 10.0;
    constexpr int sim_steps = static_cast<int>(duration / dt);

    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        .A = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(), .B = (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(), .C = Eigen::Matrix2d::Identity(), .D = Eigen::Matrix<double, 2, 1>::Zero()};

    ctrlpp::mpc_config<double, NX, NU> cfg{.horizon = horizon,
                                           .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
                                           .R = Eigen::Matrix<double, 1, 1>::Constant(0.1),
                                           .Qf = std::nullopt,
                                           .u_min = Eigen::Matrix<double, 1, 1>::Constant(-2.0),
                                           .u_max = Eigen::Matrix<double, 1, 1>::Constant(2.0)};

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> controller(sys, cfg);

    // Pre-compute reference trajectory: ramp from 0 to 5 over 5s, then hold
    constexpr double ramp_duration = 5.0;
    constexpr double target_pos = 5.0;
    int total_refs = sim_steps + horizon + 1;
    std::vector<Eigen::Vector2d> trajectory(static_cast<std::size_t>(total_refs));

    for(int k = 0; k < total_refs; ++k)
    {
        double tk = k * dt;
        double pos = (tk < ramp_duration) ? target_pos * tk / ramp_duration : target_pos;
        double vel = (tk < ramp_duration) ? target_pos / ramp_duration : 0.0;
        trajectory[static_cast<std::size_t>(k)] = Eigen::Vector2d(pos, vel);
    }

    Eigen::Vector2d x = Eigen::Vector2d::Zero();

    std::cout << "time,x_0,x_1,x_ref_0,control\n";

    for(int k = 0; k < sim_steps; ++k)
    {
        auto ref_start = static_cast<std::size_t>(k);
        std::span<const Eigen::Vector2d> ref_span(trajectory.data() + ref_start, static_cast<std::size_t>(horizon + 1));

        auto u_opt = controller.solve(x, ref_span);
        if(!u_opt)
        {
            std::cerr << "MPC solve failed at step " << k << "\n";
            return EXIT_FAILURE;
        }

        Eigen::Matrix<double, 1, 1> u = *u_opt;
        double t = k * dt;

        std::cout << std::fixed << std::setprecision(4) << t << "," << x[0] << "," << x[1] << "," << trajectory[ref_start][0] << "," << u[0] << "\n";

        x = ctrlpp::propagate(sys, x, u);
    }
}
