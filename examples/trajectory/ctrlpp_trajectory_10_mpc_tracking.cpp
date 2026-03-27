// Usage: ./ctrlpp_trajectory_10_mpc_tracking | gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '-' using 1:2 with lines title 'actual', '' using 1:4 with lines title 'reference'"
// Redirect: ./ctrlpp_trajectory_10_mpc_tracking > output.csv

#include <ctrlpp/trajectory/trapezoidal_trajectory.h>
#include <ctrlpp/mpc.h>
#include <ctrlpp/mpc/osqp_solver.h>
#include <ctrlpp/model/propagate.h>

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

    // Double-integrator: x = [position, velocity], u = force
    ctrlpp::discrete_state_space<double, NX, NU, NX> sys{
        .A = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished(),
        .B = (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished(),
        .C = Eigen::Matrix2d::Identity(),
        .D = Eigen::Matrix<double, 2, 1>::Zero()};

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon = horizon,
        .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
        .R = Eigen::Matrix<double, 1, 1>::Constant(0.1),
        .Qf = std::nullopt,
        .u_min = Eigen::Matrix<double, 1, 1>::Constant(-2.0),
        .u_max = Eigen::Matrix<double, 1, 1>::Constant(2.0)};

    ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver> controller(sys, cfg);

    // Trapezoidal trajectory as reference generator (replaces manual ramp)
    ctrlpp::trapezoidal_trajectory<double> traj({.q0 = 0.0, .q1 = 5.0, .v_max = 2.0, .a_max = 4.0});

    // Pre-compute reference trajectory for MPC lookahead
    auto const total_refs = static_cast<std::size_t>(sim_steps + horizon + 1);
    std::vector<Eigen::Vector2d> trajectory(total_refs);

    auto const traj_dur = traj.duration();
    for (std::size_t k = 0; k < total_refs; ++k)
    {
        double const tk = static_cast<double>(k) * dt;
        // Clamp to trajectory duration for times beyond profile end
        double const t_eval = std::min(tk, traj_dur);
        auto const pt = traj.evaluate(t_eval);
        trajectory[k] = Eigen::Vector2d(pt.position[0], pt.velocity[0]);
    }

    Eigen::Vector2d x = Eigen::Vector2d::Zero();

    std::cout << "time,x_pos,x_vel,ref_pos,control\n";

    for (int k = 0; k < sim_steps; ++k)
    {
        auto const ref_start = static_cast<std::size_t>(k);
        std::span<const Eigen::Vector2d> ref_span(
            trajectory.data() + ref_start,
            static_cast<std::size_t>(horizon + 1));

        auto u_opt = controller.solve(x, ref_span);
        if (!u_opt)
        {
            std::cerr << "MPC solve failed at step " << k << "\n";
            return EXIT_FAILURE;
        }

        Eigen::Matrix<double, 1, 1> u = *u_opt;
        double const t = static_cast<double>(k) * dt;

        std::cout << std::fixed << std::setprecision(4) << t << "," << x[0] << "," << x[1] << "," << trajectory[ref_start][0] << "," << u[0] << "\n";

        x = ctrlpp::propagate(sys, x, u);
    }
}
