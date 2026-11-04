// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_mpc_02_observer_mpc' using 1:2 with lines title 'true x0', '' using 1:4 with lines title 'est x0', '' using 1:6 with lines title 'control'"
// Redirect: ./ctrlpp_mpc_02_observer_mpc > output.csv

#include "ctrlpp/mpc.h"

#include "ctrlpp/estimation/kalman.h"

#include "ctrlpp/model/propagate.h"

#include "ctrlpp/mpc/osqp_solver.h"

#include <Eigen/Dense>

#include <limits>
#include <cstdlib>
#include <iomanip>
#include <iostream>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 1;
    constexpr std::size_t NY_OBS = 1;
    constexpr double dt = 0.1;

    Eigen::Matrix2d Ad = (Eigen::Matrix2d() << 1.0, dt, 0.0, 1.0).finished();
    Eigen::Vector2d Bd = (Eigen::Vector2d() << 0.5 * dt * dt, dt).finished();

    // The observer measures position only, while the MPC needs the full state as its output.
    ctrlpp::discrete_state_space<double, NX, NU, NY_OBS> obs_sys{Ad, Bd, (Eigen::Matrix<double, 1, 2>() << 1.0, 0.0).finished(), Eigen::Matrix<double, 1, 1>::Zero()};

    ctrlpp::discrete_state_space<double, NX, NU, NX> mpc_sys{Ad, Bd, Eigen::Matrix2d::Identity(), Eigen::Matrix<double, 2, 1>::Zero()};

    auto kf_result = ctrlpp::kalman_filter<double, NX, NU, NY_OBS>::create(
            obs_sys, {.Q = Eigen::Matrix2d::Identity() * 0.01, .R = Eigen::Matrix<double, 1, 1>::Constant(0.1), .x0 = Eigen::Vector2d::Zero(), .P0 = Eigen::Matrix2d::Identity()});
    if(!kf_result.has_value())
    {
        std::cerr << "invalid Kalman filter configuration\n";
        return 1;
    }
    auto& kf = *kf_result;

    ctrlpp::mpc_config<double, NX, NU> cfg{.horizon = 20,
            .Q = Eigen::Vector2d(10.0, 1.0).asDiagonal(),
            .R = Eigen::Matrix<double, 1, 1>::Identity(),
            .Qf = std::nullopt,
            .u_min = Eigen::Matrix<double, 1, 1>::Constant(-1.0),
            .u_max = Eigen::Matrix<double, 1, 1>::Constant(1.0),
            .x_min = Eigen::Vector2d(-std::numeric_limits<double>::infinity(), -2.0),
            .x_max = Eigen::Vector2d(std::numeric_limits<double>::infinity(), 2.0)};

    auto controller_result = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(mpc_sys, cfg);
    if(!controller_result.has_value())
    {
        std::cerr << "invalid MPC configuration\n";
        return EXIT_FAILURE;
    }
    auto& controller = *controller_result;

    Eigen::Vector2d x_true(5.0, 0.0);
    Eigen::Matrix<double, 1, 1> u = Eigen::Matrix<double, 1, 1>::Zero();
    constexpr double duration = 10.0;

    std::cout << "time,x_true_0,x_true_1,x_est_0,x_est_1,control\n";

    for(double t = 0.0; t < duration; t += dt)
    {
        kf.predict(u);

        x_true = ctrlpp::propagate(mpc_sys, x_true, u);

        Eigen::Matrix<double, 1, 1> z = obs_sys.C * x_true;

        // update reports whether it ran, so a non-finite sample cannot silently poison the estimate the MPC
        // then solves against.
        if(const auto stepped = kf.update(z); !stepped)
        {
            std::cerr << "Kalman filter rejected the measurement at t=" << t << "\n";
            return EXIT_FAILURE;
        }

        auto u_opt = controller.solve(kf.state());
        if(!u_opt)
        {
            std::cerr << "MPC solve failed at t=" << t << "\n";
            return EXIT_FAILURE;
        }
        u = u_opt->input;

        std::cout << std::fixed << std::setprecision(4) << t << "," << x_true[0] << "," << x_true[1] << "," << kf.state()[0] << "," << kf.state()[1] << "," << u[0] << "\n";
    }
}
