// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_mpc_05_terminal_set' using 1:2 with lines title 'x0', '' using 1:3 with lines title 'x1', '' using 1:4 with lines title 'control'"
// Redirect: ./ctrlpp_mpc_05_terminal_set > output.csv

#include "ctrlpp/mpc.h"
#include "ctrlpp/mpc/osqp_solver.h"
#include "ctrlpp/mpc/invariant_set.h"
#include "ctrlpp/model/propagate.h"

#include <Eigen/Dense>

#include <cstdlib>
#include <iomanip>
#include <iostream>

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

    Eigen::Matrix2d const Q = Eigen::Matrix2d::Identity();
    Eigen::Matrix<double, 1, 1> const R = Eigen::Matrix<double, 1, 1>::Identity();
    Eigen::Matrix<double, 1, 1> const u_max = Eigen::Matrix<double, 1, 1>::Constant(1.0);
    Eigen::Vector2d const x_max = Eigen::Vector2d::Constant(10.0);

    // Qf is the Riccati solution, and invariant_set is the maximal output admissible set of the
    // LQR law: the terminal set under which recursive feasibility and stability hold with hard
    // state constraints. The ellipsoid in ingredients->set is not invariant once encoded in the QP.
    auto ingredients = ctrlpp::terminal_ingredients<double, NX, NU>(sys.A, sys.B, Q, R, -u_max, u_max, -x_max, x_max);
    if(!ingredients || !ingredients->invariant_set)
    {
        std::cerr << "terminal ingredients refused\n";
        return EXIT_FAILURE;
    }
    std::cerr << "terminal set rows: " << ingredients->invariant_set->H.rows() << "\n";

    ctrlpp::mpc_config<double, NX, NU> cfg{
        .horizon                 = 30,
        .Q                       = Q,
        .R                       = R,
        .Qf                      = ingredients->Qf,
        .u_min                   = -u_max,
        .u_max                   = u_max,
        .x_min                   = -x_max,
        .x_max                   = x_max,
        .terminal_constraint_set = ctrlpp::terminal_set<double, NX>{*ingredients->invariant_set},
        .hard_state_constraints  = true
    };

    auto controller_result = ctrlpp::mpc<double, NX, NU, ctrlpp::osqp_solver>::create(sys, cfg);
    if(!controller_result.has_value())
    {
        std::cerr << "invalid MPC configuration\n";
        return EXIT_FAILURE;
    }
    auto& controller = *controller_result;

    Eigen::Vector2d x(4.0, 0.0);
    constexpr double duration = 8.0;

    std::cout << "time,x_0,x_1,control\n";

    for(double t = 0.0; t < duration; t += dt)
    {
        auto u_opt = controller.solve(x);
        if(!u_opt)
        {
            std::cerr << "MPC solve failed at t=" << t << "\n";
            return EXIT_FAILURE;
        }

        Eigen::Matrix<double, 1, 1> u = u_opt->input;
        std::cout << std::fixed << std::setprecision(4) << t << "," << x[0] << "," << x[1] << "," << u[0] << "\n";

        x = ctrlpp::propagate(sys, x, u);
    }
}
