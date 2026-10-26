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
#include <string_view>

namespace
{

std::string_view name_of(ctrlpp::terminal_set_error e)
{
    switch(e)
    {
        case ctrlpp::terminal_set_error::input_zero_not_interior:
            return "input_zero_not_interior";
        case ctrlpp::terminal_set_error::empty_terminal_set:
            return "empty_terminal_set";
        case ctrlpp::terminal_set_error::dare_failed:
            return "dare_failed";
        case ctrlpp::terminal_set_error::halfplanes_truncated:
            return "halfplanes_truncated";
        case ctrlpp::terminal_set_error::not_converged:
            return "not_converged";
        case ctrlpp::terminal_set_error::closed_loop_not_stable:
            return "closed_loop_not_stable";
        case ctrlpp::terminal_set_error::unbounded_terminal_set:
            return "unbounded_terminal_set";
        case ctrlpp::terminal_set_error::origin_not_interior:
            return "origin_not_interior";
        case ctrlpp::terminal_set_error::dimension_unsupported:
            return "dimension_unsupported";
        case ctrlpp::terminal_set_error::certification_failed:
            return "certification_failed";
        case ctrlpp::terminal_set_error::non_finite_input:
            return "non_finite_input";
    }
    return "unknown";
}

}

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

    // Recursive feasibility and stability need hard state constraints and a terminal set invariant under the LQR law.
    // That is invariant_set; the ellipsoid in ingredients->set is not invariant once encoded in the QP.
    auto ingredients = ctrlpp::terminal_ingredients<double, NX, NU>(sys.A, sys.B, Q, R, -u_max, u_max, -x_max, x_max);
    if(!ingredients)
    {
        std::cerr << "terminal ingredients refused: " << name_of(ingredients.error()) << "\n";
        return EXIT_FAILURE;
    }
    if(!ingredients->invariant_set)
    {
        std::cerr << "invariant terminal set refused: " << name_of(ingredients->invariant_set.error()) << "\n";
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
