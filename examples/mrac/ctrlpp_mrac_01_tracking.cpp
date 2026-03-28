// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_mrac_01_tracking' using 1:2 with lines title 'reference', '' using 1:3 with lines title 'plant', '' using 1:4 with lines title 'model'"
// Redirect: ./ctrlpp_mrac_01_tracking > output.csv

/// @file ctrlpp_mrac_01_tracking.cpp
/// @brief MRAC step tracking with Lyapunov adaptation.
///
/// Demonstrates SISO MRAC tracking a step reference.
/// Unknown plant: x_{k+1} = 0.8 x_k + 0.5 u_k
/// Reference model: x_{k+1} = 0.9 x_k + 0.1 r_k
///
/// Output: CSV to stdout (pipe to gnuplot or redirect to file)

#include "ctrlpp/control/mrac.h"

#include <iomanip>
#include <iostream>

int main()
{
    using controller = ctrlpp::mrac_controller<double>;
    using config = controller::config_type;

    ctrlpp::discrete_state_space<double, 1, 1, 1> ref_model{};
    ref_model.A(0, 0) = 0.9;
    ref_model.B(0, 0) = 0.1;
    ref_model.C(0, 0) = 1.0;
    ref_model.D(0, 0) = 0.0;

    config cfg{};
    cfg.reference_model = ref_model;
    cfg.gamma_x << 0.5;
    cfg.gamma_r << 0.5;

    controller ctrl(cfg);

    constexpr double a_p = 0.8;
    constexpr double b_p = 0.5;
    constexpr int steps = 500;
    constexpr double r_val = 1.0;

    double x_plant = 0.0;

    std::cout << "# time,reference,plant_state,model_state,control,theta_x,theta_r\n";

    ctrlpp::Vector<double, 1> r;
    r[0] = r_val;

    for(int k = 0; k < steps; ++k)
    {
        ctrlpp::Vector<double, 1> x;
        x[0] = x_plant;

        auto u = ctrl.evaluate(x, r);

        std::cout << std::fixed << std::setprecision(6)
                  << k << ","
                  << r_val << ","
                  << x_plant << ","
                  << ctrl.x_model()[0] << ","
                  << u[0] << ","
                  << ctrl.theta_x()(0, 0) << ","
                  << ctrl.theta_r()(0, 0) << "\n";

        x_plant = a_p * x_plant + b_p * u[0];
    }
}
