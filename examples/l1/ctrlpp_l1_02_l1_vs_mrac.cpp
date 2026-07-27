// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_l1_02_l1_vs_mrac' using 1:2 with lines title 'reference', '' using 1:3 with lines title 'L1 plant', '' using 1:4 with lines title 'MRAC plant'"
// Redirect: ./ctrlpp_l1_02_l1_vs_mrac > output.csv

/// @file ctrlpp_l1_02_l1_vs_mrac.cpp
/// @brief L1 vs MRAC comparison on identical plant.
///
/// Both controllers use the same unknown plant (x=0.8x+0.5u) and
/// reference/predictor model (pole 0.9). Large adaptation gains
/// show that L1's transient is bounded by the filter bandwidth
/// while MRAC's transient depends on the gain magnitude.
///
/// Output: CSV to stdout (pipe to gnuplot or redirect to file)

#include "ctrlpp/control/l1.h"
#include "ctrlpp/control/mrac.h"

#include <iomanip>
#include <iostream>

int main()
{
    // Shared reference/predictor model: pole at 0.9
    ctrlpp::discrete_state_space<double, 1, 1, 1> model{};
    model.A(0, 0) = 0.9;
    model.B(0, 0) = 0.1;
    model.C(0, 0) = 1.0;
    model.D(0, 0) = 0.0;

    // L1 config: large adaptation gain (safe due to filter)
    ctrlpp::l1_config<double, 1, 1> l1_cfg{};
    l1_cfg.predictor_model = model;
    l1_cfg.gamma(0, 0) = 5000.0;
    l1_cfg.theta_min[0] = -10.0;
    l1_cfg.theta_max[0] = 10.0;

    // create validates the filter design and the predictor model, reporting a
    // rejection through ctrlpp::expected<l1_controller, l1_error>.
    auto l1_result = ctrlpp::l1_controller<double>::create(l1_cfg, 5.0, 100.0);
    if(!l1_result.has_value())
    {
        std::cerr << "invalid L1 configuration\n";
        return 1;
    }
    auto& l1_ctrl = *l1_result;

    // MRAC config: large adaptation gains (causes oscillation)
    ctrlpp::mrac_config<double, 1, 1> mrac_cfg{};
    mrac_cfg.reference_model = model;
    mrac_cfg.gamma_x(0, 0) = 5.0;
    mrac_cfg.gamma_r(0, 0) = 5.0;

    ctrlpp::mrac_controller<double> mrac_ctrl(mrac_cfg);

    constexpr double a_p = 0.8;
    constexpr double b_p = 0.5;
    constexpr int steps = 500;
    constexpr double r_val = 1.0;

    double x_l1 = 0.0;
    double x_mrac = 0.0;

    std::cout << "# time,reference,l1_plant,mrac_plant,l1_control,mrac_control,"
                 "l1_sigma_hat,mrac_theta_x\n";

    ctrlpp::Vector<double, 1> r;
    r[0] = r_val;

    for(int k = 0; k < steps; ++k)
    {
        ctrlpp::Vector<double, 1> x_l1_vec;
        x_l1_vec[0] = x_l1;

        ctrlpp::Vector<double, 1> x_mrac_vec;
        x_mrac_vec[0] = x_mrac;

        auto u_l1 = l1_ctrl.evaluate(x_l1_vec, r);
        auto u_mrac = mrac_ctrl.evaluate(x_mrac_vec, r);

        std::cout << std::fixed << std::setprecision(6)
                  << k << ","
                  << r_val << ","
                  << x_l1 << ","
                  << x_mrac << ","
                  << u_l1[0] << ","
                  << u_mrac[0] << ","
                  << l1_ctrl.sigma_hat()[0] << ","
                  << mrac_ctrl.theta_x()(0, 0) << "\n";

        x_l1 = a_p * x_l1 + b_p * u_l1[0];
        x_mrac = a_p * x_mrac + b_p * u_mrac[0];
    }
}
