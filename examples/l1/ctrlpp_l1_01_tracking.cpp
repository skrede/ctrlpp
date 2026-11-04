// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_l1_01_tracking' using 1:2 with lines title 'reference', '' using 1:3 with lines title 'plant', '' using 1:4 with lines title 'predictor'"
// Redirect: ./ctrlpp_l1_01_tracking > output.csv

#include "ctrlpp/control/l1.h"

#include <cstdint>
#include <iomanip>
#include <iostream>

int main()
{
    using controller = ctrlpp::l1_controller<double>;
    using config = controller::config_type;

    ctrlpp::discrete_state_space<double, 1, 1, 1> pred_model{};
    pred_model.A(0, 0) = 0.9;
    pred_model.B(0, 0) = 0.1;
    pred_model.C(0, 0) = 1.0;
    pred_model.D(0, 0) = 0.0;

    config cfg{};
    cfg.predictor_model = pred_model;
    cfg.gamma(0, 0) = 1000.0;
    cfg.theta_min[0] = -10.0;
    cfg.theta_max[0] = 10.0;

    // A second-order Butterworth low-pass at 5 Hz for a 100 Hz sample rate.
    auto ctrl_result = controller::create(cfg, 5.0, 100.0);
    if(!ctrl_result.has_value())
    {
        std::cerr << "invalid L1 configuration\n";
        return 1;
    }
    auto& ctrl = *ctrl_result;

    constexpr double a_p = 0.8;
    constexpr double b_p = 0.5;
    constexpr std::int32_t steps = 500;
    constexpr double r_val = 1.0;

    double x_plant = 0.0;

    std::cout << "# time,reference,plant_state,predictor_state,control,sigma_hat\n";

    ctrlpp::Vector<double, 1> r;
    r[0] = r_val;

    for(std::int32_t k = 0; k < steps; ++k)
    {
        ctrlpp::Vector<double, 1> x;
        x[0] = x_plant;

        // A refused cycle leaves the adaptive parameters untouched, since one non-finite sample would destroy
        // them for good. This example stops; a real caller decides whether the actuator holds the last command,
        // drives a safe value or fails over.
        auto step = ctrl.evaluate(x, r);
        if(!step.has_value())
        {
            std::cerr << "l1 refused the cycle at step " << k << "\n";
            return 1;
        }
        const auto& u = *step;

        std::cout << std::fixed << std::setprecision(6) << k << "," << r_val << "," << x_plant << "," << ctrl.x_hat()[0] << "," << u[0] << "," << ctrl.sigma_hat()[0] << "\n";

        x_plant = a_p * x_plant + b_p * u[0];
    }
}
