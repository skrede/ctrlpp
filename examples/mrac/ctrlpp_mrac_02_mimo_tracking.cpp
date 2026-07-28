// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_mrac_02_mimo_tracking' using 1:2 with lines title 'ref_ch1', '' using 1:3 with lines title 'plant_ch1', '' using 1:4 with lines title 'model_ch1', '' using 1:5 with lines title 'ref_ch2', '' using 1:6 with lines title 'plant_ch2', '' using 1:7 with lines title 'model_ch2'"
// Redirect: ./ctrlpp_mrac_02_mimo_tracking > output.csv

/// @file ctrlpp_mrac_02_mimo_tracking.cpp
/// @brief MIMO MRAC step tracking with Lyapunov adaptation.
///
/// Demonstrates 2-channel MRAC tracking with matrix adaptation gains.
/// Unknown plant: x_{k+1} = A_p * x_k + B_p * u_k
///   A_p = [[0.8, 0.1], [0.0, 0.7]], B_p = [[0.5, 0.0], [0.0, 0.4]]
/// Reference model: A_m = [[0.9, 0.0], [0.0, 0.85]], B_m = [[0.1, 0.0], [0.0, 0.15]]
///
/// Output: CSV to stdout (pipe to gnuplot or redirect to file)

#include "ctrlpp/control/mrac.h"

#include <iomanip>
#include <iostream>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 2;

    using controller = ctrlpp::mrac_controller<double, NX, NU>;
    using config = controller::config_type;

    ctrlpp::discrete_state_space<double, NX, NU, NX> ref_model{};
    ref_model.A << 0.9, 0.0,
                   0.0, 0.85;
    ref_model.B << 0.1, 0.0,
                   0.0, 0.15;
    ref_model.C = ctrlpp::Matrix<double, NX, NX>::Identity();

    config cfg{};
    cfg.reference_model = ref_model;
    cfg.gamma_x = 0.3 * ctrlpp::Matrix<double, NX, NX>::Identity();
    cfg.gamma_r = 0.3 * ctrlpp::Matrix<double, NU, NU>::Identity();

    controller ctrl(cfg);

    ctrlpp::Matrix<double, NX, NX> A_p;
    A_p << 0.8, 0.1,
           0.0, 0.7;
    ctrlpp::Matrix<double, NX, NU> B_p;
    B_p << 0.5, 0.0,
           0.0, 0.4;

    ctrlpp::Vector<double, NX> x_plant = ctrlpp::Vector<double, NX>::Zero();
    ctrlpp::Vector<double, NU> r;
    r << 1.0, 0.5;

    constexpr int steps = 500;

    std::cout << "# time,ref_ch1,plant_ch1,model_ch1,ref_ch2,plant_ch2,model_ch2,u_ch1,u_ch2\n";

    for(int k = 0; k < steps; ++k)
    {
        // A refused cycle produced no command AND left the adaptive parameters
        // untouched, which is the point of the guard: admitting one non-finite
        // sample would destroy them permanently, because nothing re-derives
        // them. This example stops; a real caller must decide what the actuator
        // does -- hold the last command, drive a configured safe value, or fail
        // over.
        auto step = ctrl.evaluate(x_plant, r);
        if(!step.has_value())
        {
            std::cerr << "mrac refused the cycle at step " << k << "\n";
            return 1;
        }
        const auto& u = *step;

        std::cout << std::fixed << std::setprecision(6)
                  << k << ","
                  << r[0] << "," << x_plant[0] << "," << ctrl.x_model()[0] << ","
                  << r[1] << "," << x_plant[1] << "," << ctrl.x_model()[1] << ","
                  << u[0] << "," << u[1] << "\n";

        x_plant = A_p * x_plant + B_p * u;
    }
}
