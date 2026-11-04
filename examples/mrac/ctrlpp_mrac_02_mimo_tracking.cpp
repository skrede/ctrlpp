// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_mrac_02_mimo_tracking' using 1:2 with lines title 'ref_ch1', '' using 1:3 with lines title 'plant_ch1', '' using 1:4 with lines title 'model_ch1', '' using 1:5 with lines title 'ref_ch2', '' using 1:6 with lines title 'plant_ch2', '' using 1:7 with lines title 'model_ch2'"
// Redirect: ./ctrlpp_mrac_02_mimo_tracking > output.csv

#include "ctrlpp/control/mrac.h"

#include <cstdint>
#include <iomanip>
#include <iostream>

int main()
{
    constexpr std::size_t NX = 2;
    constexpr std::size_t NU = 2;

    using controller = ctrlpp::mrac_controller<double, NX, NU>;
    using config = controller::config_type;

    ctrlpp::discrete_state_space<double, NX, NU, NX> ref_model{};
    ref_model.A << 0.9, 0.0, 0.0, 0.85;
    ref_model.B << 0.1, 0.0, 0.0, 0.15;
    ref_model.C = ctrlpp::Matrix<double, NX, NX>::Identity();
    ref_model.D = ctrlpp::Matrix<double, NX, NU>::Zero();

    config cfg{};
    cfg.reference_model = ref_model;
    cfg.gamma_x = 0.3 * ctrlpp::Matrix<double, NX, NX>::Identity();
    cfg.gamma_r = 0.3 * ctrlpp::Matrix<double, NU, NU>::Identity();

    controller ctrl(cfg);

    ctrlpp::Matrix<double, NX, NX> A_p;
    A_p << 0.8, 0.1, 0.0, 0.7;
    ctrlpp::Matrix<double, NX, NU> B_p;
    B_p << 0.5, 0.0, 0.0, 0.4;

    ctrlpp::Vector<double, NX> x_plant = ctrlpp::Vector<double, NX>::Zero();
    ctrlpp::Vector<double, NU> r;
    r << 1.0, 0.5;

    constexpr std::int32_t steps = 500;

    std::cout << "# time,ref_ch1,plant_ch1,model_ch1,ref_ch2,plant_ch2,model_ch2,u_ch1,u_ch2\n";

    for(std::int32_t k = 0; k < steps; ++k)
    {
        // A refused cycle leaves the adaptive parameters untouched, since one non-finite sample would destroy
        // them for good. This example stops; a real caller decides whether the actuator holds the last command,
        // drives a safe value or fails over.
        auto step = ctrl.evaluate(x_plant, r);
        if(!step.has_value())
        {
            std::cerr << "mrac refused the cycle at step " << k << "\n";
            return 1;
        }
        const auto& u = *step;

        std::cout << std::fixed << std::setprecision(6) << k << "," << r[0] << "," << x_plant[0] << "," << ctrl.x_model()[0] << "," << r[1] << "," << x_plant[1] << ","
                  << ctrl.x_model()[1] << "," << u[0] << "," << u[1] << "\n";

        x_plant = A_p * x_plant + B_p * u;
    }
}
