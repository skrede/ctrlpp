// Usage: ./ctrlpp_advanced_pid_console | gnuplot -p -e "plot '-' using 1:3 with lines"
// Redirect: ./ctrlpp_advanced_pid_console > output.csv

#include "ctrlpp/pid.h"
#include "ctrlpp/pid_policies.h"
#include "ctrlpp/eigen_linalg.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Policy = ctrlpp::EigenLinalgPolicy;
    using Vec = Eigen::Matrix<double, 1, 1>;

    using OuterPid = ctrlpp::Pid<double, 1, 1, 1, Policy,
        ctrlpp::AntiWindup<ctrlpp::BackCalc>,
        ctrlpp::PerfAssessment<ctrlpp::IAE>>;

    OuterPid::config_type outer_cfg{};
    outer_cfg.kp = Vec::Constant(1.5);
    outer_cfg.ki = Vec::Constant(0.3);
    outer_cfg.kd = Vec::Constant(0.0);
    outer_cfg.output_min = Vec::Constant(0.0);
    outer_cfg.output_max = Vec::Constant(100.0);

    OuterPid outer(outer_cfg);

    using InnerPid = ctrlpp::Pid<double, 1, 1, 1, Policy,
        ctrlpp::DerivFilter>;

    InnerPid::config_type inner_cfg{};
    inner_cfg.kp = Vec::Constant(3.0);
    inner_cfg.ki = Vec::Constant(1.0);
    inner_cfg.kd = Vec::Constant(0.2);
    inner_cfg.output_min = Vec::Constant(0.0);
    inner_cfg.output_max = Vec::Constant(100.0);
    inner_cfg.template policy<ctrlpp::DerivFilter>().n = {10.0};

    InnerPid inner(inner_cfg);

    double temperature = 20.0;
    double heater_power = 0.0;
    constexpr double temp_setpoint = 60.0;
    constexpr double dt = 0.01;
    constexpr double duration = 30.0;

    std::cout << "time,setpoint,measurement,control" << std::endl;

    for (double t = 0.0; t < duration; t += dt) {
        Vec sp_outer = Vec::Constant(temp_setpoint);
        Vec meas_outer = Vec::Constant(temperature);
        Vec inner_output = Vec::Constant(heater_power);

        auto heater_sp_vec = outer.compute(sp_outer, meas_outer, dt, inner_output);
        double heater_sp = heater_sp_vec[0];

        Vec sp_inner = Vec::Constant(heater_sp);
        Vec meas_inner = Vec::Constant(heater_power);
        auto control = inner.compute(sp_inner, meas_inner, dt);

        heater_power = 0.8 * heater_power + 0.2 * control[0];
        temperature = 0.99 * temperature + 0.01 * heater_power;

        std::cout << std::fixed << std::setprecision(4)
                  << t << "," << temp_setpoint << "," << temperature << "," << control[0] << std::endl;
    }
}
