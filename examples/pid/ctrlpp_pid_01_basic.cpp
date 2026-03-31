// Usage: gnuplot -p -e "set datafile separator ','; set key autotitle columnheader; plot '<./ctrlpp_pid_01_basic' using 1:2 with lines title 'setpoint', '' using 1:3 with lines title 'measurement', '' using 1:4 with lines title 'control'"
// Redirect: ./ctrlpp_pid_01_basic > output.csv

#include "ctrlpp/control/pid.h"

#include <iomanip>
#include <iostream>

int main()
{
    using Pid = ctrlpp::pid<double, 1, 1, 1>;
    using Vec = Pid::vector_t;

    Pid::config_type cfg{};
    cfg.kp = Vec::Constant(2.0);
    cfg.ki = Vec::Constant(1.0);
    cfg.kd = Vec::Constant(0.0);
    cfg.output_min = Vec::Constant(-10.0);
    cfg.output_max = Vec::Constant(10.0);

    Pid ctrl(cfg);

    double y = 0.0;
    constexpr double setpoint = 1.0;
    constexpr double dt = 0.01;
    constexpr double a = 0.9;
    constexpr double duration = 10.0;

    std::cout << "time,setpoint,measurement,control\n";

    for(double t = 0.0; t < duration; t += dt)
    {
        auto sp = Vec::Constant(setpoint);
        auto meas = Vec::Constant(y);
        auto u = ctrl.compute(sp, meas, dt);
        y = a * y + (1.0 - a) * u[0];

        std::cout << std::fixed << std::setprecision(4) << t << "," << setpoint << "," << y << "," << u[0] << "\n";
    }
}
